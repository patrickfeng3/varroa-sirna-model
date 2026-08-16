"use strict";

const METRICS = new Set(["layer1", "layer2", "layer3", "total"]);
const FEATURE_COLORS = { "5UTR": "rgba(64,145,108,0.10)", "CDS": "rgba(42,111,151,0.08)", "3UTR": "rgba(224,122,95,0.10)" };
let payload;
let scoreIndex;

function safeState(searchParams) {
  const supported = payload.supported_guide_lengths;
  const requestedGuide = Number.parseInt(searchParams.get("guide"), 10);
  const guide = supported.includes(requestedGuide) ? requestedGuide : supported[0];
  const defaultRegion = Math.min(Math.max(96, guide), payload.target.transcript_length_nt);
  const requestedRegion = Number.parseInt(searchParams.get("region"), 10);
  const region = Number.isInteger(requestedRegion) && requestedRegion >= guide && requestedRegion <= payload.target.transcript_length_nt
    ? requestedRegion : defaultRegion;
  const requestedMetric = searchParams.get("metric");
  const metric = METRICS.has(requestedMetric) ? requestedMetric : "total";
  return { guide, region, metric };
}

function buildIndex() {
  const index = {};
  for (const guide of payload.supported_guide_lengths) {
    const rows = payload.candidates.filter(row => row.candidate_length_nt === guide)
      .sort((a, b) => a.target_start_1based - b.target_start_1based);
    index[guide] = { rows, prefix: {} };
    for (const metric of METRICS) {
      const prefix = [0];
      for (const row of rows) prefix.push(prefix[prefix.length - 1] + row[metric]);
      index[guide].prefix[metric] = prefix;
    }
  }
  return index;
}

function renderGuideChoices(supported, selected) {
  const container = document.getElementById("guide-length-options");
  container.innerHTML = supported.map(length =>
    `<label><input type="radio" name="guide-length" value="${length}" ${length === selected ? "checked" : ""}><span>${length} nt</span></label>`
  ).join("");
}

function selectedGuideLength() {
  const selected = document.querySelector('input[name="guide-length"]:checked');
  return selected ? Number.parseInt(selected.value, 10) : payload.supported_guide_lengths[0];
}

function startFeature(position) {
  const annotation = payload.target.annotations.find(item => position >= item.start_1based && position <= item.end_1based);
  return annotation ? annotation.feature : null;
}

function calculateRegions(guide, regionLength, metric) {
  const transcriptLength = payload.target.transcript_length_nt;
  if (!payload.supported_guide_lengths.includes(guide) || !METRICS.has(metric) || regionLength < guide || regionLength > transcriptLength) return [];
  const contained = regionLength - guide + 1;
  const prefix = scoreIndex[guide].prefix[metric];
  const regions = [];
  for (let start = 1; start <= transcriptLength - regionLength + 1; start += 1) {
    const end = start + regionLength - 1;
    const first = start - 1;
    const score = contained === 1
      ? scoreIndex[guide].rows[first][metric]
      : (prefix[first + contained] - prefix[first]) / contained;
    regions.push({
      start,
      end,
      feature: startFeature(start),
      guide,
      regionLength,
      metric,
      score,
      contained,
      sequence: payload.target.transcript_sequence.slice(start - 1, end),
    });
  }
  return regions;
}

function updateUrl(state) {
  const query = new URLSearchParams({ guide: state.guide, region: state.region, metric: state.metric });
  history.replaceState(null, "", `${location.pathname}?${query}${location.hash}`);
}

function renderPlot(regions, state) {
  const plot = document.getElementById("region-plot");
  if (!regions.length) {
    Plotly.purge(plot);
    plot.innerHTML = '<p class="na">N/A — region length must be between the selected guide length and transcript length.</p>';
    return;
  }
  const custom = regions.map(row => [row.end, row.feature || "N/A", row.guide, row.regionLength, row.metric, row.contained]);
  const shapes = payload.target.annotations.map(item => ({
    type: "rect", xref: "x", yref: "paper", x0: item.start_1based, x1: item.end_1based,
    y0: 0, y1: 1, fillcolor: FEATURE_COLORS[item.feature] || "rgba(0,0,0,.04)", line: { width: 0 }, layer: "below",
  }));
  Plotly.react(plot, [{
    x: regions.map(row => row.start), y: regions.map(row => row.score), customdata: custom,
    mode: "lines", line: { color: "#087f72", width: 2 },
    hovertemplate: "Start %{x}<br>End %{customdata[0]}<br>Feature %{customdata[1]}<br>Guide %{customdata[2]} nt<br>Region %{customdata[3]} nt<br>Metric %{customdata[4]}<br>Mean %{y:.6f}<br>Windows %{customdata[5]}<extra></extra>",
  }], {
    title: `${payload.target.display_name} · ${state.metric} · ${state.guide}-nt guides · ${state.region}-nt regions`,
    xaxis: { title: "Region start position (nt)" }, yaxis: { title: "Mean contained-window score", range: [0, 1] },
    shapes, margin: { t: 70, r: 25, b: 60, l: 65 }, paper_bgcolor: "#fff", plot_bgcolor: "#fff",
  }, { responsive: true, displaylogo: false });
}

function rankedRows(regions, feature, descending) {
  return regions.filter(row => row.feature === feature).sort((a, b) => {
    const scoreOrder = descending ? b.score - a.score : a.score - b.score;
    return scoreOrder || a.start - b.start;
  }).slice(0, 5);
}

function tableRows(rows) {
  const padded = rows.concat(Array(Math.max(0, 5 - rows.length)).fill(null));
  return padded.map((row, index) => row
    ? `<tr><td>${index + 1}</td><td>${row.start}–${row.end}</td><td>${row.score.toFixed(6)}</td><td class="sequence">${row.sequence}</td></tr>`
    : `<tr><td>${index + 1}</td><td colspan="3" class="na">N/A</td></tr>`).join("");
}

function renderTables(regions) {
  const container = document.getElementById("feature-tables");
  container.innerHTML = ["5UTR", "CDS", "3UTR"].map(feature => {
    const top = rankedRows(regions, feature, true);
    const bottom = rankedRows(regions, feature, false);
    return `<article class="feature-card"><h3>${feature}</h3>
      <h4>Top 5</h4><table><thead><tr><th>#</th><th>Region</th><th>Mean</th><th>Sequence</th></tr></thead><tbody>${tableRows(top)}</tbody></table>
      <h4>Bottom 5</h4><table><thead><tr><th>#</th><th>Region</th><th>Mean</th><th>Sequence</th></tr></thead><tbody>${tableRows(bottom)}</tbody></table>
    </article>`;
  }).join("");
}

function update() {
  const guide = selectedGuideLength();
  const region = Number.parseInt(document.getElementById("region-length").value, 10);
  const metric = document.getElementById("metric-mode").value;
  const state = { guide, region, metric };
  const regions = calculateRegions(guide, region, metric);
  document.getElementById("state-message").textContent = regions.length
    ? `${regions.length} valid starts; ${region - guide + 1} fully contained ${guide}-nt windows per region.`
    : `N/A — choose a region length from ${guide} to ${payload.target.transcript_length_nt} nt.`;
  renderPlot(regions, state);
  renderTables(regions);
  if (regions.length) updateUrl(state);
}

async function init() {
  if (window.STAGE11_PAYLOAD) {
    payload = window.STAGE11_PAYLOAD;
  } else {
    const source = document.body.dataset.stage11Source;
    const response = await fetch(source);
    if (!response.ok) throw new Error(`Unable to load Stage 11 data: ${response.status}`);
    payload = await response.json();
  }
  scoreIndex = buildIndex();
  const state = safeState(new URLSearchParams(location.search));
  renderGuideChoices(payload.supported_guide_lengths, state.guide);
  document.getElementById("region-length").value = state.region;
  document.getElementById("metric-mode").value = state.metric;
  document.getElementById("target-summary").textContent = `${payload.target.display_name} · ${payload.target.transcript_length_nt} nt · ${payload.schema_version}`;
  for (const control of document.querySelectorAll("select, input")) control.addEventListener("change", update);
  document.getElementById("region-length").addEventListener("input", update);
  update();
}

init().catch(error => {
  document.getElementById("state-message").textContent = error.message;
  console.error(error);
});
