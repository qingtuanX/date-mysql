const state = {
  datasets: [],
  comparisons: [],
  points: [],
  tableRows: [],
  pointIndexById: new Map(),
  activePointId: null,
  pinnedPointId: null
};

const refs = {
  datasetSelect: document.getElementById("datasetSelect"),
  comparisonSelect: document.getElementById("comparisonSelect"),
  geneSearch: document.getElementById("geneSearch"),
  clearPinBtn: document.getElementById("clearPinBtn"),
  pointInfo: document.getElementById("pointInfo"),
  volcanoChart: document.getElementById("volcanoChart"),
  resultTableBody: document.getElementById("resultTableBody"),
  pointCounter: document.getElementById("pointCounter")
};

function escapeHtml(value) {
  return String(value)
    .replaceAll("&", "&amp;")
    .replaceAll("<", "&lt;")
    .replaceAll(">", "&gt;")
    .replaceAll('"', "&quot;")
    .replaceAll("'", "&#039;");
}

function regulationColor(regulation) {
  if (regulation === "Up") return "#c0392b";
  if (regulation === "Down") return "#1f8a70";
  return "#8a96a3";
}

function fmt(number, digits = 4) {
  if (number === null || number === undefined || Number.isNaN(Number(number))) return "-";
  return Number(number).toFixed(digits);
}

async function requestJson(url) {
  const resp = await fetch(url);
  if (!resp.ok) {
    const msg = await resp.text();
    throw new Error(msg || `Request failed: ${resp.status}`);
  }
  return resp.json();
}

async function init() {
  bindEvents();
  await loadDatasets();
}

function bindEvents() {
  refs.datasetSelect.addEventListener("change", async () => {
    await loadComparisons();
    await loadAndRenderAll();
  });

  refs.comparisonSelect.addEventListener("change", async () => {
    await loadAndRenderAll();
  });

  refs.geneSearch.addEventListener("input", debounce(async () => {
    renderTable();
  }, 250));

  refs.clearPinBtn.addEventListener("click", () => {
    state.pinnedPointId = null;
    if (state.activePointId) {
      focusPointById(state.activePointId, false);
    }
  });
}

async function loadDatasets() {
  state.datasets = await requestJson("/api/datasets");
  refs.datasetSelect.innerHTML = state.datasets
    .map((dataset) => `<option value="${escapeHtml(dataset)}">${escapeHtml(dataset)}</option>`)
    .join("");

  await loadComparisons();
  await loadAndRenderAll();
}

async function loadComparisons() {
  const dataset = refs.datasetSelect.value;
  state.comparisons = await requestJson(`/api/comparisons?dataset=${encodeURIComponent(dataset)}`);
  refs.comparisonSelect.innerHTML = state.comparisons
    .map((comparison) => `<option value="${escapeHtml(comparison)}">${escapeHtml(comparison)}</option>`)
    .join("");
}

async function loadAndRenderAll() {
  const dataset = refs.datasetSelect.value;
  const comparison = refs.comparisonSelect.value;
  if (!dataset || !comparison) return;

  const [points, rows] = await Promise.all([
    requestJson(
      `/api/points?dataset=${encodeURIComponent(dataset)}&comparison=${encodeURIComponent(comparison)}`
    ),
    requestJson(
      `/api/top?dataset=${encodeURIComponent(dataset)}&comparison=${encodeURIComponent(comparison)}&limit=200`
    )
  ]);

  state.points = points;
  state.tableRows = rows;
  state.pointIndexById.clear();
  state.points.forEach((point, idx) => {
    state.pointIndexById.set(point.id, idx);
  });
  state.activePointId = null;
  state.pinnedPointId = null;

  refs.pointCounter.textContent = `共 ${state.points.length} 个点`;
  renderChart();
  renderTable();
  updatePointInfo(null);
}

function renderChart() {
  const up = [];
  const down = [];
  const ns = [];
  for (const point of state.points) {
    if (point.regulation === "Up") up.push(point);
    else if (point.regulation === "Down") down.push(point);
    else ns.push(point);
  }

  const traces = [
    buildTrace("Down", down),
    buildTrace("NotSig", ns),
    buildTrace("Up", up)
  ];

  const layout = {
    margin: { l: 58, r: 20, t: 20, b: 58 },
    paper_bgcolor: "#ffffff",
    plot_bgcolor: "#ffffff",
    xaxis: {
      title: "log2 Fold Change (X)",
      zeroline: true,
      zerolinecolor: "#b7c6d5",
      gridcolor: "#e6eef7"
    },
    yaxis: {
      title: "-log10(FDR) (Y)",
      gridcolor: "#e6eef7"
    },
    legend: { orientation: "h", y: 1.1 },
    hovermode: "closest",
    showlegend: true
  };

  Plotly.newPlot(refs.volcanoChart, traces, layout, { responsive: true, displayModeBar: true });

  refs.volcanoChart.on("plotly_hover", (eventData) => {
    if (state.pinnedPointId) return;
    const id = eventData?.points?.[0]?.customdata?.id;
    if (!id) return;
    focusPointById(id, false);
  });

  refs.volcanoChart.on("plotly_click", (eventData) => {
    const id = eventData?.points?.[0]?.customdata?.id;
    if (!id) return;
    state.pinnedPointId = id;
    focusPointById(id, true);
  });

  refs.volcanoChart.on("plotly_unhover", () => {
    if (!state.pinnedPointId) {
      state.activePointId = null;
      highlightRow(null);
      updatePointInfo(null);
      clearChartHighlight();
    }
  });
}

function buildTrace(name, points) {
  return {
    name,
    type: "scattergl",
    mode: "markers",
    x: points.map((p) => Number(p.log2fc)),
    y: points.map((p) => Number(p.negLog10Fdr)),
    marker: {
      size: points.map((p) => (p.id === state.activePointId ? 12 : 7)),
      color: regulationColor(name),
      opacity: 0.8,
      line: {
        width: points.map((p) => (p.id === state.activePointId ? 2 : 0)),
        color: "#0f6cab"
      }
    },
    customdata: points.map((p) => ({ id: p.id, gene: p.gene, regulation: p.regulation })),
    hovertemplate:
      "<b>%{customdata.gene}</b><br>Regulation: %{customdata.regulation}<br>log2FC: %{x:.4f}<br>-log10(FDR): %{y:.4f}<extra></extra>"
  };
}

function renderTable() {
  const keyword = refs.geneSearch.value.trim().toLowerCase();
  const rows = keyword
    ? state.tableRows.filter((row) => String(row.gene).toLowerCase().includes(keyword))
    : state.tableRows;

  refs.resultTableBody.innerHTML = rows
    .map(
      (row) => `
      <tr data-id="${row.id}">
        <td>${escapeHtml(row.gene)}</td>
        <td>${fmt(row.log2fc)}</td>
        <td>${fmt(row.negLog10Fdr)}</td>
        <td>${fmt(row.fdr, 6)}</td>
        <td><span style="color:${regulationColor(row.regulation)};font-weight:600;">${row.regulation}</span></td>
      </tr>
      `
    )
    .join("");

  const trs = refs.resultTableBody.querySelectorAll("tr");
  for (const tr of trs) {
    tr.addEventListener("mouseenter", () => {
      if (state.pinnedPointId) return;
      const id = Number(tr.dataset.id);
      focusPointById(id, false);
    });

    tr.addEventListener("click", () => {
      const id = Number(tr.dataset.id);
      state.pinnedPointId = id;
      focusPointById(id, true);
    });
  }
}

function focusPointById(id, pinned) {
  const point = state.points[state.pointIndexById.get(id)];
  if (!point) return;
  state.activePointId = id;

  updatePointInfo(point, pinned);
  highlightRow(id);
  refreshHighlightOnChart(id);
}

function refreshHighlightOnChart(activeId) {
  const gd = refs.volcanoChart;
  const updates = {
    "marker.size": [],
    "marker.line.width": []
  };

  for (let i = 0; i < gd.data.length; i += 1) {
    const trace = gd.data[i];
    const sizes = [];
    const widths = [];
    for (const cd of trace.customdata) {
      const isActive = cd.id === activeId;
      sizes.push(isActive ? 12 : 7);
      widths.push(isActive ? 2 : 0);
    }
    updates["marker.size"].push(sizes);
    updates["marker.line.width"].push(widths);
  }

  Plotly.restyle(gd, updates);

  const activePoint = state.points[state.pointIndexById.get(activeId)];
  if (!activePoint) return;
  Plotly.relayout(gd, {
    annotations: [
      {
        x: Number(activePoint.log2fc),
        y: Number(activePoint.negLog10Fdr),
        text: `${activePoint.gene} (${fmt(activePoint.log2fc)}, ${fmt(activePoint.negLog10Fdr)})`,
        showarrow: true,
        arrowhead: 2,
        ax: 20,
        ay: -30,
        bgcolor: "rgba(15,108,171,0.86)",
        bordercolor: "#fff",
        font: { color: "#fff", size: 12 }
      }
    ]
  });
}

function clearChartHighlight() {
  const gd = refs.volcanoChart;
  const updates = {
    "marker.size": [],
    "marker.line.width": []
  };

  for (let i = 0; i < gd.data.length; i += 1) {
    const trace = gd.data[i];
    updates["marker.size"].push(trace.customdata.map(() => 7));
    updates["marker.line.width"].push(trace.customdata.map(() => 0));
  }

  Plotly.restyle(gd, updates);
  Plotly.relayout(gd, { annotations: [] });
}

function highlightRow(id) {
  const trs = refs.resultTableBody.querySelectorAll("tr");
  for (const tr of trs) {
    tr.classList.toggle("active", Number(tr.dataset.id) === Number(id));
  }
}

function updatePointInfo(point, pinned = false) {
  if (!point) {
    refs.pointInfo.innerHTML = "<p>未选择数据点</p>";
    return;
  }

  refs.pointInfo.innerHTML = `
    <div><b>Gene:</b> ${escapeHtml(point.gene)}</div>
    <div><b>Dataset:</b> ${escapeHtml(point.dataset)}</div>
    <div><b>Comparison:</b> ${escapeHtml(point.comparison)}</div>
    <div><b>X (log2FC):</b> ${fmt(point.log2fc)}</div>
    <div><b>Y (-log10FDR):</b> ${fmt(point.negLog10Fdr)}</div>
    <div><b>FDR:</b> ${fmt(point.fdr, 6)}</div>
    <div><b>Status:</b> <span style="color:${regulationColor(point.regulation)}">${point.regulation}</span></div>
    <div><b>模式:</b> ${pinned ? "已固定(点击取消固定按钮解除)" : "悬浮联动"}</div>
  `;
}

function debounce(fn, wait) {
  let timer = null;
  return (...args) => {
    if (timer) clearTimeout(timer);
    timer = setTimeout(() => fn(...args), wait);
  };
}

init().catch((error) => {
  refs.pointInfo.innerHTML = `<p>加载失败: ${escapeHtml(error.message)}</p>`;
});
