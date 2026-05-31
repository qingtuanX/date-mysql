const state = {
  datasets: [], comparisons: [], points: [], tableRows: [],
  pointIndexById: new Map(), activePointId: null, pinnedPointId: null,
  immuneColumns: [], immuneRows: [], immuneGroups: [], immuneGroup: "All",
  immuneTruncated: false, plotlyReady: false, immuneActiveCell: null
};

const refs = {
  datasetSelect: document.getElementById("datasetSelect"),
  comparisonSelect: document.getElementById("comparisonSelect"),
  geneSearch: document.getElementById("geneSearch"),
  clearPinBtn: document.getElementById("clearPinBtn"),
  pointInfo: document.getElementById("pointInfo"),
  volcanoChart: document.getElementById("volcanoChart"),
  resultTableBody: document.getElementById("resultTableBody"),
  pointCounter: document.getElementById("pointCounter"),
  immuneGroupSelect: document.getElementById("immuneGroupSelect"),
  immuneSearch: document.getElementById("immuneSearch"),
  immuneMetricSearch: document.getElementById("immuneMetricSearch"),
  immuneHeatmap: document.getElementById("immuneHeatmap"),
  immuneTableHead: document.getElementById("immuneTableHead"),
  immuneTableBody: document.getElementById("immuneTableBody"),
  immuneCounter: document.getElementById("immuneCounter"),
  immunePointInfo: document.getElementById("immunePointInfo"),
  viewDiffBtn: document.getElementById("viewDiffBtn"),
  viewImmuneBtn: document.getElementById("viewImmuneBtn"),
  diffSection: document.getElementById("diffSection"),
  immuneSection: document.getElementById("immuneSection")
};

function escapeHtml(v) { return String(v).replaceAll("&","&amp;").replaceAll("<","&lt;").replaceAll(">","&gt;").replaceAll('"',"&quot;").replaceAll("'","&#039;"); }
function regulationColor(r) { return r==="Up"?"#c0392b":r==="Down"?"#1f8a70":"#8a96a3"; }
function fmt(n,d) { d=d||4; if(n===null||n===undefined||isNaN(Number(n)))return"-"; return Number(n).toFixed(d); }
function fi(v) { if(v===null||v===undefined||v==="")return"-"; if(typeof v==="number")return fmt(v,4); return escapeHtml(v); }

let _plotlyQueue = [];
function whenPlotly(cb) { if(state.plotlyReady) cb(); else _plotlyQueue.push(cb); }

async function requestJson(url) { const r=await fetch(url); if(!r.ok)throw new Error(await r.text()); return r.json(); }

function debounce(fn,w) { let t; return (...a)=>{clearTimeout(t);t=setTimeout(()=>fn(...a),w)}; }

// ===== INIT =====
async function init() {
  window.addEventListener("plotly-loaded",()=>{ state.plotlyReady=true; _plotlyQueue.forEach(f=>f()); _plotlyQueue=[] });
  bindEvents(); setActiveView("diff"); await loadDatasets();
}

function bindEvents() {
  refs.viewDiffBtn.addEventListener("click",()=>setActiveView("diff"));
  refs.viewImmuneBtn.addEventListener("click",()=>setActiveView("immune"));
  refs.datasetSelect.addEventListener("change",async()=>{await loadComparisons();await loadAndRenderAll();await loadImmuneData()});
  refs.comparisonSelect.addEventListener("change",()=>loadAndRenderAll());
  refs.geneSearch.addEventListener("input",debounce(renderTable,250));
  refs.immuneGroupSelect.addEventListener("change",()=>{state.immuneGroup=refs.immuneGroupSelect.value;renderImmuneHeatmap();renderImmuneTable();updateImmuneInfo(null)});
  refs.immuneSearch.addEventListener("input",debounce(()=>{renderImmuneHeatmap();renderImmuneTable();updateImmuneInfo(null)},250));
  refs.immuneMetricSearch.addEventListener("input",debounce(()=>{renderImmuneHeatmap();renderImmuneTable();updateImmuneInfo(null)},250));
  refs.clearPinBtn.addEventListener("click",()=>{state.pinnedPointId=null;if(state.activePointId)focusPointById(state.activePointId,false)});
}

function setActiveView(v) {
  const im = v==="immune";
  refs.diffSection.classList.toggle("is-hidden",im);
  refs.immuneSection.classList.toggle("is-hidden",!im);
  refs.viewDiffBtn.classList.toggle("is-active",!im);
  refs.viewImmuneBtn.classList.toggle("is-active",im);
}

// ===== DATASETS / COMPARISONS =====
async function loadDatasets() {
  state.datasets = await requestJson("/api/datasets");
  refs.datasetSelect.innerHTML = state.datasets.map(d=>`<option>${escapeHtml(d)}</option>`).join("");
  await loadComparisons(); await loadAndRenderAll(); await loadImmuneData();
}

async function loadComparisons() {
  const d = refs.datasetSelect.value;
  state.comparisons = await requestJson(`/api/comparisons?dataset=${encodeURIComponent(d)}`);
  refs.comparisonSelect.innerHTML = state.comparisons.map(c=>`<option>${escapeHtml(c)}</option>`).join("");
}

// ===== DIFF ANALYSIS =====
async function loadAndRenderAll() {
  const d=refs.datasetSelect.value, c=refs.comparisonSelect.value;
  if(!d||!c)return;
  const [points,rows] = await Promise.all([
    requestJson(`/api/points?dataset=${encodeURIComponent(d)}&comparison=${encodeURIComponent(c)}`),
    requestJson(`/api/top?dataset=${encodeURIComponent(d)}&comparison=${encodeURIComponent(c)}&limit=200`)
  ]);
  state.points=points; state.tableRows=rows; state.pointIndexById.clear();
  points.forEach((p,i)=>state.pointIndexById.set(p.id,i));
  state.activePointId=state.pinnedPointId=null;
  refs.pointCounter.textContent=`共 ${points.length} 个点`;
  whenPlotly(renderChart); renderTable(); updatePointInfo(null);
}

function renderChart() {
  const up=[],down=[],ns=[];
  for(const p of state.points){ if(p.regulation==="Up")up.push(p);else if(p.regulation==="Down")down.push(p);else ns.push(p); }
  const traces = [ bt("Down",down), bt("NotSig",ns), bt("Up",up) ];
  const layout={
    margin:{l:58,r:20,t:20,b:58},paper_bgcolor:"#fff",plot_bgcolor:"#fff",
    xaxis:{title:"log2 Fold Change",zeroline:true,zerolinecolor:"#b7c6d5",gridcolor:"#e6eef7"},
    yaxis:{title:"-log10(FDR)",gridcolor:"#e6eef7"},
    legend:{orientation:"h",y:1.1},hovermode:"closest",showlegend:true
  };
  Plotly.newPlot(refs.volcanoChart,traces,layout,{responsive:true,displayModeBar:false});
  refs.volcanoChart.on("plotly_hover",ev=>{if(state.pinnedPointId)return;const id=ev?.points?.[0]?.customdata?.id;if(id)focusPointById(id,false)});
  refs.volcanoChart.on("plotly_click",ev=>{const id=ev?.points?.[0]?.customdata?.id;if(id){state.pinnedPointId=id;focusPointById(id,true)}});
  refs.volcanoChart.on("plotly_unhover",()=>{if(!state.pinnedPointId){state.activePointId=null;highlightRow(null);updatePointInfo(null);clearChartHighlight()}});
}

function bt(name,pts) {
  return {name,type:"scattergl",mode:"markers",x:pts.map(p=>Number(p.log2fc)),y:pts.map(p=>Number(p.negLog10Fdr)),
    marker:{size:pts.map(p=>p.id===state.activePointId?12:7),color:regulationColor(name),opacity:.8,line:{width:pts.map(p=>p.id===state.activePointId?2:0),color:"#0f6cab"}},
    customdata:pts.map(p=>({id:p.id,gene:p.gene,regulation:p.regulation})),
    hovertemplate:"<b>%{customdata.gene}</b><br>log2FC: %{x:.4f}<br>-log10(FDR): %{y:.4f}<extra></extra>"};
}

function renderTable() {
  const kw=refs.geneSearch.value.trim().toLowerCase();
  const rows=kw?state.tableRows.filter(r=>String(r.gene).toLowerCase().includes(kw)):state.tableRows;
  refs.resultTableBody.innerHTML=rows.map(r=>
    `<tr data-id="${r.id}"><td>${escapeHtml(r.gene)}</td><td>${fmt(r.log2fc)}</td><td>${fmt(r.negLog10Fdr)}</td><td>${fmt(r.fdr,6)}</td><td style="color:${regulationColor(r.regulation)};font-weight:600">${r.regulation}</td></tr>`
  ).join("");
  refs.resultTableBody.querySelectorAll("tr").forEach(tr=>{
    tr.addEventListener("mouseenter",()=>{if(state.pinnedPointId)return;focusPointById(Number(tr.dataset.id),false)});
    tr.addEventListener("click",()=>{const id=Number(tr.dataset.id);state.pinnedPointId=id;focusPointById(id,true)});
  });
}

function focusPointById(id,pinned) {
  const p=state.points[state.pointIndexById.get(id)]; if(!p)return;
  state.activePointId=id; updatePointInfo(p,pinned); highlightRow(id); refreshHighlightOnChart(id);
}

function refreshHighlightOnChart(activeId) {
  const gd=refs.volcanoChart,up={["marker.size"]:[],["marker.line.width"]:[]};
  for(let i=0;i<gd.data.length;i++){const t=gd.data[i],sz=[],wd=[];for(const cd of t.customdata){const a=cd.id===activeId;sz.push(a?12:7);wd.push(a?2:0)}up["marker.size"].push(sz);up["marker.line.width"].push(wd)}
  Plotly.restyle(gd,up); const ap=state.points[state.pointIndexById.get(activeId)]; if(!ap)return;
  Plotly.relayout(gd,{annotations:[{x:Number(ap.log2fc),y:Number(ap.negLog10Fdr),text:`${ap.gene} (${fmt(ap.log2fc)}, ${fmt(ap.negLog10Fdr)})`,showarrow:true,arrowhead:2,ax:20,ay:-30,bgcolor:"rgba(15,108,171,.86)",bordercolor:"#fff",font:{color:"#fff",size:12}}]});
}

function clearChartHighlight() { const gd=refs.volcanoChart,up={["marker.size"]:[],["marker.line.width"]:[]};for(let i=0;i<gd.data.length;i++){up["marker.size"].push(gd.data[i].customdata.map(()=>7));up["marker.line.width"].push(gd.data[i].customdata.map(()=>0))}Plotly.restyle(gd,up);Plotly.relayout(gd,{annotations:[]})}
function highlightRow(id) { refs.resultTableBody.querySelectorAll("tr").forEach(t=>t.classList.toggle("active",Number(t.dataset.id)===Number(id))) }
function updatePointInfo(pt,pinned) {
  if(!pt){refs.pointInfo.innerHTML="<p>未选择数据点</p>";return}
  refs.pointInfo.innerHTML=`<div><b>Gene:</b> ${escapeHtml(pt.gene)}</div><div><b>log2FC:</b> ${fmt(pt.log2fc)}</div><div><b>-log10(FDR):</b> ${fmt(pt.negLog10Fdr)}</div><div><b>FDR:</b> ${fmt(pt.fdr,6)}</div><div><b>Status:</b> <span style="color:${regulationColor(pt.regulation)}">${pt.regulation}</span></div>${pinned?'<div style="color:var(--accent)">📍 已固定</div>':''}`;
}

// ===== IMMUNE =====
function inferImmuneGroup(sid) {
  if(!sid)return"Unknown"; const t=String(sid),parts=t.split("_");
  if(parts.length>1){const tail=parts[parts.length-1];if(/^[A-Za-z]+$/.test(tail))return tail}
  const m=t.match(/^[A-Za-z]+/); return m?m[0]:"Unknown";
}

function buildImmuneGroups(rows) { const s=new Set(); rows.forEach(r=>s.add(inferImmuneGroup(r.ID))); return [...s].sort((a,b)=>a.localeCompare(b)); }

function getFilteredImmuneRows() {
  const kw=refs.immuneSearch.value.trim().toLowerCase();
  return state.immuneRows.filter(r=>{
    if(state.immuneGroup!=="All"&&inferImmuneGroup(r.ID)!==state.immuneGroup)return false;
    return kw?String(r.ID).toLowerCase().includes(kw):true;
  });
}

function getFilteredImmuneMetrics() {
  const kw=refs.immuneMetricSearch.value.trim().toLowerCase();
  const m=state.immuneColumns.filter(c=>c!=="ID");
  return kw?m.filter(c=>c.toLowerCase().includes(kw)):m;
}

function nRow(vals) {
  const nums=vals.filter(v=>typeof v==="number"&&isFinite(v));
  if(!nums.length)return vals.map(()=>0);
  const mean=nums.reduce((a,b)=>a+b,0)/nums.length;
  const vr=nums.reduce((a,b)=>a+(b-mean)**2,0)/nums.length;
  const sd=Math.sqrt(vr)||1;
  return vals.map(v=>(typeof v==="number"&&isFinite(v))?(v-mean)/sd:0);
}

async function loadImmuneData() {
  const d=refs.datasetSelect.value; if(!d)return;
  try {
    const p=await requestJson(`/api/immune?dataset=${encodeURIComponent(d)}&limit=2000`);
    state.immuneColumns=p.columns||[]; state.immuneRows=p.rows||[]; state.immuneTruncated=!!p.truncated;
    state.immuneGroups=buildImmuneGroups(state.immuneRows);
    if(!state.immuneGroups.includes(state.immuneGroup))state.immuneGroup="All";
    updateImmuneGroupOptions(); renderImmuneHeatmap(); renderImmuneTable(); updateImmuneInfo(null);
  } catch(e) {
    state.immuneColumns=[]; state.immuneRows=[]; state.immuneGroups=[]; state.immuneGroup="All"; state.immuneTruncated=false;
    refs.immuneTableHead.innerHTML=""; refs.immuneTableBody.innerHTML=""; refs.immuneCounter.textContent=`加载失败: ${e.message}`;
    refs.immuneHeatmap.innerHTML=""; updateImmuneInfo(null);
  }
}

function updateImmuneGroupOptions() {
  const opts=["All",...state.immuneGroups];
  refs.immuneGroupSelect.innerHTML=opts.map(g=>`<option value="${escapeHtml(g)}">${escapeHtml(g)}</option>`).join("");
  refs.immuneGroupSelect.value=state.immuneGroup;
}

function renderImmuneHeatmap() {
  const metrics=getFilteredImmuneMetrics();
  if(!state.immuneRows.length||!metrics.length){refs.immuneHeatmap.innerHTML="<p class='loading'>暂无可视化数据</p>";return}
  const rows=getFilteredImmuneRows(); if(!rows.length){refs.immuneHeatmap.innerHTML="<p>未找到样本</p>";return}

  const grouped=new Map(); rows.forEach(r=>{const g=inferImmuneGroup(r.ID);if(!grouped.has(g))grouped.set(g,[]);grouped.get(g).push(r)});
  const gOrder=state.immuneGroup==="All"?state.immuneGroups:[state.immuneGroup];
  const ordRows=[], segs=[]; let cur=0;
  gOrder.forEach(g=>{const items=grouped.get(g)||[];items.sort((a,b)=>String(a.ID).localeCompare(String(b.ID)));if(items.length){segs.push({group:g,start:cur,end:cur+items.length-1});ordRows.push(...items);cur+=items.length}});

  const xLab=ordRows.map(r=>r.ID);
  const z=metrics.map(m=>nRow(ordRows.map(r=>r[m])));
  const trace={type:"heatmap",x:xLab,y:metrics,z,colorscale:[[0,"#2c7bb6"],[.5,"#f7f7f7"],[1,"#d7191c"]],zmid:0,hovertemplate:"<b>%{x}</b><br>%{y}: %{z:.3f}<extra></extra>"};

  const n=xLab.length, ann=[]; segs.forEach(s=>{ann.push({xref:"paper",yref:"paper",x:(s.start+s.end+1)/(2*n),y:1.06,text:s.group,showarrow:false,font:{size:12,color:"#0f6cab"}})});

  whenPlotly(()=>{
    Plotly.newPlot(refs.immuneHeatmap,[trace],{margin:{l:180,r:20,t:70,b:120},paper_bgcolor:"#fff",plot_bgcolor:"#fff",xaxis:{tickangle:-45,automargin:true},yaxis:{automargin:true},annotations:ann},{responsive:true,displayModeBar:false});
    refs.immuneHeatmap.on("plotly_click",ev=>{
      const pt=ev?.points?.[0]; if(!pt)return;
      const sid=pt.x, metric=pt.y, val=pt.z;
      updateImmuneInfo({sample:sid,metric:metric,value:val});
    });
  });
}

function updateImmuneInfo(info) {
  if(!info){refs.immunePointInfo.innerHTML="<p>点击热力图查看指标</p>";return}
  const row=state.immuneRows.find(r=>String(r.ID)===info.sample);
  let html=`<div><b>样本:</b> ${escapeHtml(info.sample)}</div><div><b>指标:</b> ${escapeHtml(info.metric)}</div><div><b>Z-score:</b> ${fmt(info.value,3)}</div>`;
  if(row){html+=`<hr style="margin:8px 0;border-color:var(--line)">`;getFilteredImmuneMetrics().slice(0,10).forEach(m=>{html+=`<div style="font-size:.82em${m===info.metric?';color:var(--accent);font-weight:600':''}">${escapeHtml(m)}: ${fi(row[m])}</div>`})}
  refs.immunePointInfo.innerHTML=html;
}

function renderImmuneTable() {
  const metrics=getFilteredImmuneMetrics();
  const cols=["Group","ID",...metrics];
  if(!cols.length){refs.immuneTableHead.innerHTML="";refs.immuneTableBody.innerHTML="";return}
  const rows=getFilteredImmuneRows(),limit=200,view=rows.slice(0,limit);
  const tn=(state.immuneTruncated?" (API截断)":"")+(rows.length>limit?" (显示前200)":"");
  refs.immuneCounter.textContent=`共 ${rows.length} 个样本${tn}`;
  refs.immuneTableHead.innerHTML=cols.map(c=>`<th>${escapeHtml(c)}</th>`).join("");
  refs.immuneTableBody.innerHTML=view.map(r=>{
    const g=inferImmuneGroup(r.ID),cells=[g,r.ID,...metrics.map(m=>r[m])];
    return `<tr data-sid="${escapeHtml(r.ID)}">${cells.map(v=>`<td>${fi(v)}</td>`).join("")}</tr>`;
  }).join("");
  refs.immuneTableBody.querySelectorAll("tr").forEach(tr=>{
    tr.addEventListener("click",()=>{
      const sid=tr.dataset.sid;
      updateImmuneInfo({sample:sid,metric:"(查看表格)",value:null});
    });
  });
}

init();
