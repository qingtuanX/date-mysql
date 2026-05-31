require("dotenv").config();

const express = require("express");
const compression = require("compression");
const fs = require("fs");
const path = require("path");
const readline = require("readline");
const pool = require("./db");

const app = express();
const port = Number(process.env.PORT || 3000);
const projectRoot = path.resolve(__dirname, "..");

function parseCsvLine(line) {
  const values = [];
  let current = "";
  let inQuotes = false;
  for (let i = 0; i < line.length; i += 1) {
    const char = line[i];
    if (char === '"') {
      if (inQuotes && line[i + 1] === '"') { current += '"'; i += 1; }
      else { inQuotes = !inQuotes; }
      continue;
    }
    if (char === "," && !inQuotes) { values.push(current.trim()); current = ""; continue; }
    current += char;
  }
  values.push(current.trim());
  return values.map((value) => value.replace(/^"|"$/g, ""));
}

function toNumberOrNull(value) {
  if (value === undefined || value === null) return null;
  const trimmed = String(value).trim();
  if (trimmed === "" || trimmed.toUpperCase() === "NA") return null;
  const n = Number(trimmed);
  return Number.isFinite(n) ? n : trimmed;
}

app.use(compression({ threshold: 512 }));
app.use(express.json());
app.use(express.static(path.join(__dirname, "public"), { maxAge: "1d" }));

app.get("/api/datasets", async (_req, res) => {
  try {
    const [rows] = await pool.query("SELECT DISTINCT dataset FROM limma_results ORDER BY dataset");
    res.json(rows.map((row) => row.dataset));
  } catch (error) {
    res.status(500).json({ message: "Failed", error: error.message });
  }
});

app.get("/api/comparisons", async (req, res) => {
  const { dataset } = req.query;
  if (!dataset) { res.status(400).json({ message: "dataset required" }); return; }
  try {
    const [rows] = await pool.query("SELECT DISTINCT comparison FROM limma_results WHERE dataset = ? ORDER BY comparison", [dataset]);
    res.json(rows.map((row) => row.comparison));
  } catch (error) {
    res.status(500).json({ message: "Failed", error: error.message });
  }
});

app.get("/api/points", async (req, res) => {
  const { dataset, comparison } = req.query;
  if (!dataset || !comparison) { res.status(400).json({ message: "dataset and comparison required" }); return; }
  let sql = `SELECT id,dataset,comparison,gene,log2fc,p_value AS pValue,fdr,neg_log10_fdr AS negLog10Fdr,
    CASE WHEN fdr < 0.05 AND log2fc >= 1 THEN 'Up' WHEN fdr < 0.05 AND log2fc <= -1 THEN 'Down' ELSE 'NotSig' END AS regulation
    FROM limma_results WHERE dataset = ? AND comparison = ? ORDER BY ABS(log2fc) DESC, fdr ASC`;
  try {
    const [rows] = await pool.query(sql, [dataset, comparison]);
    res.json(rows);
  } catch (error) {
    res.status(500).json({ message: "Failed", error: error.message });
  }
});

app.get("/api/top", async (req, res) => {
  const { dataset, comparison, limit } = req.query;
  if (!dataset || !comparison) { res.status(400).json({ message: "dataset and comparison required" }); return; }
  const topN = Math.max(10, Math.min(500, Number(limit || 200)));
  const sql = `SELECT id,gene,log2fc,p_value AS pValue,fdr,neg_log10_fdr AS negLog10Fdr,
    CASE WHEN fdr < 0.05 AND log2fc >= 1 THEN 'Up' WHEN fdr < 0.05 AND log2fc <= -1 THEN 'Down' ELSE 'NotSig' END AS regulation
    FROM limma_results WHERE dataset = ? AND comparison = ? ORDER BY fdr ASC, ABS(log2fc) DESC LIMIT ?`;
  try {
    const [rows] = await pool.query(sql, [dataset, comparison, topN]);
    res.json(rows);
  } catch (error) {
    res.status(500).json({ message: "Failed", error: error.message });
  }
});

app.get("/api/immune", async (req, res) => {
  const { dataset } = req.query;
  if (!dataset) { res.status(400).json({ message: "dataset required" }); return; }
  const filePath = path.join(projectRoot, "mRNA", "免疫浸润", dataset, "tme_combine.csv");
  try { await fs.promises.access(filePath, fs.constants.R_OK); } catch (_) {
    res.status(404).json({ message: "Immune file not found" }); return;
  }
  const maxRows = 2000;
  try {
    const stream = fs.createReadStream(filePath, { encoding: "utf8" });
    const rl = readline.createInterface({ input: stream, crlfDelay: Infinity });
    const rows = [], columns = []; let lineNo = 0;
    for await (const line of rl) {
      if (!line.trim()) continue; lineNo += 1;
      if (lineNo === 1) { columns.push(...parseCsvLine(line)); continue; }
      const values = parseCsvLine(line); if (!values.length) continue;
      const row = {}; columns.forEach((col, idx) => { const raw = values[idx] ?? ""; row[col] = col === "ID" ? raw : toNumberOrNull(raw); });
      rows.push(row);
      if (rows.length >= maxRows) break;
    }
    res.json({ columns, rows, truncated: rows.length >= maxRows });
  } catch (error) {
    res.status(500).json({ message: "Failed", error: error.message });
  }
});

app.listen(port, () => { console.log(`Server running at http://localhost:${port}`); });
