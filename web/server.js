require("dotenv").config();

const express = require("express");
const path = require("path");
const pool = require("./db");

const app = express();
const port = Number(process.env.PORT || 3000);

app.use(express.json());
app.use(express.static(path.join(__dirname, "public")));

app.get("/api/datasets", async (_req, res) => {
  try {
    const [rows] = await pool.query(
      "SELECT DISTINCT dataset FROM limma_results ORDER BY dataset"
    );
    res.json(rows.map((row) => row.dataset));
  } catch (error) {
    res.status(500).json({ message: "Failed to fetch datasets", error: error.message });
  }
});

app.get("/api/comparisons", async (req, res) => {
  const { dataset } = req.query;
  if (!dataset) {
    res.status(400).json({ message: "dataset is required" });
    return;
  }

  try {
    const [rows] = await pool.query(
      "SELECT DISTINCT comparison FROM limma_results WHERE dataset = ? ORDER BY comparison",
      [dataset]
    );
    res.json(rows.map((row) => row.comparison));
  } catch (error) {
    res.status(500).json({ message: "Failed to fetch comparisons", error: error.message });
  }
});

app.get("/api/points", async (req, res) => {
  const { dataset, comparison, q } = req.query;
  if (!dataset || !comparison) {
    res.status(400).json({ message: "dataset and comparison are required" });
    return;
  }

  let sql = `
    SELECT
      id,
      dataset,
      comparison,
      gene,
      log2fc,
      p_value AS pValue,
      fdr,
      neg_log10_fdr AS negLog10Fdr,
      CASE
        WHEN fdr < 0.05 AND log2fc >= 1 THEN 'Up'
        WHEN fdr < 0.05 AND log2fc <= -1 THEN 'Down'
        ELSE 'NotSig'
      END AS regulation
    FROM limma_results
    WHERE dataset = ? AND comparison = ?
  `;
  const params = [dataset, comparison];

  if (q) {
    sql += " AND gene LIKE ?";
    params.push(`%${q}%`);
  }

  sql += " ORDER BY ABS(log2fc) DESC, fdr ASC";

  try {
    const [rows] = await pool.query(sql, params);
    res.json(rows);
  } catch (error) {
    res.status(500).json({ message: "Failed to fetch points", error: error.message });
  }
});

app.get("/api/top", async (req, res) => {
  const { dataset, comparison, limit } = req.query;
  if (!dataset || !comparison) {
    res.status(400).json({ message: "dataset and comparison are required" });
    return;
  }

  const topN = Math.max(10, Math.min(500, Number(limit || 200)));
  const sql = `
    SELECT
      id,
      gene,
      log2fc,
      p_value AS pValue,
      fdr,
      neg_log10_fdr AS negLog10Fdr,
      CASE
        WHEN fdr < 0.05 AND log2fc >= 1 THEN 'Up'
        WHEN fdr < 0.05 AND log2fc <= -1 THEN 'Down'
        ELSE 'NotSig'
      END AS regulation
    FROM limma_results
    WHERE dataset = ? AND comparison = ?
    ORDER BY fdr ASC, ABS(log2fc) DESC
    LIMIT ?
  `;

  try {
    const [rows] = await pool.query(sql, [dataset, comparison, topN]);
    res.json(rows);
  } catch (error) {
    res.status(500).json({ message: "Failed to fetch table rows", error: error.message });
  }
});

app.get("/api/point/:id", async (req, res) => {
  const { id } = req.params;
  try {
    const [rows] = await pool.query(
      `
      SELECT
        id,
        dataset,
        comparison,
        gene,
        log2fc,
        p_value AS pValue,
        fdr,
        neg_log10_fdr AS negLog10Fdr
      FROM limma_results
      WHERE id = ?
      `,
      [id]
    );
    if (rows.length === 0) {
      res.status(404).json({ message: "Point not found" });
      return;
    }
    res.json(rows[0]);
  } catch (error) {
    res.status(500).json({ message: "Failed to fetch point", error: error.message });
  }
});

app.listen(port, () => {
  console.log(`Server running at http://localhost:${port}`);
});
