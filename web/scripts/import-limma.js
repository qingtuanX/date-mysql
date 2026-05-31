require("dotenv").config();

const fs = require("fs");
const path = require("path");
const readline = require("readline");
const pool = require("../db");

const projectRoot = path.resolve(__dirname, "..", "..");
const mrnaRoot = path.join(projectRoot, "mRNA");
const datasetPattern = /^GSE\d+$/;
const resultPattern = /^(.+)_limma_test_result\.(.+)\.txt$/;

async function findResultFiles() {
  const entries = await fs.promises.readdir(mrnaRoot, { withFileTypes: true });
  const files = [];

  // Look in 差异化分析 subdirectory (corrected location)
  const diffDir = path.join(mrnaRoot, "差异化分析");
  let diffEntries = [];
  try {
    diffEntries = await fs.promises.readdir(diffDir, { withFileTypes: true });
  } catch (e) {
    diffEntries = [];
  }

  for (const entry of diffEntries) {
    if (!entry.isDirectory() || !datasetPattern.test(entry.name)) continue;
    const datasetDir = path.join(diffDir, entry.name);
    const children = await fs.promises.readdir(datasetDir, { withFileTypes: true });
    for (const child of children) {
      if (!child.isFile()) continue;
      const m = child.name.match(resultPattern);
      if (!m) continue;
      files.push({
        dataset: entry.name,
        comparison: m[2],
        filePath: path.join(datasetDir, child.name)
      });
    }
  }

  return files;
}

function toNumberOrNull(value) {
  if (value === undefined || value === null || value === "") return null;
  const n = Number(value);
  return Number.isFinite(n) ? n : null;
}

function calcNegLog10(fdr) {
  if (fdr === null || fdr === undefined) return null;
  const n = Number(fdr);
  if (!Number.isFinite(n) || n <= 0) return 320;
  return -Math.log10(n);
}

async function parseResultFile(meta) {
  const rows = [];
  const stream = fs.createReadStream(meta.filePath, { encoding: "utf8" });
  const rl = readline.createInterface({ input: stream, crlfDelay: Infinity });

  let lineNo = 0;
  for await (const line of rl) {
    lineNo += 1;
    if (!line.trim()) continue;
    if (lineNo === 1 && line.includes("log2FC")) continue;

    const cols = line.split("\t");
    if (cols.length < 6) continue;

    const gene = cols[0].trim();
    const log2fc = toNumberOrNull(cols[1]);
    const tValue = toNumberOrNull(cols[2]);
    const bValue = toNumberOrNull(cols[3]);
    const pValue = toNumberOrNull(cols[4]);
    const fdr = toNumberOrNull(cols[5]);

    if (!gene || log2fc === null || fdr === null) continue;

    rows.push([
      meta.dataset,
      meta.comparison,
      gene,
      log2fc,
      tValue,
      bValue,
      pValue,
      fdr,
      calcNegLog10(fdr)
    ]);
  }

  return rows;
}

async function importFile(meta) {
  const rows = await parseResultFile(meta);
  if (rows.length === 0) {
    console.log(`[SKIP] ${meta.dataset} ${meta.comparison}: no data`);
    return;
  }

  await pool.query("DELETE FROM limma_results WHERE dataset = ? AND comparison = ?", [
    meta.dataset,
    meta.comparison
  ]);

  const chunkSize = 1000;
  for (let i = 0; i < rows.length; i += chunkSize) {
    const chunk = rows.slice(i, i + chunkSize);
    await pool.query(
      `
      INSERT INTO limma_results
      (dataset, comparison, gene, log2fc, t_value, b_value, p_value, fdr, neg_log10_fdr)
      VALUES ?
      `,
      [chunk]
    );
  }

  console.log(`[OK] ${meta.dataset} ${meta.comparison}: ${rows.length} rows`);
}

async function main() {
  const files = await findResultFiles();
  if (files.length === 0) {
    console.log("No limma result files found.");
    return;
  }

  console.log(`Found ${files.length} result files.`);
  for (const file of files) {
    await importFile(file);
  }
}

main()
  .catch((error) => {
    console.error(error);
    process.exitCode = 1;
  })
  .finally(async () => {
    await pool.end();
  });
