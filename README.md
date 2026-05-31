# date-mysql — 结核病宿主 mRNA 转录组分析平台

> Tuberculosis Host mRNA Transcriptomics — Differential Expression, Immune Infiltration & Interactive Visualization

---

## 在线访问 | Live Demo

| 站点 | 地址 |
|------|------|
| **交互式分析平台** (部署服务器) | http://49.235.115.164 |
| **图表展示页** (GitHub Pages) | https://qingtuanx.github.io/date-mysql/ |

---

## 项目简介 | Overview

面向结核病/肺病 mRNA 转录组数据，整合 **R limma-voom 差异分析**、**IOBR 免疫浸润去卷积** 和 **MySQL + Node.js 交互可视化**。

| 层级 | 技术栈 |
|------|--------|
| **统计分析** | R (≥4.0), limma, edgeR, IOBR, ggplot2, pheatmap |
| **后端** | Node.js, Express, MySQL2 |
| **前端** | 原生 HTML/CSS/JS, Plotly.js |
| **数据库** | MySQL 8.0 |

---

## 项目结构 | Structure

```
date-mysql/
├── mRNA/
│   ├── GSE148036/                 # 原始数据 | Lung tissue (Normal/TB/AD/SA)
│   ├── GSE158767/                 # 原始数据 | Lung tissue (H/L groups)
│   ├── GSE99374/                  # 原始数据 | PBMC (Healthy/LTBI)
│   ├── 差异化分析/                  # limma-voom/trend 差异分析结果
│   │   ├── GSE148036/             #   → *_limma_test_result.*.txt + volcano/heatmap plots
│   │   ├── GSE158767/
│   │   └── GSE99374/
│   ├── 免疫浸润/                    # TPM → IOBR 免疫细胞去卷积结果
│   │   ├── GSE148036/             #   → 5-7 methods CSV + tme_combine + heatmap
│   │   ├── GSE158767/
│   │   └── GSE99374/
│   ├── SCI绘图模板/                 # 多组/单组火山图 + limma R Markdown 模板
│   ├── README_diff.md             # 差异分析详细文档
│   ├── README_immune.md           # 免疫浸润详细文档
│   └── index.html                 # 静态图表浏览页面
├── web/                           # Web 交互应用
│   ├── public/                    # 前端 (HTML/CSS/JS + Plotly 火山图)
│   ├── scripts/import-limma.js    # limma 结果 → MySQL 导入
│   ├── sql/init.sql               # 数据库建表 DDL
│   └── server.js                  # Express API (datasets/comparisons/points/immune)
├── result/                        # 独立配对比较结果 (DESeq2 等)
└── README.md
```

---

## 差异分析 | Differential Expression

**方法**: limma-voom (counts) / limma-trend (FPKM)

| 数据集 | 样本数 | 分组 | 方法 | 分析基因数 |
|--------|--------|------|------|-----------|
| GSE148036 | 20 | Normal, TB, Adenocarcinoma, Sarcoidosis | voom | 18,235 |
| GSE158767 | 10 | H, L | limma-trend | 3,518 |
| GSE99374 | 40 | Healthy, Latent TB infection | voom | 14,097 |

对比设计: **Target group vs All Others** (non-pairwise)

- `log2(counts+1)` → **voom (mean-variance modeling)** ✓
- 低表达过滤: CPM > 1 (≥2 样本) / mean FPKM ≥ 1
- FDR: Benjamini-Hochberg 校正

📖 [详细文档](mRNA/README_diff.md)

---

## 免疫浸润 | Immune Infiltration

**方法**: counts/FPKM → **TPM 转换** → IOBR 多算法去卷积

| 方法 | 说明 | GSE148036 | GSE158767 | GSE99374 |
|------|------|:---:|:---:|:---:|
| MCPcounter | Marker gene 绝对分数 | ✓ | ✓ | ✓ |
| ESTIMATE | 基质/免疫综合评分 | ✓ | ✓ | ✓ |
| EPIC | 免疫+基质细胞比例 | ✓ | ✓ | ✓ |
| IPS | 免疫表型分数 | ✓ | ✓ | ✓ |
| CIBERSORT | 22 种免疫细胞 LM22 | ✓ | ✓ | ✓ |
| TIMER | 6 种免疫细胞 (LUAD) | ✓ | ✓ | — |
| quanTIseq | 10 种细胞绝对比例 | ✓ | — | — |

- TPM 转换: `RPK = counts / (gene_length/1000)`, `TPM = RPK / sum(RPK) × 10⁶`
- 基因长度来源: `Human.GRCh38.p13.annot.tsv` / `org.Hs.eg.db`
- 已删除 xCell (64 种细胞类型噪声大，含与肺病无关的细胞)

📖 [详细文档](mRNA/README_immune.md)

---

## Web 应用 | Web App

交互式火山图 + 免疫浸润热力图，MySQL 驱动联动查询。

```
cd web
npm install
npm start               # → http://localhost:3000
npm run import:data     # 差异分析结果 → MySQL 导入
```

**API 端点**:
| 端点 | 说明 |
|------|------|
| `GET /api/datasets` | 获取数据集列表 |
| `GET /api/comparisons?dataset=` | 获取某数据集比较组 |
| `GET /api/points?dataset=&comparison=` | 获取火山图散点 |
| `GET /api/top?dataset=&comparison=` | Top N 差异基因表 |
| `GET /api/immune?dataset=` | 免疫浸润数据 + 分组热力图 |

---

## 快速开始 | Quick Start

**1. 差异分析**（需要有 R + limma/edgeR）:
```bash
cd mRNA/差异化分析
Rscript GSE148036/run_limma_gse148036.R
Rscript GSE158767/run_limma_gse158767.R
Rscript GSE99374/run_limma_gse99374.R
```

**2. 免疫浸润**（需要有 R + IOBR）:
```bash
Rscript generate_all_plots.R    # 重新生成全部图表
```

**3. 数据库导入 + Web 服务**:
```bash
cd web
mysql -u root -p < sql/init.sql
npm run import:data
npm start
```

---

## 数据来源 | Data Sources

| 数据集 | GEO ID | 平台 | 组织 |
|--------|--------|------|------|
| GSE148036 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148036) | GPL21290 Illumina HiSeq 3000 | Lung tissue |
| GSE158767 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE158767) | Illumina HiSeq | Lung tissue |
| GSE99374 | [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE99374) | GPL16791 Illumina HiSeq 2500 | PBMC |

---

## 版本历史 | Changelog

| 版本 | 更新 |
|------|------|
| v2.0 | ✅ 差异分析改用 **voom/trend**, 免疫浸润改用 **TPM 转换**, 图表全部重生成, 文档补全 |
| v1.0 | 初始: `log2(counts+1)` limma, `deconvo_tme` with raw counts |
