# date-mysql — mRNA 转录组差异分析与交互可视化平台

> mRNA Transcriptomics Differential Expression Analysis & Interactive Visualization Platform

---

## 📖 项目简介 | Overview

**date-mysql** 是一个面向生物信息学研究的 mRNA 转录组差异表达分析（DEA）全流程工具，整合了 R 语言统计分析、Python/Node.js 数据处理以及 Web 交互可视化。项目以 GEO 公共数据集为核心，使用 `limma` 进行多组差异分析，并通过 MySQL 数据库驱动的 Web 应用实现火山图的动态交互浏览。

**date-mysql** is a full-stack bioinformatics pipeline for mRNA transcriptomics differential expression analysis (DEA). It integrates R-based statistical analysis, Node.js data processing, and an interactive web visualization platform. Using public GEO datasets, it performs multi-group differential analysis with `limma` and provides a MySQL-driven web interface for dynamic volcano plot exploration.

---

## 🏗 项目结构 | Project Structure

```
date-mysql/
├── mRNA/                          # 核心数据分析模块 | Core analysis module
│   ├── GSE148036/                 # 结核病/肺癌数据集 | Tuberculosis/Lung Cancer
│   ├── GSE158767/                 # 黑色素瘤 MITF 数据集 | Melanoma MITF
│   ├── GSE99374/                  # 潜伏结核感染数据集 | Latent TB Infection
│   ├── 【199】多组火山图/           # 多组火山图绘制模板（SCI 标准）
│   ├── 【206】火山图/               # 单组火山图绘制模板（SCI 标准）
│   └── 【378】limma多组差异分析/     # limma 多组差异分析 R Markdown 脚本
├── web/                           # Web 交互应用 | Web application
│   ├── public/                    # 前端静态资源（HTML/CSS/JS）
│   ├── scripts/                   # 数据导入脚本（limma 结果 → MySQL）
│   ├── sql/                       # 数据库初始化 SQL
│   ├── server.js                  # Express 主服务入口
│   └── db.js                      # MySQL 连接池配置
└── README.md
```

---

## 🚀 分析工作流 | Analysis Workflow

### 1. 数据获取 | Data Acquisition
从 NCBI GEO 下载 GSE series 数据集，包括临床元数据和表达矩阵（FPKM/TPM）。

### 2. 预处理 | Preprocessing
表达矩阵的标准化（Normalization）和 log 转换（Log-transformation）。

### 3. 差异分析 | Differential Analysis
使用 R 包 `limma` 计算 Fold Change 和显著性 p-value，支持多组比较（如 `immune vs Others`、`MITF-low vs Others`）。

### 4. 免疫浸润分析 | Immune Infiltration
整合 CIBERSORT、EPIC、ESTIMATE、MCP-counter 等算法评估肿瘤微环境免疫细胞组分。

### 5. 可视化 | Visualization
- **火山图（Volcano Plot）**：log2FC → -log10(adj.P.Val)
- **热图（Heatmap）**：差异表达基因聚类
- **Web 交互面板**：MySQL 驱动，支持表格 ↔ 图表双向联动

---

## 🛠 技术栈 | Tech Stack

| 层级 | 技术 |
|------|------|
| 统计分析 | R (≥4.0), `limma`, `ggplot2`, `dplyr`, `Biobase` |
| 后端 | Node.js, Express, MySQL2 |
| 前端 | 原生 HTML/CSS/JS, Canvas API |
| 数据库 | MySQL 8.0 |

---

## 🖥 Web 应用快速启动 | Web App Quick Start

### 前置要求 | Prerequisites

- **Node.js** ≥ 16
- **MySQL** ≥ 8.0
- **R** ≥ 4.0（仅数据分析需）| (analysis only)

### 1. 初始化数据库 | Initialize Database

在 MySQL 中执行 `web/sql/init.sql`，创建 `date_mysql` 数据库和 `limma_results` 表。

```bash
mysql -u root -p < web/sql/init.sql
```

### 2. 配置环境变量 | Configure Environment

```bash
cp web/.env.example web/.env
# 编辑 web/.env，填入 MySQL 连接信息
```

### 3. 安装依赖 | Install Dependencies

```bash
cd web
npm install
```

### 4. 导入数据 | Import Data

将 `mRNA/GSE*/` 下的 limma 分析结果批量导入 MySQL：

```bash
npm run import:data
```

### 5. 启动服务 | Start Server

```bash
npm start
```

浏览器打开 `http://localhost:3000`。

---

## 🖱 Web 交互说明 | Interaction Guide

| 操作 | 效果 |
|------|------|
| **表格行悬停** | 火山图中临时高亮对应点并显示坐标 |
| **表格行点击** | 固定选中点，显示详细信息 |
| **图中点点击** | 反向联动，更新坐标信息面板 |
| **取消固定** | 恢复 hover 联动模式 |

---

## 📊 数据集说明 | Datasets

| 数据集 | 来源 | 样本数 | 研究方向 |
|--------|------|--------|----------|
| GSE148036 | GEO | 116 | 结核病 vs 腺癌 vs 结节病 vs 正常 |
| GSE158767 | GEO | 47 | 黑色素瘤 MITF 高低表达分组 |
| GSE99374 | GEO | 179 | 潜伏结核感染 vs 健康对照 |

---

## ⚠️ 注意事项 | Notes

- 部分 `.zip` / `.gz` 文件为原始数据和注释文件的压缩包，超过 50MB 建议使用 Git LFS 管理。
- `web/scripts/import-limma.js` 需在配置好数据库后运行。
- `.env` 文件包含敏感数据库密码，已加入 `.gitignore`，请勿提交。

---

## 📄 License

This project is for academic research purposes.
