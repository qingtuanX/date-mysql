# Volcano Web (MySQL)

该子项目实现了你需要的网站能力：

- 从 MySQL 读取差异分析结果；
- 展示火山图（X=`log2FC`, Y=`-log10(FDR)`）；
- **悬浮或点击数据库表格中的值时，动态高亮火山图对应点并显示坐标**；
- 在图上点击点时，反向联动到坐标信息面板。

## 1) 初始化数据库

先在 MySQL 中执行：

`web/sql/init.sql`

它会创建：

- 数据库：`date_mysql`
- 表：`limma_results`

## 2) 配置环境变量

复制：

`web/.env.example` -> `web/.env`

按你的 MySQL 实际账号修改。

## 3) 安装依赖

在 `web` 目录执行：

```bash
npm install
```

## 4) 导入现有 limma 结果到 MySQL

```bash
npm run import:data
```

脚本会自动扫描 `mRNA/GSE*/` 下形如：

`*_limma_test_result.*.txt`

并写入 MySQL。

## 5) 启动网站

```bash
npm start
```

打开：

`http://localhost:3000`

## 交互说明

- 表格行 **hover**：临时高亮图中点并显示坐标
- 表格行 **click**：固定该点，显示固定状态
- 图中点 **click**：同样固定并刷新坐标信息
- 点击“取消固定点”：恢复 hover 联动模式
