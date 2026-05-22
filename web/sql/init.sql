CREATE DATABASE IF NOT EXISTS date_mysql
  DEFAULT CHARACTER SET utf8mb4
  DEFAULT COLLATE utf8mb4_unicode_ci;

USE date_mysql;

CREATE TABLE IF NOT EXISTS limma_results (
  id BIGINT UNSIGNED NOT NULL AUTO_INCREMENT,
  dataset VARCHAR(64) NOT NULL,
  comparison VARCHAR(128) NOT NULL,
  gene VARCHAR(255) NOT NULL,
  log2fc DOUBLE NOT NULL,
  t_value DOUBLE NULL,
  b_value DOUBLE NULL,
  p_value DOUBLE NULL,
  fdr DOUBLE NOT NULL,
  neg_log10_fdr DOUBLE NOT NULL,
  created_at TIMESTAMP NOT NULL DEFAULT CURRENT_TIMESTAMP,
  PRIMARY KEY (id),
  KEY idx_dataset_comparison (dataset, comparison),
  KEY idx_gene (gene),
  KEY idx_fdr (fdr)
) ENGINE=InnoDB DEFAULT CHARSET=utf8mb4 COLLATE=utf8mb4_unicode_ci;
