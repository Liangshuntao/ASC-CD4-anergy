# Import reviewed labels only; this function performs no cell-type inference.
map_reviewed_labels <- function(meta, mapping) {
  required <- c("cluster", "celltype_v2", "major_celltype_v2", "n_cells")
  if (!"seurat_clusters" %in% names(meta) || !all(required %in% names(mapping))) stop("Missing reviewed annotation fields")
  clusters <- as.character(meta$seurat_clusters)
  keys <- as.character(mapping$cluster)
  if (anyNA(clusters) || anyNA(keys) || anyDuplicated(keys) || !all(clusters %in% keys)) stop("Invalid or incomplete cluster map")
  n <- as.integer(table(factor(clusters, levels = keys)))
  if (anyNA(mapping$n_cells) || !identical(n, as.integer(mapping$n_cells))) stop("Reviewed map cell counts do not match checkpoint")
  idx <- match(clusters, keys)
  labels <- mapping[idx, c("celltype_v2", "major_celltype_v2"), drop = FALSE]
  if (anyNA(labels) || any(!nzchar(as.matrix(labels)))) stop("Missing reviewed labels")
  rownames(labels) <- rownames(meta)
  labels
}
