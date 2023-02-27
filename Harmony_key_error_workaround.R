library(harmony)

RunHarmony.Seurat_CM <- function(
  object,
  group.by.vars,
  reduction = 'pca',
  dims.use = NULL,
  theta = NULL,
  lambda = NULL,
  sigma = 0.1,
  nclust = NULL,
  tau = 0,
  block.size = 0.05,
  max.iter.harmony = 10,
  max.iter.cluster = 20,
  epsilon.cluster = 1e-5,
  epsilon.harmony = 1e-4,
  plot_convergence = FALSE,
  verbose = TRUE,
  reference_values = NULL,
  reduction.save = "harmony",
  assay.use = NULL,
  project.dim = TRUE,
  ...
) {
  if (!requireNamespace('Seurat', quietly = TRUE)) {
    stop("Running Harmony on a Seurat object requires Seurat")
  }
  assay.use <- assay.use %||% Seurat::DefaultAssay(object)
  if (reduction == "pca" && !reduction %in% Seurat::Reductions(object = object)) {
    if (isTRUE(x = verbose)) {
      message("Harmony needs PCA. Trying to run PCA now.")
    }
    object <- tryCatch(
      expr = Seurat::RunPCA(
        object = object,
        assay = assay.use,
        verbose = verbose,
        reduction.name = reduction
      ),
      error = function(...) {
        stop("Harmony needs PCA. Tried to run PCA and failed.")
      }
    )
  }
  if (!reduction %in% Seurat::Reductions(object = object)) {
    stop("Requested dimension reduction is not present in the Seurat object")
  }
  embedding <- Seurat::Embeddings(object, reduction = reduction)
  if (is.null(dims.use)) {
    dims.use <- seq_len(ncol(embedding))
  }
  dims_avail <- seq_len(ncol(embedding))
  if (!all(dims.use %in% dims_avail)) {
    stop("trying to use more dimensions than computed. Rereun dimension reduction
         with more dimensions or run Harmony with fewer dimensions")
  }
  if (length(dims.use) == 1) {
    stop("only specified one dimension in dims.use")
  }
  metavars_df <- Seurat::FetchData(
    object,
    group.by.vars,
    cells = Seurat::Cells(x = object[[reduction]])
  )

  harmonyEmbed <- HarmonyMatrix(
    embedding,
    metavars_df,
    group.by.vars,
    FALSE,
    0,
    theta,
    lambda,
    sigma,
    nclust,
    tau,
    block.size,
    max.iter.harmony,
    max.iter.cluster,
    epsilon.cluster,
    epsilon.harmony,
    plot_convergence,
    FALSE,
    verbose,
    reference_values
  )

  #reduction.key <- Seurat::Key(reduction.save, quiet = TRUE)
  reduction.key <- "harmony_"
  rownames(harmonyEmbed) <- rownames(embedding)
  colnames(harmonyEmbed) <- paste0(reduction.key, seq_len(ncol(harmonyEmbed)))

  object[[reduction.save]] <- Seurat::CreateDimReducObject(
    embeddings = harmonyEmbed,
    stdev = as.numeric(apply(harmonyEmbed, 2, stats::sd)),
    assay = Seurat::DefaultAssay(object = object[[reduction]]),
    key = reduction.key
  )
  if (project.dim) {
    object <- Seurat::ProjectDim(
      object,
      reduction = reduction.save,
      overwrite = TRUE,
      verbose = FALSE
    )
  }
  return(object)
}
