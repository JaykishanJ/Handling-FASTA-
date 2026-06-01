#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

cran_packages <- c("data.table", "RobustRankAggreg")
bioc_packages <- c("GEOquery", "limma", "WGCNA")

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (bioc) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
          install.packages("BiocManager", repos = "https://cloud.r-project.org")
        }
        BiocManager::install(pkg, ask = FALSE, update = FALSE)
      } else {
        install.packages(pkg, repos = "https://cloud.r-project.org")
      }
    }
  }
}

install_if_missing(cran_packages, bioc = FALSE)
install_if_missing(bioc_packages, bioc = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(limma)
  library(RobustRankAggreg)
  library(WGCNA)
})

WGCNA::allowWGCNAThreads()

make_dataset_key <- function(gse, gpl) {
  paste(gse, gpl, sep = "_")
}

output_dir <- "outputs"
deg_dir <- file.path(output_dir, "deg")
sample_dir <- file.path(output_dir, "sample_annotations")
rra_dir <- file.path(output_dir, "rra")
wgcna_dir <- file.path(output_dir, "wgcna")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(rra_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(wgcna_dir, showWarnings = FALSE, recursive = TRUE)

datasets <- data.frame(
  gse = c(
    "GSE36376", "GSE39791", "GSE45114", "GSE57957", "GSE60502",
    "GSE76297", "GSE76427", "GSE84005", "GSE14520", "GSE14520"
  ),
  gpl = c(
    "GPL10558", "GPL10558", "GPL5918", "GPL10558", "GPL96",
    "GPL17586", "GPL10558", "GPL5175", "GPL571", "GPL3921"
  ),
  stringsAsFactors = FALSE
)
datasets$key <- make_dataset_key(datasets$gse, datasets$gpl)

strict_grouping <- FALSE
get_gpl <- TRUE
rra_top_n <- 200
log2_transform_threshold <- 50
wgcna_dataset_key <- "GSE14520_GPL571"
wgcna_max_genes <- 5000
wgcna_min_module_size <- 30
wgcna_merge_cut_height <- 0.25
wgcna_default_soft_power <- 6
wgcna_top_hub_genes <- 50
get_symbol_column_name <- function(fdata) {
  candidates <- c(
    "Gene Symbol", "Gene symbol", "GENE_SYMBOL", "Symbol", "SYMBOL",
    "gene_assignment", "Gene Assignment", "GENE", "GENES", "ILMN_Gene"
  )
  cols_lower <- tolower(colnames(fdata))
  idx <- match(tolower(candidates), cols_lower)
  idx <- idx[!is.na(idx)][1]
  if (is.na(idx)) {
    return(NA_character_)
  }
  colnames(fdata)[idx]
}

clean_symbols <- function(symbols) {
  symbols <- as.character(symbols)
  symbols <- gsub("\\s*///\\s*", ";", symbols)
  symbols <- gsub("\\s*//\\s*", ";", symbols)
  symbols <- gsub("\\s*\\|\\s*", ";", symbols)
  symbols <- gsub("\\s*,\\s*", ";", symbols)
  tokens <- strsplit(symbols, ";")
  cleaned <- vapply(tokens, function(x) {
    x <- trimws(x)
    x <- x[nzchar(x)]
    if (length(x) == 0) {
      return(NA_character_)
    }
    x[1]
  }, character(1))
  cleaned
}

assign_tumor_normal_groups <- function(pheno, key) {
  combined <- apply(pheno, 1, function(x) paste(x, collapse = " ; "))
  combined_lower <- tolower(combined)
  combined_lower <- gsub("tumour", "tumor", combined_lower)
  is_normal <- grepl("normal|adjacent|non[\\s-]?tumor|control|healthy", combined_lower)
  is_tumor <- grepl("tumor|carcinoma|hcc|hepatocellular", combined_lower)
  group <- ifelse(is_tumor & !is_normal, "Tumor",
    ifelse(is_normal, "Normal", NA_character_)
  )
  if (any(is.na(group))) {
    out_path <- file.path(sample_dir, paste0(key, "_pheno.csv"))
    write.csv(
      data.frame(sample = rownames(pheno), group = group, pheno = combined),
      out_path,
      row.names = FALSE
    )
    if (strict_grouping) {
      stop("Missing sample groups for ", key, ". Fill in: ", out_path)
    }
    message("Skipping ", key, " due to missing sample groups. See: ", out_path)
    return(NULL)
  }
  factor(group, levels = c("Normal", "Tumor"))
}

prepare_eset <- function(eset, key) {
  expr <- exprs(eset)
  if (max(expr, na.rm = TRUE) > log2_transform_threshold) {
    expr <- log2(expr + 1)
  }

  fdata <- fData(eset)
  symbol_col <- get_symbol_column_name(fdata)
  if (is.na(symbol_col)) {
    stop("No gene symbol column found for ", key)
  }
  symbols <- clean_symbols(fdata[[symbol_col]])
  keep <- !is.na(symbols) & symbols != ""
  expr <- expr[keep, , drop = FALSE]
  symbols <- symbols[keep]
  expr <- limma::avereps(expr, ID = symbols)

  pheno <- pData(eset)
  group <- assign_tumor_normal_groups(pheno, key)
  if (is.null(group)) {
    return(NULL)
  }
  list(expr = expr, group = group, pheno = pheno)
}

run_limma <- function(expr, group) {
  design <- model.matrix(~0 + group)
  colnames(design) <- levels(group)
  fit <- lmFit(expr, design)
  contrast <- makeContrasts(Tumor - Normal, levels = design)
  fit2 <- eBayes(contrasts.fit(fit, contrast))
  res <- topTable(fit2, number = Inf, sort.by = "P")
  res$Gene <- rownames(res)
  res <- res[, c("Gene", setdiff(colnames(res), "Gene"))]
  res
}

extract_ranked_up_down_genes <- function(res) {
  res <- res[!is.na(res$logFC) & !is.na(res$P.Value), ]
  up <- res[res$logFC > 0, ]
  down <- res[res$logFC < 0, ]
  up <- up[order(up$P.Value, -up$logFC), ]
  down <- down[order(down$P.Value, down$logFC), ]
  list(
    up = unique(up$Gene),
    down = unique(down$Gene)
  )
}

get_gene_names <- function(rra) {
  if (is.null(rra)) {
    return(character())
  }
  if ("Name" %in% names(rra)) {
    return(rra$Name)
  }
  if ("name" %in% names(rra)) {
    return(rra$name)
  }
  warning(
    sprintf(
      "Unexpected RRA column names (found: %s); using first column as gene names.",
      paste(names(rra), collapse = ", ")
    )
  )
  rra[[1]]
}

fetch_expression_sets <- function(gse) {
  geo <- getGEO(gse, GSEMatrix = TRUE, getGPL = get_gpl)
  if (inherits(geo, "ExpressionSet")) {
    return(list(geo))
  }
  geo
}

deg_results <- list()
ranked_gene_lists_upregulated <- list()
ranked_gene_lists_downregulated <- list()
wgcna_prepared <- NULL

for (i in seq_len(nrow(datasets))) {
  gse <- datasets$gse[i]
  gpl <- datasets$gpl[i]
  key <- datasets$key[i]
  message("Processing ", key)
  esets <- fetch_expression_sets(gse)
  eset <- NULL
  for (candidate in esets) {
    if (annotation(candidate) == gpl) {
      eset <- candidate
      break
    }
  }
  if (is.null(eset)) {
    message("No platform match for ", key)
    next
  }

  prepared <- prepare_eset(eset, key)
  if (is.null(prepared)) {
    next
  }
  if (key == wgcna_dataset_key) {
    wgcna_prepared <- prepared
  }
  deg <- run_limma(prepared$expr, prepared$group)
  deg_results[[key]] <- deg
  fwrite(deg, file = file.path(deg_dir, paste0(key, "_limma.csv")))

  ranks <- extract_ranked_up_down_genes(deg)
  ranked_gene_lists_upregulated[[key]] <- ranks$up
  ranked_gene_lists_downregulated[[key]] <- ranks$down
}

if (length(deg_results) == 0) {
  stop("No differential expression results produced.")
}

if (length(ranked_gene_lists_upregulated) > 1) {
  rra_up <- aggregateRanks(ranked_gene_lists_upregulated)
  fwrite(rra_up, file = file.path(rra_dir, "rra_up.csv"))
} else {
  rra_up <- NULL
}

if (length(ranked_gene_lists_downregulated) > 1) {
  rra_down <- aggregateRanks(ranked_gene_lists_downregulated)
  fwrite(rra_down, file = file.path(rra_dir, "rra_down.csv"))
} else {
  rra_down <- NULL
}

if (!is.null(rra_up)) {
  rra_up_top <- head(get_gene_names(rra_up), rra_top_n)
} else if (length(ranked_gene_lists_upregulated) > 0) {
  rra_up_top <- unique(unlist(lapply(ranked_gene_lists_upregulated, head, rra_top_n)))
} else {
  rra_up_top <- character()
}

if (!is.null(rra_down)) {
  rra_down_top <- head(get_gene_names(rra_down), rra_top_n)
} else if (length(ranked_gene_lists_downregulated) > 0) {
  rra_down_top <- unique(unlist(lapply(ranked_gene_lists_downregulated, head, rra_top_n)))
} else {
  rra_down_top <- character()
}

rra_union_top <- unique(c(rra_up_top, rra_down_top))

if (is.null(wgcna_prepared)) {
  stop("WGCNA dataset not available: ", wgcna_dataset_key)
}

message("Running WGCNA on ", wgcna_dataset_key)
wgcna_data <- wgcna_prepared
expr <- wgcna_data$expr
group <- wgcna_data$group

if (nrow(expr) > wgcna_max_genes) {
  variances <- apply(expr, 1, var)
  keep <- order(variances, decreasing = TRUE)[seq_len(wgcna_max_genes)]
  expr <- expr[keep, ]
}

datExpr <- t(expr)
wgcna_qc_results <- goodSamplesGenes(datExpr, verbose = 0)
if (!wgcna_qc_results$allOK) {
  removed_samples <- sum(!wgcna_qc_results$goodSamples)
  removed_genes <- sum(!wgcna_qc_results$goodGenes)
  message("WGCNA QC removed ", removed_samples, " samples and ", removed_genes, " genes.")
  datExpr <- datExpr[wgcna_qc_results$goodSamples, wgcna_qc_results$goodGenes]
}

powers <- 1:20
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0)
soft_power <- sft$powerEstimate
if (is.na(soft_power)) {
  soft_power <- wgcna_default_soft_power
  message("Using default WGCNA soft power (", soft_power, ") for a signed network.")
}

net <- blockwiseModules(
  datExpr,
  power = soft_power,
  TOMType = "signed",
  minModuleSize = wgcna_min_module_size,
  mergeCutHeight = wgcna_merge_cut_height,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  verbose = 0
)

module_colors <- labels2colors(net$colors)
MEs <- net$MEs
tumor_status <- data.frame(Tumor = as.numeric(group == "Tumor"))
module_trait_cor <- cor(MEs, tumor_status, use = "p")
module_trait_p <- corPvalueStudent(module_trait_cor, nrow(datExpr))

module_names <- rownames(module_trait_cor)
if (all(startsWith(module_names, "ME"))) {
  module_names <- sub("^ME", "", module_names)
}
module_trait <- data.frame(
  module = module_names,
  correlation = module_trait_cor[, 1],
  p_value = module_trait_p[, 1]
)
fwrite(module_trait, file = file.path(wgcna_dir, paste0(wgcna_dataset_key, "_module_trait.csv")))

target_idx <- which.max(abs(module_trait$correlation))
target_module <- module_trait$module[target_idx]
target_color <- target_module

module_genes <- names(module_colors)[module_colors == target_color]
if (length(module_genes) == 0) {
  stop("No genes found for module ", target_color)
}

kme <- cor(datExpr, MEs, use = "p")
target_me <- paste0("ME", target_color)
kme_target <- kme[, target_me]
gene_significance <- cor(datExpr, tumor_status$Tumor, use = "p")

module_summary <- data.frame(
  gene = module_genes,
  kME = kme_target[module_genes],
  geneSignificance = gene_significance[module_genes]
)
module_summary <- module_summary[order(-abs(module_summary$kME)), ]
fwrite(module_summary, file = file.path(wgcna_dir, paste0(wgcna_dataset_key, "_", target_color, "_module_genes.csv")))

hub_genes <- head(module_summary, wgcna_top_hub_genes)
fwrite(hub_genes, file = file.path(wgcna_dir, paste0(wgcna_dataset_key, "_", target_color, "_hub_genes.csv")))

candidate_hubs <- intersect(module_genes, rra_union_top)
fwrite(
  data.frame(gene = candidate_hubs),
  file = file.path(wgcna_dir, paste0(wgcna_dataset_key, "_", target_color, "_rra_intersect.csv"))
)

message("Pipeline complete. Outputs saved in: ", output_dir)
