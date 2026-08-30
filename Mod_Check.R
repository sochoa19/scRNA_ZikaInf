#Single cell module performance
library(Seurat)
library(hdWGCNA)
library(fs)
library(ggplot2)

set.seed(1234)
B <- 300

Target_name<-"new_UCC_int"
Attempt_number<-"Attempt_1"

dataout<-paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/Module_Covariance/")
dir_create(dataout)

obj <- readRDS(paste0("/data/scRNA/HMC3_ZSC/Seurat_OUT/",Target_name,"/PIC_Only_WGCNA/",Attempt_number,"/",Target_name,"WGCNA_final_PIC_WT_projected_allPIC.rds" ))
expr <- GetAssayData(obj, assay = "RNA", layer = "data")
meta <- obj[[]]

# Example: module_genes is a named list of gene vectors
module_genes<-GetModules(obj)# module_genes <- biclustMembership
mod_list<-list()
for (i in setdiff(levels(factor(module_genes$module)),"grey")){
  mod_list[[i]]<-unlist(module_genes[module_genes$module==i,1])
}

module_genes<-mod_list

eligible_genes <- intersect(rownames(expr), unique(unlist(module_genes)))

module_variance_test <- function(module_genes, expr, meta,
                                 sample_col = "Sample",
                                 n_perm = 1000) {
  module_genes <- intersect(unique(module_genes), rownames(expr))
  background <- rownames(expr)[Matrix::rowSums(expr != 0) > 0]
  
  if (length(module_genes) < 2) return(NULL)
  
  results <- lapply(unique(meta[[sample_col]]), function(sample_id) {
    cells <- rownames(meta)[meta[[sample_col]] == sample_id]
    cells <- intersect(cells, colnames(expr))
    
    if (length(cells) < 3) return(NULL)
    
    # Gene-level variance across cells within this sample
    gene_var <- apply(expr[, cells, drop = FALSE], 1, var, na.rm = TRUE)
    
    observed <- mean(gene_var[module_genes], na.rm = TRUE)
    
    #Take a sampling of random genes from any one cell, (the number of genes samples is the size of the module being tested)
    #then find the mean variance between those genes. repeats this for n_perms times.
    null_values <- replicate(n_perm, {
      #take a random sampling of genes from expr
      random_genes <- sample(background, length(module_genes))
      
      mean(gene_var[random_genes], na.rm = TRUE)
    })
    
    data.frame(
      sample = sample_id,
      n_cells = length(cells),
      n_module_genes = length(module_genes),
      observed_variance = observed,
      null_mean = mean(null_values),
      null_sd = sd(null_values),
      variance_z = ( mean(null_values) - observed) / sd(null_values),
      empirical_p = (1 + sum(null_values <= observed)) / (n_perm + 1)
    )
  })
  
  do.call(rbind, results)
}

variance_results <- do.call(
  rbind,
  lapply(names(module_genes), function(module_name) {
    result <- module_variance_test(module_genes[[module_name]], expr, meta)
    if (!is.null(result)) {
      result$module <- module_name
    }
    print(module_name)
    result
  })
)

# Mean pairwise covariance of genes across cells in one sample.
# The diagonal is excluded because it is gene variance, not coexpression.
mean_pairwise_covariance <- function(gene_names, expr, cells, min_cells = 3) {
  gene_names <- intersect(unique(gene_names), rownames(expr))
  cells <- intersect(cells, colnames(expr))

  if (length(gene_names) < 2 || length(cells) < min_cells) return(NA_real_)

  expression_matrix <- as.matrix(expr[gene_names, cells, drop = FALSE])
  expression_matrix <- expression_matrix[
    apply(expression_matrix, 1, function(values) sum(is.finite(values)) >= min_cells),
    , drop = FALSE]
  if (nrow(expression_matrix) < 2) return(NA_real_)

  covariance_matrix <- stats::cov(t(expression_matrix),
                           use = "pairwise.complete.obs")
  covariance_matrix[lower.tri(covariance_matrix, diag = TRUE)] <- NA_real_
  covariance_values <- covariance_matrix[upper.tri(covariance_matrix)]
  if (!any(is.finite(covariance_values))) return(NA_real_)
  mean(covariance_values, na.rm = TRUE)
}

module_covariance_test <- function(module_genes, expr, meta,
                                   sample_col = "Sample",
                                   background_col = "Background",
                                   reference_background = "C",
                                   n_perm = 1000,
                                   sample_size = NULL,
                                   min_cells = 3) {
  module_genes <- intersect(unique(module_genes), rownames(expr))
  background_genes <- rownames(expr)[Matrix::rowSums(expr != 0) > 0]

  if (length(module_genes) < 2) return(NULL)
  if (is.null(sample_size)) sample_size <- length(module_genes)
  sample_size <- min(sample_size, length(module_genes), length(background_genes))
  if (sample_size < 2) return(NULL)

  sample_results <- lapply(unique(meta[[sample_col]]), function(sample_id) {
    cells <- rownames(meta)[meta[[sample_col]] == sample_id]
    cells <- intersect(cells, colnames(expr))
    if (length(cells) < min_cells) return(NULL)

    observed_covariance <- mean_pairwise_covariance(module_genes, expr, cells,
                                                  min_cells = min_cells)
    if (!is.finite(observed_covariance)) return(NULL)

    # Match the module and null gene-set sizes on every iteration.
    module_check <- replicate(n_perm, {
      selected_genes <- sample(module_genes, sample_size)
      mean_pairwise_covariance(selected_genes, expr, cells,
                              min_cells = min_cells)
    })
    Null_check <- replicate(n_perm, {
      selected_genes <- sample(background_genes, sample_size)
      mean_pairwise_covariance(selected_genes, expr, cells,
                              min_cells = min_cells)
    })

    module_mean <- mean(module_check, na.rm = TRUE)
    null_mean <- mean(Null_check, na.rm = TRUE)
    null_sd <- sd(Null_check, na.rm = TRUE)
    module_p <- (1 + sum(Null_check >= observed_covariance, na.rm = TRUE)) /
      (1 + sum(!is.na(Null_check)))

    data.frame(
      sample = sample_id,
      background = paste(unique(meta[rownames(meta) %in% cells, background_col]),
             collapse = ";"),
      n_cells = length(cells),
      n_module_genes = length(module_genes),
      module_mean_covariance = observed_covariance,
      null_mean_covariance = null_mean,
      null_sd_covariance = null_sd,
      covariance_z = (observed_covariance - null_mean) / null_sd,
      empirical_p_greater_than_null = module_p,
      Null_check = I(list(Null_check)),
      Module_check = I(list(module_check))
    )
  })

  do.call(rbind, sample_results)
}

covariance_results <- do.call(
  rbind,
  lapply(names(module_genes), function(module_name) {
    result <- module_covariance_test(module_genes[[module_name]], expr, meta,
                                     n_perm = B,sample_size = 100)
    if (!is.null(result)) result$module <- module_name
    print(module_name)
    result
  })
)

saveRDS(covariance_results,paste0(dataout,"Covariance_df.rds"))
# Compare module covariance against C for each module. These are comparisons
# of sample-level means, rather than treating permutation draws as replicates.
reference_background <- "C"
c_reference <- covariance_results[
  covariance_results$background == reference_background, , drop = FALSE]
covariance_results$c_vs_reference_p <- NA_real_
for (module_name in unique(covariance_results$module)) {
  reference_rows <- which(c_reference$module == module_name)
  reference_values <- unlist(c_reference$Module_check[reference_rows])
  query_rows <- which(covariance_results$module == module_name &
                      covariance_results$background != reference_background)
  if (length(reference_values) > 0 && length(query_rows) > 0) {
    covariance_results$c_vs_reference_p[query_rows] <- vapply(
      query_rows,
      function(row) t.test(unlist(covariance_results$Module_check[[row]]),
                            y = reference_values)$p.value,
      numeric(1)
    )
  }
}

write.csv(covariance_results,paste0(dataout,"module_covariance_results.csv") , row.names = FALSE)

plot_module_covariance <- function(result, output_dir = "module_covariance_plots") {
  output_dir<-paste0(dataout,output_dir)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  for (row in seq_len(nrow(result))) {
    null_values <- unlist(result$Null_check[[row]])
    module_mean <- result$module_mean_covariance[row]
    png(file.path(output_dir, paste0(result$module[row], "_",
                                     result$sample[row], ".png")))
    hist(null_values, breaks = 30,
         main = paste(result$module[row], result$sample[row]),
         xlab = "Mean pairwise gene covariance",
         xlim = c(0,module_mean+module_mean/100),
         log2="x")
    abline(v = module_mean, col = "red", lwd = 2)
    dev.off()
  }
}

plot_module_covariance(covariance_results)

# Compare module coexpression distributions across backgrounds. With one
# biological sample per background, these violins describe permutation draws.
plot_background_coexpression <- function(result,
                                          output_file = "module_coexpression_by_background.png") {
  plot_values <- do.call(rbind, lapply(seq_len(nrow(result)), function(row) {
    data.frame(
      module = result$module[row],
      sample = result$sample[row],
      background = result$background[row],
      covariance = unlist(result$Module_check[[row]])
    )
  }))

  observed_values <- unique(result[c("module", "sample", "background",
                                     "module_mean_covariance")])

  p <- ggplot(plot_values, aes(x = background, y = covariance,
                               fill = background)) +
    geom_violin(trim = FALSE, scale = "width") +
    geom_point(data = observed_values,
               aes(x = background, y = module_mean_covariance),
               color = "black", size = 1.5, inherit.aes = FALSE,
               position = position_jitter(width = 0.08)) +
    facet_wrap(~ module, scales = "free_y") +
    labs(x = "Background", y = "Mean pairwise gene covariance") +
    theme_classic() +
    theme(legend.position = "none")

  ggsave(paste0(dataout,output_file), p, width = 12, height = 8, dpi = 300)
  p
}

background_coexpression_plot <- plot_background_coexpression(covariance_results)

# Exploratory metacell-level analysis.
# Instead of testing on the pooled metacell expression, subset to the original
# cells merged into each metacell and treat that cell subset as its own
# pseudo-sample. Then summarize the metacell-level results back to the larger
# biological sample for interpretation.
metacell_obj <- MetacellsByGroups(
  seurat_obj = obj,
  group.by = "Sample",
  reduction = "pca",
  k = 60,
  max_shared = 15,
  ident.group = "Sample"
)
metacell_obj <- NormalizeMetacells(metacell_obj)
mc <- GetMetacellObject(metacell_obj)
mc_meta <- mc[[]]

if (!"cells_merged" %in% colnames(mc_meta)) {
  stop("Metacell metadata must contain cells_merged")
}

mc_meta$metacell_id <- rownames(mc_meta)
mc_meta$cell_count <- lengths(strsplit(mc_meta$cells_merged, split = ","))
mc_meta$Background <- gsub("ZSC|P", "", mc_meta$Sample)

####RUN METACELL COVARIANCE TEST
metacell_covariance_results <- do.call(
  rbind,
  lapply(seq_len(nrow(mc_meta)), function(i) {
    #START BY SUBSETTING THE ENTIRE SEURAT OBJECT WITH THE CELLS INSIDE METACELL i
    cells_in_metacell <- strsplit(mc_meta$cells_merged[i], split = ",")[[1]]
    cells_in_metacell <- intersect(cells_in_metacell, colnames(expr))
    if (length(cells_in_metacell) < 2) return(NULL)
    message("Working on metacell #",i, " out of ",nrow(mc_meta))
    
    metacell_expr <- expr[, cells_in_metacell, drop = FALSE]
    metacell_meta <- meta[rownames(meta) %in% cells_in_metacell, , drop = FALSE]
    metacell_meta$Sample <- paste0(metacell_meta$Sample, "_", i)
    metacell_meta$Background <- unique(metacell_meta$Background)
####NOW ITERATE THE COVARIANCE TEST THROUGH EACH MODULE
    do.call(rbind, lapply(names(module_genes), function(module_name) {
      module_gene_vec <- module_genes[[module_name]]
      if (length(module_gene_vec) < 2) return(NULL)

      result <- module_covariance_test(
        module_genes = module_gene_vec,
        expr = metacell_expr,
        meta = metacell_meta,
        sample_col = "Sample",
        background_col = "Background",
        n_perm = B,
        sample_size = min(100, length(module_gene_vec))
      )

      if (is.null(result)) return(NULL)

      result$module <- module_name
      result$metacell_id <- mc_meta$metacell_id[i]
      result$sample_group <- as.character(mc_meta$Sample[i])
      result$sample <- result$sample_group
      result
    }))
  })
)

if (!is.null(metacell_covariance_results) && nrow(metacell_covariance_results) > 0) {
  metacell_sample_summary <- do.call(
    rbind,
    lapply(unique(metacell_covariance_results$module), function(module_name) {
      module_rows <- metacell_covariance_results[
        metacell_covariance_results$module == module_name,
        , drop = FALSE]
      if (nrow(module_rows) == 0) return(NULL)

      do.call(rbind, lapply(unique(module_rows$sample_group), function(sample_group_name) {
        rows <- module_rows[module_rows$sample_group == sample_group_name, , drop = FALSE]
        data.frame(
          module = module_name,
          sample = sample_group_name,
          background = paste(unique(rows$background), collapse = ";"),
          n_metacells = nrow(rows),
          module_mean_covariance = median(rows$module_mean_covariance, na.rm = TRUE),
          null_mean_covariance = median(rows$null_mean_covariance, na.rm = TRUE),
          null_sd_covariance = median(rows$null_sd_covariance, na.rm = TRUE),
          covariance_z = median(rows$covariance_z, na.rm = TRUE),
          empirical_p_greater_than_null = min(rows$empirical_p_greater_than_null, na.rm = TRUE),
          metacell_covariance_distribution = I(list(rows$module_mean_covariance)),
          metacell_ids = I(list(rows$metacell_id))
        )
      }))
    })
  )
} else {
  metacell_sample_summary <- NULL
}

write.csv(metacell_covariance_results,
          paste0(dataout, "metacell_covariance_results.csv"),
          row.names = FALSE)
if (!is.null(metacell_sample_summary)) {
  write.csv(metacell_sample_summary,
            paste0(dataout, "metacell_covariance_summary.csv"),
            row.names = FALSE)
}
plot_module_covariance(metacell_covariance_results,
                       output_dir = "metacell_covariance_plots")
plot_background_coexpression(
  metacell_covariance_results,
  output_file = "metacell_coexpression_by_background.png"
)





