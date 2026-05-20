library(readr)
library(parallel)
library(MendelianRandomization)
library(gsmr2)

args <- commandArgs(trailingOnly = TRUE)
input_path <- args[1]
output_path <- args[2]
xqtl_data <- args[3]
feature_col <- args[4]
qtl_sample_size <- as.numeric(args[5])
jobs <- as.numeric(args[6])

process_one <- function(file, xqtl_data, qtl_sample_size, feature_col) {
  file_prefix <- sub("\\.tsv$", "", file)
  file_name <- basename(file)
  parts <- strsplit(file_name, "\\.", fixed = FALSE)[[1]]
  trait_name <- if (length(parts) >= 2) parts[2] else NA_character_

  make_result_row <- function(feature, ivw_beta = NA_real_, ivw_se = NA_real_,
                              ivw_p = NA_real_, ivw_n_total = NA_integer_,
                              error = NA_character_) {
    data.frame(
      Dataset = xqtl_data,
      Protein = feature,
      Trait = trait_name,
      IVW_Beta = ivw_beta,
      IVW_SE = ivw_se,
      IVW_P = ivw_p,
      IVW_N_Total = ivw_n_total,
      Error = error,
      stringsAsFactors = FALSE
    )
  }

  snp_coeff_id <- tryCatch({
    scan(sprintf("%s.xmat.gz", file_prefix), what = "", nlines = 1, quiet = TRUE)
  }, error = function(e) character())

  snp_coeff <- tryCatch({
    read.table(sprintf("%s.xmat.gz", file_prefix), header = FALSE, skip = 2, check.names = FALSE)
  }, error = function(e) NULL)

  mr_data_all <- tryCatch({
    read_tsv(sprintf("%s.gz", file_prefix), show_col_types = FALSE, progress = FALSE)
  }, error = function(e) NULL)

  if (is.null(mr_data_all) || !feature_col %in% names(mr_data_all)) {
    return(make_result_row(feature = NA, error = "Failed to read input files"))
  }

  feature_vals <- unique(mr_data_all[[feature_col]])
  if (length(feature_vals) == 0L) {
    return(make_result_row(feature = NA, ivw_n_total = 0L,
                           error = sprintf("No values found in feature column: %s", feature_col)))
  }

  res_list <- vector("list", length(feature_vals))

  for (i in seq_along(feature_vals)) {
    val <- feature_vals[i]
    feature_idx <- if (is.na(val)) {
      is.na(mr_data_all[[feature_col]])
    } else {
      !is.na(mr_data_all[[feature_col]]) & mr_data_all[[feature_col]] == val
    }
    mr_data <- mr_data_all[feature_idx, , drop = FALSE]

    if (is.null(snp_coeff) || length(snp_coeff_id) == 0 || nrow(mr_data) == 0L) {
      res_list[[i]] <- make_result_row(feature = val, error = "Failed to read input files")
      next
    }

    snp_id <- Reduce(intersect, list(mr_data$snp, snp_coeff_id))
    if (length(snp_id) == 0L) {
      res_list[[i]] <- make_result_row(feature = val, ivw_n_total = 0L,
                                       error = "No overlapping SNPs")
      next
    }

    mr_data <- mr_data[match(snp_id, mr_data$snp), , drop = FALSE]
    snp_order <- match(snp_id, snp_coeff_id)
    snp_coeff_id_use <- snp_coeff_id[snp_order]
    snp_coeff_use <- snp_coeff[, snp_order, drop = FALSE]

    mr_data$n_xqtl <- qtl_sample_size
    std_zx <- gsmr2::std_effect(mr_data$A1F, mr_data$beta_qtl, mr_data$se_qtl, mr_data$n_xqtl)
    mr_data$beta_qtl <- std_zx$b
    mr_data$se_qtl <- std_zx$se

    ldrho <- tryCatch({
      R <- cor(snp_coeff_use)
      colnames(R) <- rownames(R) <- snp_coeff_id_use
      R
    }, error = function(e) NULL)

    if (is.null(ldrho)) {
      res_list[[i]] <- make_result_row(feature = val, error = "Failed to compute LD matrix")
      next
    }

    ivw_object <- tryCatch({
      mr_input_object <- mr_input(
        bx = mr_data$beta_qtl,
        bxse = mr_data$se_qtl,
        by = mr_data$beta,
        byse = mr_data$se,
        correlation = ldrho
      )
      mr_ivw(
        mr_input_object,
        model = "default",
        robust = FALSE,
        penalized = FALSE,
        correl = TRUE,
        weights = "simple",
        psi = 0,
        distribution = "normal",
        alpha = 0.05
      )
    }, error = function(e) e)

    if (inherits(ivw_object, "error")) {
      res_list[[i]] <- make_result_row(feature = val, error = "IVW error")
    } else {
      res_list[[i]] <- make_result_row(
        feature = val,
        ivw_beta = ivw_object@Estimate,
        ivw_se = ivw_object@StdError,
        ivw_p = ivw_object@Pvalue,
        ivw_n_total = length(mr_data$beta_qtl)
      )
    }
  }

  do.call(rbind, res_list)
}

files <- list.files(
  path = input_path,
  pattern = paste0("^", xqtl_data, ".*\\.tsv$"),
  full.names = TRUE
)

if (length(files) == 0L) {
  stop(sprintf("No LD input .tsv files found for dataset %s in %s", xqtl_data, input_path))
}

results_list <- mclapply(
  files,
  function(f) process_one(
    file = f,
    xqtl_data = xqtl_data,
    qtl_sample_size = qtl_sample_size,
    feature_col = feature_col
  ),
  mc.cores = jobs
)

results_df <- do.call(rbind, results_list)
dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(output_path, sprintf("IVW_%s.csv", xqtl_data))
write.csv(results_df, out_file, row.names = FALSE)
