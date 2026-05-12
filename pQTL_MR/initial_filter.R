library(ieugwasr)
library(readr)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
output_file <- args[2]
feature_col <- args[3]
r2_thresh <- as.numeric(args[4])
pval_thresh <- as.numeric(args[5])
plink_path <- args[6]
bfile_path <- args[7]

X <- read_tsv(input_file, show_col_types = FALSE)

feature_vals <- unique(X[[feature_col]])
res_list <- vector("list", length(feature_vals))

for (i in seq_along(feature_vals)) {
    val <- feature_vals[i]
    df <- X %>% filter(.data[[feature_col]] == val)

    X_clump <- tryCatch(
        {
            df %>%
                rename(rsid = snp, pval = p_qtl) %>%
                ieugwasr::ld_clump(
                    dat = .,
                    clump_r2 = r2_thresh,
                    clump_p = pval_thresh,
                    plink_bin = plink_path,
                    bfile = bfile_path
                )
        },
        error = function(e) {
            NULL
        }
    )

    if (is.null(X_clump) || nrow(X_clump) == 0) {
        res_list[[i]] <- df[0, ]
    } else {
        top_vars <- unique(X_clump$rsid)
        res_list[[i]] <- df %>% filter(snp %in% top_vars)
    }
}

X_subset <- bind_rows(res_list)
write_tsv(X_subset, output_file)
