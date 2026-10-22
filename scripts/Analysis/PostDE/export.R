postde_common_path <- function() {
    frames <- sys.frames()
    for (fr in frames) {
        if (exists("postde_lib", envir = fr, inherits = FALSE)) {
            return(file.path(dirname(get("postde_lib", envir = fr)), "common.R"))
        }
    }
    NULL
}

postde_common <- postde_common_path()
if (!is.null(postde_common) && file.exists(postde_common)) {
    source(postde_common)
} else {
    postde_contrast <- function(formula, metadata, A, B) {
        md <- as.data.frame(metadata)
        if (!"condition" %in% colnames(md)) {
            stop("postde_contrast: metadata has no condition column")
        }
        md$condition <- relevel(md$condition, ref = B)
        mdA <- md
        mdB <- md
        mdA$condition <- factor(rep(A, nrow(md)), levels = levels(md$condition))
        mdB$condition <- factor(rep(B, nrow(md)), levels = levels(md$condition))
        mA <- model.matrix(formula, data = mdA)
        mB <- model.matrix(formula, data = mdB)
        setNames(colMeans(mA) - colMeans(mB), colnames(mA))
    }
}

postde_bundles <- list()
postde_entry_ids <- character(0)
if (!exists("postde_cache_dir")) {
    postde_cache_dir <- NULL
}

postde_results_deseq2 <- function(res) {
    out <- data.frame(
        gene_id = rownames(res),
        logFC = as.numeric(res$log2FoldChange),
        pvalue = as.numeric(res$pvalue),
        padj = as.numeric(res$padj),
        stat = as.numeric(res$stat),
        mean = as.numeric(res$baseMean),
        stringsAsFactors = FALSE
    )
    rownames(out) <- rownames(res)
    out
}

postde_results_edger <- function(qlf) {
    F <- qlf$table$F
    if (any(is.finite(F) & F < 0)) {
        stop("postde_results_edger: negative F statistics found")
    }
    df_test <- if (!is.null(qlf$df.test)) qlf$df.test else 1
    one_df <- isTRUE(all(df_test == 1))
    stat <- ifelse(is.finite(F) & F >= 0 & one_df, sign(qlf$table$logFC) * sqrt(F), NA_real_)
    out <- data.frame(
        gene_id = rownames(qlf$table),
        logFC = as.numeric(qlf$table$logFC),
        pvalue = as.numeric(qlf$table$PValue),
        padj = p.adjust(qlf$table$PValue, method = "BH"),
        stat = stat,
        mean = as.numeric(qlf$table$logCPM),
        stringsAsFactors = FALSE
    )
    rownames(out) <- rownames(qlf$table)
    out
}

postde_capture <- function(engine, id, A, B, normalized, metadata, counts, expression, formula, results, mean_scale, library_normalized = TRUE) {
    metadata <- as.data.frame(metadata)
    counts <- as.matrix(counts)
    expression <- as.matrix(expression)
    results <- as.data.frame(results)
    if (is.null(rownames(results)) || is.null(rownames(expression))) {
        stop(paste0("postde: results and expression must have rownames for entry '", id, "'"))
    }
    if (!identical(rownames(results), rownames(expression))) {
        stop(paste0("postde: results rownames do not match expression rownames for entry '", id, "'"))
    }
    if (anyDuplicated(rownames(results))) {
        stop(paste0("postde: results contain duplicate gene ids for entry '", id, "'"))
    }
    if (!identical(rownames(metadata), colnames(counts))) {
        stop(paste0("postde: metadata rownames do not match counts colnames for entry '", id, "'"))
    }
    if (!identical(colnames(counts), colnames(expression))) {
        stop(paste0("postde: counts colnames do not match expression colnames for entry '", id, "'"))
    }
    if (!all(rownames(results) %in% rownames(counts))) {
        stop(paste0("postde: results contain genes absent from counts for entry '", id, "'"))
    }
    if (!all(rownames(results) == results$gene_id)) {
        stop(paste0("postde: results rownames do not match gene_id column for entry '", id, "'"))
    }
    md <- metadata
    md$condition <- relevel(md$condition, ref = B)
    design <- model.matrix(formula, data = md)
    contrast <- postde_contrast(formula, metadata, A, B)
    postde_validate_design(design, contrast)
    entry <- list(
        id = id,
        A = A,
        B = B,
        normalized = isTRUE(normalized),
        library_normalized = isTRUE(library_normalized),
        metadata = metadata,
        counts = counts,
        expression = expression,
        design = design,
        contrast = contrast,
        results = results,
        mean_scale = mean_scale
    )
    # Spill each entry to its own file immediately rather than keeping every
    # contrast's counts/expression matrices resident for the whole run: with
    # dozens of contrasts (and a _norm variant for each) the accumulated list
    # has been large enough to exhaust job memory before the final save. Only
    # the small id and file path are kept in memory; postde_write() streams
    # the files back together one at a time for the single combined saveRDS.
    if (!is.null(postde_cache_dir)) {
        san_id <- postde_sanitize_id(id)
        entry_path <- file.path(postde_cache_dir, paste0(san_id, ".rds"))
        saveRDS(entry, entry_path)
        postde_entry_ids <<- c(postde_entry_ids, id)
        postde_bundles[[length(postde_bundles) + 1L]] <<- entry_path
        rm(entry)
    } else {
        postde_entry_ids <<- c(postde_entry_ids, id)
        postde_bundles[[length(postde_bundles) + 1L]] <<- entry
    }
    invisible(NULL)
}

postde_write <- function(outdir, combi, engine) {
    if (length(postde_bundles) == 0L) {
        stop("MONSDA_POSTDE=1 but no contrast entries were captured; refusing to write an empty bundle")
    }
    ids <- postde_entry_ids
    if (anyDuplicated(ids)) {
        stop(paste0("postde: duplicate entry ids: ", paste(ids[duplicated(ids)], collapse = ", ")))
    }
    postde_check_sanitized_ids(ids)
    contrasts <- setNames(vector("list", length(postde_bundles)), ids)
    for (i in seq_along(postde_bundles)) {
        item <- postde_bundles[[i]]
        contrasts[[i]] <- if (is.character(item)) readRDS(item) else item
    }
    bundle <- list(
        schema_version = 1L,
        engine = engine,
        contrasts = contrasts
    )
    path <- file.path(outdir, paste0("DE_", engine, "_", combi, "_postde.rds"))
    saveRDS(bundle, path)
    if (!is.null(postde_cache_dir)) {
        unlink(postde_cache_dir, recursive = TRUE)
    }
    invisible(path)
}
