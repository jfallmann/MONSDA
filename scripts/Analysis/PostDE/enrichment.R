postde_read_term2gene <- function(path) {
    df <- postde_read_tsv(path, header = TRUE)
    if (!identical(sort(colnames(df)), c("gene", "term"))) {
        stop("term2gene must have exactly columns 'term' and 'gene'")
    }
    df <- df[, c("term", "gene")]
    df$term <- as.character(df$term)
    df$gene <- as.character(df$gene)
    if (any(is.na(df$term)) || any(!nzchar(df$term)) || any(is.na(df$gene)) || any(!nzchar(df$gene))) {
        stop("term2gene contains empty or NA ids")
    }
    df <- df[!duplicated(df), , drop = FALSE]
    df
}

postde_read_term2name <- function(path) {
    df <- postde_read_tsv(path, header = TRUE)
    if (!identical(sort(colnames(df)), c("name", "term"))) {
        stop("term2name must have exactly columns 'term' and 'name'")
    }
    df <- df[, c("term", "name")]
    df$term <- as.character(df$term)
    df$name <- as.character(df$name)
    if (any(is.na(df$term)) || any(!nzchar(df$term)) || any(is.na(df$name)) || any(!nzchar(df$name))) {
        stop("term2name contains empty or NA ids")
    }
    df <- df[!duplicated(df), , drop = FALSE]
    n_names <- tapply(df$name, df$term, function(x) length(unique(x)))
    bad <- names(n_names)[n_names > 1]
    if (length(bad) > 0) {
        stop(paste0("term2name has multiple distinct names for term ids: ", paste(bad, collapse = ", ")))
    }
    df
}

postde_build_gomap <- function(term2gene) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("clusterProfiler is required for go_expand but is not installed")
    }
    if (!requireNamespace("GO.db", quietly = TRUE)) {
        stop("GO.db is required for go_expand but is not installed")
    }
    if (!all(c("term", "gene") %in% colnames(term2gene))) {
        stop("term2gene must have columns 'term' and 'gene'")
    }
    term2gene <- term2gene[!duplicated(term2gene), , drop = FALSE]
    known <- ls(GO.db::GOTERM)
    unknown <- setdiff(unique(as.character(term2gene$term)), known)
    if (length(unknown) > 0) {
        stop(paste0("term2gene contains unknown GO terms: ", paste(head(unknown, 5), collapse = ", ")))
    }
    clusterProfiler::buildGOmap(term2gene)
}

postde_flatten_list_cols <- function(df) {
    for (col in colnames(df)) {
        if (is.list(df[[col]])) {
            df[[col]] <- vapply(df[[col]], function(x) paste(unlist(x), collapse = ";"), character(1))
        }
    }
    df
}

postde_run_gprofiler <- function(entry, config, outdir, gost_fun = NULL) {
    if (!requireNamespace("gprofiler2", quietly = TRUE)) {
        stop("gprofiler enabled but package gprofiler2 is not installed")
    }
    if (is.null(gost_fun)) {
        gost_fun <- gprofiler2::gost
    }
    gp <- config$gprofiler
    if (is.null(gp$organism) || !nzchar(gp$organism)) {
        stop("gprofiler: organism must be set explicitly (no default)")
    }
    organism <- gp$organism
    domain_scope <- if (is.null(gp$domain_scope)) "custom" else gp$domain_scope
    if (!domain_scope %in% c("annotated", "known", "custom", "custom_annotated")) {
        stop(paste0("gprofiler: invalid domain_scope '", domain_scope, "' (use annotated, known, custom or custom_annotated)"))
    }
    if (!isTRUE(gp$allow_network)) {
        stop("gprofiler: gost is an online call and requires gprofiler.allow_network=true")
    }
    sources <- unlist(gp$sources)
    if (length(sources) == 0) {
        sources <- NULL
    }
    correction_method <- if (is.null(gp$correction_method)) "g_SCS" else gp$correction_method
    if (!correction_method %in% c("g_SCS", "bonferroni", "fdr", "false_discovery_rate", "gSCS", "analytical")) {
        stop(paste0("gprofiler: invalid correction_method '", correction_method, "'"))
    }
    universe <- postde_ora_universe(entry$results)
    custom_bg <- universe
    if (!is.null(gp$background)) {
        bg <- postde_read_tsv(gp$background, header = FALSE)[[1]]
        custom_bg <- intersect(bg, universe)
        if (length(custom_bg) == 0) {
            stop("gprofiler: background does not intersect the tested universe")
        }
    }
    padj <- config$padj
    lfc <- config$lfc
    out <- list()
    for (direction in c("all", "up", "down")) {
        query <- postde_significant(entry$results, direction = direction, padj = padj, lfc = lfc)
        res_file <- file.path(outdir, paste0("gprofiler_", direction, "_result.tsv"))
        meta_file <- file.path(outdir, paste0("gprofiler_", direction, "_meta.rds"))
        query_file <- file.path(outdir, paste0("gprofiler_", direction, "_query.rds"))
        if (length(query) == 0) {
            postde_write_empty(res_file, "no significant genes")
            out[[direction]] <- list(status = "no significant genes", result = res_file)
            next
        }
        if (length(intersect(query, custom_bg)) == 0) {
            stop(paste0("gprofiler: query for direction '", direction, "' is disjoint from the background"))
        }
        saveRDS(list(query = query, universe = custom_bg), query_file)
        gost_args <- list(
            query = query,
            organism = organism,
            domain_scope = domain_scope,
            sources = sources,
            significant = FALSE,
            user_threshold = padj,
            correction_method = correction_method
        )
        if (domain_scope %in% c("custom", "custom_annotated")) {
            gost_args$custom_bg <- custom_bg
        }
        gost_res <- tryCatch(
            do.call(gost_fun, gost_args),
            error = function(e) {
                stop(paste0("gprofiler: gost call failed for direction '", direction, "': ", conditionMessage(e)))
            }
        )
        if (is.null(gost_res) || is.null(gost_res$result)) {
            postde_write_empty(res_file, "no enriched terms")
            out[[direction]] <- list(status = "no enriched terms", result = res_file)
            next
        }
        res_df <- postde_flatten_list_cols(as.data.frame(gost_res$result))
        write.table(res_df, res_file, sep = "\t", row.names = FALSE, quote = FALSE)
        saveRDS(gost_res$meta, meta_file)
        saveRDS(gost_res$result, file.path(outdir, paste0("gprofiler_", direction, "_raw.rds")))
        out[[direction]] <- list(status = "ok", result = res_file, meta = meta_file, query = query_file, raw = file.path(outdir, paste0("gprofiler_", direction, "_raw.rds")))
    }
    out
}

postde_run_enricher <- function(genes, universe, term2gene, term2name, min_size, max_size, outdir, direction) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("clusterprofiler enabled but package clusterProfiler is not installed")
    }
    res_file <- file.path(outdir, paste0("enricher_", direction, "_result.tsv"))
    header <- c("ID", "Description", "GeneRatio", "BgRatio", "pvalue", "p.adjust", "qvalue", "geneID", "Count")
    if (length(genes) == 0) {
        postde_write_empty(res_file, "no significant genes", header = header)
        return(list(status = "no significant genes", result = res_file))
    }
    if (anyDuplicated(genes)) {
        stop("significant gene list contains duplicate ids")
    }
    if (length(intersect(genes, term2gene$gene)) == 0) {
        stop("no overlap between significant genes and term2gene genes")
    }
    enr_args <- list(gene = genes, universe = universe, TERM2GENE = term2gene, pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = min_size, maxGSSize = max_size)
    if (is.null(term2name)) {
        enr_args$TERM2NAME <- NA
    } else {
        enr_args$TERM2NAME <- term2name
    }
    enr <- do.call(clusterProfiler::enricher, enr_args)
    res <- as.data.frame(enr)
    if (nrow(res) == 0) {
        postde_write_empty(res_file, "no enriched terms", header = header)
        return(list(status = "no enriched terms", result = res_file))
    }
    write.table(res, res_file, sep = "\t", row.names = FALSE, quote = FALSE)
    list(status = "ok", result = res_file)
}

postde_run_gsea <- function(entry, term2gene, term2name, min_size, max_size, seed, outdir) {
    if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
        stop("clusterprofiler enabled but package clusterProfiler is not installed")
    }
    res_file <- file.path(outdir, "gsea_result.tsv")
    header <- c("ID", "Description", "setSize", "enrichmentScore", "NES", "pvalue", "p.adjust", "qvalues", "rank", "leading_edge", "core_enrichment")
    if (!"stat" %in% colnames(entry$results)) {
        postde_write_empty(res_file, "rank unavailable: results have no stat column", header = header)
        return(list(status = "rank unavailable", result = res_file))
    }
    ranked <- postde_rank_genes(entry$results, stat_col = "stat")
    gene_list <- setNames(ranked$stat, ranked$gene_id)
    saveRDS(gene_list, file.path(outdir, "gsea_ranked.rds"))
    if (length(intersect(names(gene_list), term2gene$gene)) == 0) {
        stop("no overlap between ranked genes and term2gene genes")
    }
    gsea_args <- list(geneList = gene_list, TERM2GENE = term2gene, seed = seed, pvalueCutoff = 1, minGSSize = min_size, maxGSSize = max_size)
    if (is.null(term2name)) {
        gsea_args$TERM2NAME <- NA
    } else {
        gsea_args$TERM2NAME <- term2name
    }
    if ("by" %in% names(formals(clusterProfiler::GSEA))) {
        gsea_args$by <- "fgsea"
    } else {
        gsea_args$method <- "fgsea"
    }
    gsea_res <- do.call(clusterProfiler::GSEA, gsea_args)
    res <- as.data.frame(gsea_res)
    if (nrow(res) == 0) {
        postde_write_empty(res_file, "no enriched terms", header = header)
        return(list(status = "no enriched terms", result = res_file))
    }
    write.table(res, res_file, sep = "\t", row.names = FALSE, quote = FALSE)
    list(status = "ok", result = res_file, ranked = file.path(outdir, "gsea_ranked.rds"))
}

postde_run_gsva <- function(entry, term2gene, term2name, config, outdir) {
    if (!requireNamespace("GSVA", quietly = TRUE)) {
        stop("gsva enabled but package GSVA is not installed")
    }
    gs <- config$gsva
    method <- if (is.null(gs$method)) "gsva" else gs$method
    min_size <- if (is.null(gs$min_size)) config$min_size else gs$min_size
    max_size <- if (is.null(gs$max_size)) config$max_size else gs$max_size
    gene_sets <- split(term2gene$gene, term2gene$term)
    saveRDS(gene_sets, file.path(outdir, "gsva_sets.rds"))
    saveRDS(list(term2gene = term2gene, term2name = term2name, cohort = entry$id), file.path(outdir, "gsva_provenance.rds"))
    scores <- postde_gsva_scores(entry$expression, gene_sets, method = method, min_size = min_size, max_size = max_size)
    write.table(scores, file.path(outdir, "gsva_scores.tsv"), sep = "\t", col.names = NA, quote = FALSE)
    diff <- postde_differential_scores(scores, entry$design, entry$contrast)
    write.table(diff$results, file.path(outdir, "gsva_differential.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(diff$design, file.path(outdir, "gsva_design.tsv"), sep = "\t", col.names = NA, quote = FALSE)
    write.table(data.frame(term = names(diff$contrast), coefficient = diff$contrast), file.path(outdir, "gsva_contrast.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
    list(status = "ok", scores = file.path(outdir, "gsva_scores.tsv"), differential = file.path(outdir, "gsva_differential.tsv"))
}

postde_run_enrichment <- function(entry, config, outdir) {
    cp <- config$clusterprofiler
    gs <- config$gsva
    cp_enabled <- is.list(cp) && isTRUE(cp$enabled)
    gs_enabled <- is.list(gs) && isTRUE(gs$enabled)
    if (!cp_enabled && !gs_enabled) {
        return(list())
    }
    term2gene <- NULL
    term2name <- NULL
    if (cp_enabled || gs_enabled) {
        if (is.null(cp$term2gene)) {
            stop("term2gene path missing (clusterprofiler.term2gene required for clusterprofiler/gsva)")
        }
        term2gene <- postde_read_term2gene(cp$term2gene)
        if (!is.null(cp$term2name)) {
            term2name <- postde_read_term2name(cp$term2name)
        }
        if (isTRUE(cp$go_expand)) {
            term2gene <- postde_build_gomap(term2gene)
        }
    }
    universe <- postde_ora_universe(entry$results)
    if (!is.null(cp$background)) {
        bg <- postde_read_tsv(cp$background, header = FALSE)[[1]]
        universe <- intersect(bg, universe)
        if (length(universe) == 0) {
            stop("clusterprofiler: background does not intersect the tested universe")
        }
    }
    out <- list()
    if (cp_enabled) {
        ranked <- postde_rank_genes(entry$results, stat_col = "stat")
        provenance <- list(
            term2gene = term2gene,
            term2name = term2name,
            universe = universe,
            ranked = ranked,
            checksum = if (file.exists(cp$term2gene)) unname(tools::md5sum(cp$term2gene)) else NA_character_,
            config = cp
        )
        saveRDS(provenance, file.path(outdir, "enrichment_provenance.rds"))
        for (direction in c("all", "up", "down")) {
            genes <- postde_significant(entry$results, direction = direction, padj = config$padj, lfc = config$lfc)
            out[[paste0("enricher_", direction)]] <- postde_run_enricher(genes, universe, term2gene, term2name, config$min_size, config$max_size, outdir, direction)
        }
        out$gsea <- postde_run_gsea(entry, term2gene, term2name, config$min_size, config$max_size, config$seed, outdir)
    }
    if (gs_enabled) {
        out$gsva <- postde_run_gsva(entry, term2gene, term2name, config, outdir)
    }
    out
}
