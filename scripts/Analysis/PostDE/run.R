postde_parse_cli <- function(args) {
    out <- list(bundle = NULL, config = NULL, output = NULL)
    i <- 1
    while (i <= length(args)) {
        a <- args[i]
        if (grepl("^--bundle=", a)) {
            out$bundle <- sub("^--bundle=", "", a)
        } else if (a == "--bundle") {
            if (i >= length(args)) stop("--bundle requires a value")
            out$bundle <- args[i + 1]
            i <- i + 1
        } else if (grepl("^--config=", a)) {
            out$config <- sub("^--config=", "", a)
        } else if (a == "--config") {
            if (i >= length(args)) stop("--config requires a value")
            out$config <- args[i + 1]
            i <- i + 1
        } else if (grepl("^--output=", a)) {
            out$output <- sub("^--output=", "", a)
        } else if (a == "--output") {
            if (i >= length(args)) stop("--output requires a value")
            out$output <- args[i + 1]
            i <- i + 1
        } else {
            stop(paste0("unknown argument: ", a))
        }
        i <- i + 1
    }
    if (is.null(out$bundle) || is.null(out$config) || is.null(out$output)) {
        stop("usage: Rscript run.R --bundle FILE --config JSON --output DIR")
    }
    out
}

postde_config_defaults <- function(config) {
    if (is.null(config$padj)) config$padj <- 0.05
    if (is.null(config$lfc)) config$lfc <- 1
    if (is.null(config$seed)) config$seed <- 1
    if (is.null(config$min_size)) config$min_size <- 10
    if (is.null(config$max_size)) config$max_size <- 500
    config
}

postde_config_enabled <- function(sub) {
    is.list(sub) && isTRUE(sub$enabled)
}

postde_validate_config <- function(config) {
    if (!is.null(config$padj) && (!is.numeric(config$padj) || config$padj <= 0 || config$padj >= 1)) {
        stop("config padj must be numeric in (0, 1)")
    }
    if (!is.null(config$lfc) && (!is.numeric(config$lfc) || config$lfc < 0)) {
        stop("config lfc must be numeric >= 0")
    }
    if (!is.null(config$seed) && !is.numeric(config$seed)) {
        stop("config seed must be numeric")
    }
    if (!is.null(config$min_size) && (!is.numeric(config$min_size) || config$min_size < 1)) {
        stop("config min_size must be numeric >= 1")
    }
    if (!is.null(config$max_size) && (!is.numeric(config$max_size) || config$max_size < config$min_size)) {
        stop("config max_size must be numeric >= min_size")
    }
    for (sub in c("gprofiler", "clusterprofiler", "gsva", "decoupler", "dream", "plots", "report", "shiny")) {
        if (!is.null(config[[sub]]) && !is.list(config[[sub]])) {
            stop(paste0("config ", sub, " must be an object"))
        }
        if (is.list(config[[sub]]) && !is.null(config[[sub]]$enabled) && !is.logical(config[[sub]]$enabled)) {
            stop(paste0("config ", sub, ".enabled must be a boolean"))
        }
    }
    invisible(TRUE)
}

postde_validate_bundle <- function(bundle) {
    if (!is.list(bundle) || is.null(bundle$schema_version) || is.null(bundle$engine) || is.null(bundle$contrasts)) {
        stop("bundle is not a valid postde bundle (schema_version, engine, contrasts required)")
    }
    if (!identical(bundle$schema_version, 1L) && !identical(bundle$schema_version, 1)) {
        stop(paste0("unsupported bundle schema_version: ", bundle$schema_version))
    }
    if (length(bundle$contrasts) == 0) {
        stop("bundle contains no contrast entries")
    }
    ids <- vapply(bundle$contrasts, function(e) e$id, character(1))
    if (anyDuplicated(ids)) {
        stop(paste0("bundle contains duplicate entry ids: ", paste(ids[duplicated(ids)], collapse = ", ")))
    }
    postde_check_sanitized_ids(ids)
    for (id in ids) {
        if (grepl("[/\\\\]", id) || id == "." || id == "..") {
            stop(paste0("bundle entry id is unsafe: '", id, "'"))
        }
    }
    invisible(TRUE)
}

postde_script_dir <- function() {
    args <- commandArgs(trailingOnly = FALSE)
    file_arg <- sub("^--file=", "", args[grep("^--file=", args)])
    if (length(file_arg) == 0) {
        stop("Rscript --file path not found in commandArgs")
    }
    dirname(normalizePath(file_arg))
}

postde_relpath <- function(path, base) {
    path <- normalizePath(path, mustWork = FALSE)
    base <- normalizePath(base, mustWork = FALSE)
    if (identical(path, base)) {
        return(".")
    }
    if (startsWith(path, paste0(base, .Platform$file.sep))) {
        return(substr(path, nchar(base) + 2, nchar(path)))
    }
    path
}

postde_relativize <- function(x, base) {
    if (is.character(x) && length(x) == 1 && !is.na(x) && startsWith(x, base)) {
        return(postde_relpath(x, base))
    }
    if (is.list(x)) {
        return(lapply(x, postde_relativize, base = base))
    }
    x
}

postde_main <- function() {
    cli <- postde_parse_cli(commandArgs(trailingOnly = TRUE))
    if (!file.exists(cli$bundle)) {
        stop(paste0("bundle file not found: ", cli$bundle))
    }
    if (!file.exists(cli$config)) {
        stop(paste0("config file not found: ", cli$config))
    }
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
        stop("jsonlite is required but is not installed")
    }
    script_dir <- postde_script_dir()
    source(file.path(script_dir, "common.R"))
    source(file.path(script_dir, "enrichment.R"))
    source(file.path(script_dir, "regulatory.R"))
    bundle <- readRDS(cli$bundle)
    config <- postde_config_defaults(jsonlite::fromJSON(cli$config, simplifyVector = FALSE))
    postde_validate_bundle(bundle)
    postde_validate_config(config)
    if (!dir.exists(cli$output)) {
        dir.create(cli$output, recursive = TRUE)
    }
    set.seed(config$seed)
    outputs <- list()
    for (entry_name in names(bundle$contrasts)) {
        entry <- bundle$contrasts[[entry_name]]
        if (!is.list(entry) || is.null(entry$id)) {
            stop(paste0("bundle contrast '", entry_name, "' is not a valid entry (id required)"))
        }
        entry_dir <- file.path(cli$output, postde_sanitize_id(entry$id))
        dir.create(entry_dir, recursive = TRUE, showWarnings = FALSE)
        entry_outputs <- list()
        if (postde_config_enabled(config$gprofiler)) {
            entry_outputs$gprofiler <- postde_run_gprofiler(entry, config, entry_dir)
        }
        if (postde_config_enabled(config$clusterprofiler) || postde_config_enabled(config$gsva)) {
            entry_outputs$enrichment <- postde_run_enrichment(entry, config, entry_dir)
        }
        if (postde_config_enabled(config$decoupler) || postde_config_enabled(config$dream)) {
            entry_outputs$regulatory <- postde_run_regulatory(entry, config, entry_dir)
        }
        outputs[[entry$id]] <- list(id = entry$id, dir = entry_dir, analyses = entry_outputs)
    }
    report_data <- list(config = config, bundle = bundle, outputs = outputs)
    saveRDS(report_data, file.path(cli$output, "report_data.rds"))
    manifest <- list(
        bundle = postde_relpath(cli$bundle, cli$output),
        config = postde_relpath(cli$config, cli$output),
        output = ".",
        bundle_checksum = unname(tools::md5sum(cli$bundle)),
        config_checksum = unname(tools::md5sum(cli$config)),
        entries = lapply(outputs, function(o) list(id = o$id, dir = postde_relpath(o$dir, cli$output), analyses = postde_relativize(o$analyses, cli$output))),
        sessionInfo = capture.output(utils::sessionInfo())
    )
    jsonlite::write_json(manifest, file.path(cli$output, "manifest.json"), pretty = TRUE, auto_unbox = TRUE)
    any_vis <- postde_config_enabled(config$plots) || postde_config_enabled(config$report) || postde_config_enabled(config$shiny)
    if (any_vis) {
        reporting <- file.path(script_dir, "reporting.R")
        if (file.exists(reporting)) {
            source(reporting)
        }
        if (exists("render_postde", mode = "function")) {
            render_postde(bundle, config, cli$output)
        } else {
            stop("plots/report/shiny enabled but render_postde() is not available (reporting.R not yet provided)")
        }
    }
    invisible(report_data)
}

postde_main()
