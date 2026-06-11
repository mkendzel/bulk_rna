# Auto-documents a checkpoint alongside the saved .rds.
# Each gatherer is wrapped so a failure to collect metadata cannot block
# or break saving the object.

# Packages listed in the header; everything else goes in the appendix.
.checkpoint_key_pkgs <- c(
  "DESeq2", "edgeR", "apeglm", "fgsea", "AnnotationDbi", "limma",
  "org.Hs.eg.db", "org.Mm.eg.db", "clusterProfiler", "dplyr", "purrr", "tidyr"
)

# Parses a line spec into integer line indices. Accepts dash ranges
# ("1-14,18-20"), R-style colon ranges with an optional c() wrapper
# ("c(1:10,12:40)"), or a numeric vector. Returns integer(0) for empty/NULL.
parse_line_ranges <- function(spec) {
  if (is.null(spec) || length(spec) == 0) return(integer(0))
  if (is.numeric(spec)) return(as.integer(spec))

  spec <- gsub("\\s+", "", as.character(spec))
  spec <- sub("^c\\((.*)\\)$", "\\1", spec)
  if (!nzchar(spec)) return(integer(0))

  parts <- strsplit(spec, ",", fixed = TRUE)[[1]]
  idx <- lapply(parts, function(p) {
    if (grepl("[-:]", p)) {
      ends <- as.integer(strsplit(p, "[-:]")[[1]])
      if (length(ends) != 2L || anyNA(ends)) stop("Bad line range: ", p)
      if (ends[2] < ends[1]) stop("Descending line range: ", p)
      seq.int(ends[1], ends[2])
    } else {
      n <- as.integer(p)
      if (is.na(n)) stop("Bad line number: ", p)
      n
    }
  })
  unique(unlist(idx))
}

# Collapses sorted integer indices into a compact "1-10, 12-40" string.
.format_line_spec <- function(idx) {
  idx <- sort(unique(as.integer(idx)))
  if (length(idx) == 0) return("")
  breaks <- c(0, which(diff(idx) != 1), length(idx))
  runs <- vapply(seq_len(length(breaks) - 1), function(i) {
    run <- idx[(breaks[i] + 1):breaks[i + 1]]
    if (length(run) == 1) as.character(run) else sprintf("%d-%d", run[1], run[length(run)])
  }, character(1))
  paste(runs, collapse = ", ")
}

# Reads the lines named by spec from the script focused in the RStudio editor.
# Returns a list(script, spec, code) or NULL if unavailable.
capture_script_lines <- function(spec) {
  tryCatch({
    if (is.null(spec)) return(NULL)
    if (!requireNamespace("rstudioapi", quietly = TRUE) ||
        !rstudioapi::isAvailable()) return(NULL)

    ctx <- rstudioapi::getSourceEditorContext()
    if (is.null(ctx)) return(NULL)

    contents <- ctx$contents
    idx <- parse_line_ranges(spec)
    idx <- sort(unique(idx[idx >= 1 & idx <= length(contents)]))
    if (length(idx) == 0) return(NULL)

    script <- if (nzchar(ctx$path)) basename(ctx$path) else "(unsaved editor)"
    list(script = script, spec = .format_line_spec(idx), code = contents[idx])
  }, error = function(e) NULL)
}

# Collects environment, git, and session facts. Never errors.
capture_provenance <- function() {
  prov <- list(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    user      = tryCatch(unname(Sys.info()[["user"]]),     error = function(e) NA),
    host      = tryCatch(unname(Sys.info()[["nodename"]]), error = function(e) NA),
    r_version = R.version.string,
    platform  = R.version$platform
  )

  # Active script path, when RStudio is available.
  prov$script <- tryCatch({
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      p <- rstudioapi::getSourceEditorContext()$path
      if (nzchar(p)) p else NA
    } else NA
  }, error = function(e) NA)

  # Loaded packages and their versions.
  prov$packages <- tryCatch({
    loaded <- loadedNamespaces()
    vers <- vapply(loaded, function(p) {
      tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
    }, character(1))
    vers <- vers[!is.na(vers)]
    vers[order(names(vers))]
  }, error = function(e) character(0))

  # Git commit and dirty state.
  prov$git <- tryCatch({
    sha <- system2("git", c("rev-parse", "--short", "HEAD"),
                   stdout = TRUE, stderr = FALSE)
    status <- system2("git", c("status", "--porcelain"),
                      stdout = TRUE, stderr = FALSE)
    if (length(sha) == 0 || !nzchar(sha[1])) {
      NA
    } else {
      n_dirty <- length(status)
      dirty <- if (n_dirty > 0) sprintf("dirty: %d files modified", n_dirty) else "clean"
      sprintf("%s (%s)", sha[1], dirty)
    }
  }, error = function(e) NA)

  prov
}

# Dispatches on object type and returns a character vector of bullet bodies.
describe_object <- function(obj) {
  tryCatch({
    if (methods::is(obj, "DESeqDataSet")) {
      describe_deseq(obj)
    } else if (is.matrix(obj)) {
      describe_matrix(obj)
    } else if (is.data.frame(obj)) {
      describe_data_frame(obj)
    } else if (is.list(obj)) {
      describe_list(obj)
    } else {
      describe_default(obj)
    }
  }, error = function(e) describe_default(obj))
}

# Builds the one-line type and size summary used in the header.
object_headline <- function(obj) {
  tryCatch({
    cls <- class(obj)[1]
    size <- format(utils::object.size(obj), units = "auto")
    dim_str <- if (!is.null(dim(obj))) {
      paste(dim(obj), collapse = " x ")
    } else {
      paste0("length ", length(obj))
    }
    sprintf("%s, %s, %s", cls, dim_str, size)
  }, error = function(e) class(obj)[1])
}

# Joins names, truncating to n with a count suffix when longer.
.trunc_names <- function(x, n = 12) {
  if (length(x) <= n) return(paste(x, collapse = ", "))
  paste0(paste(utils::head(x, n), collapse = ", "), ", ... (", length(x), " total)")
}

describe_matrix <- function(obj) {
  c(
    sprintf("Dimensions: %d rows x %d cols", nrow(obj), ncol(obj)),
    if (!is.null(colnames(obj)))
      sprintf("Columns/samples (%d): %s", ncol(obj), .trunc_names(colnames(obj)))
  )
}

describe_data_frame <- function(obj) {
  c(
    sprintf("Dimensions: %d rows x %d cols", nrow(obj), ncol(obj)),
    sprintf("Columns (%d): %s", ncol(obj), .trunc_names(names(obj)))
  )
}

describe_list <- function(obj) {
  nm <- names(obj)
  rows <- vapply(obj, function(x) {
    if (!is.null(nrow(x))) as.character(nrow(x)) else as.character(length(x))
  }, character(1))
  c(
    sprintf("List of %d element(s)", length(obj)),
    if (!is.null(nm))
      sprintf("Elements: %s", .trunc_names(paste0(nm, " (n=", rows, ")")))
  )
}

describe_deseq <- function(obj) {
  design_str <- tryCatch(paste(deparse(DESeq2::design(obj)), collapse = " "),
                         error = function(e) NA)
  cond <- tryCatch(SummarizedExperiment::colData(obj)$condition,
                  error = function(e) NULL)
  out <- c(
    sprintf("DESeqDataSet: %d genes x %d samples", nrow(obj), ncol(obj)),
    if (!is.na(design_str)) sprintf("Design: %s", design_str)
  )
  if (is.factor(cond)) {
    out <- c(out,
      sprintf("Reference condition: %s", levels(cond)[1]),
      sprintf("Conditions: %s", paste(levels(cond), collapse = ", "))
    )
  }
  sf <- tryCatch(!is.null(DESeq2::sizeFactors(obj)), error = function(e) FALSE)
  disp_set <- tryCatch(!all(is.na(DESeq2::dispersions(obj))), error = function(e) FALSE)
  out <- c(out,
    sprintf("Size factors estimated: %s", sf),
    sprintf("Dispersions estimated: %s", disp_set)
  )
  out
}

describe_default <- function(obj) {
  c(
    sprintf("Class: %s", paste(class(obj), collapse = ", ")),
    sprintf("Length: %d", length(obj))
  )
}

# Assembles the Markdown document and writes it to path.
# Assembles the Markdown document and writes it to path.
render_checkpoint_md <- function(path, name, version, obj,
                                 prov, code, notes = NULL) {
  L <- character(0)
  add <- function(...) L[[length(L) + 1]] <<- paste0(...)

  add("# Checkpoint: ", name, " (v", sprintf("%02d", version), ")")
  add("")
  add("- Saved:    ", prov$timestamp)
  add("- Saved by: ", prov$user, " @ ", prov$host)
  add("- Object:   ", object_headline(obj))
  add("")

  add("## Key decisions (auto-detected)")
  for (line in describe_object(obj)) add("- ", line)
  if (!is.null(notes)) add("- Notes: ", notes)
  add("")

  add("## Provenance")
  if (!is.na(prov$script)) add("- Active script: ", prov$script)
  if (!is.na(prov$git))    add("- Git commit:    ", prov$git)
  add("- R version:     ", prov$r_version)
  add("- Platform:      ", prov$platform)
  key <- prov$packages[names(prov$packages) %in% .checkpoint_key_pkgs]
  if (length(key) > 0) {
    add("- Key packages:  ",
        paste(sprintf("%s %s", names(key), key), collapse = ", "))
  }
  add("")

  add("## Code")
  if (is.null(code)) {
    add("_No code lines were specified (pass `lines = \"1-30\"` or `lines = c(1:30)`",
        " to record them)_")
  } else {
    add("Lines ", code$spec, " of ", code$script, ":")
    add("")
    add("```r")
    for (line in code$code) add(line)
    add("```")
  }
  add("")

  if (length(prov$packages) > 0) {
    add("## Session packages")
    add("")
    add("```")
    for (p in names(prov$packages)) add(p, " ", prov$packages[[p]])
    add("```")
    add("")
  }

  writeLines(unlist(L), path)
  invisible(path)
}