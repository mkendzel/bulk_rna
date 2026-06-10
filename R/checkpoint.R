save_checkpoint <- function(obj, name, dir = "data/checkpoints", compress = TRUE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  
  # Find next version number based on existing files
  existing <- list.files(dir, pattern = paste0("^", name, "_v\\d+\\.rds$"))
  version <- length(existing) + 1
  
  path <- file.path(dir, sprintf("%s_v%02d.rds", name, version))
  saveRDS(obj, path, compress = compress)
  
  message("Saved: ", basename(path))
  invisible(path)
}

load_checkpoint <- function(name, version = "latest", dir = "data/checkpoints") {
  pattern <- paste0("^", name, "_v\\d+\\.rds$")
  files <- sort(list.files(dir, pattern = pattern, full.names = TRUE))
  
  if (length(files) == 0) stop("No checkpoints found for: ", name)
  
  path <- if (identical(version, "latest")) tail(files, 1) else files[version]
  
  if (is.na(path)) stop("Version ", version, " does not exist for: ", name)
  
  message("Loading: ", basename(path))
  readRDS(path)
}

list_checkpoints <- function(name = NULL, dir = "data/checkpoints") {
  pattern <- if (is.null(name)) "_v\\d+\\.rds$" else paste0("^", name, "_v\\d+\\.rds$")
  files <- sort(list.files(dir, pattern = pattern, full.names = TRUE))
  
  if (length(files) == 0) {
    message("No checkpoints found")
    return(invisible(NULL))
  }
  
  info <- file.info(files)
  data.frame(
    file    = basename(files),
    version = seq_along(files),
    size_mb = round(info$size / 1e6, 2),
    saved   = format(info$mtime, "%Y-%m-%d %H:%M"),
    row.names = NULL
  )
}