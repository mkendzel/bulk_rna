# Helper to write a figure once per format in FIG_FORMATS
FIG_FORMATS <- c("pdf", "tiff")

save_fig <- function(p, name, width = 7, height = 6, subdir = NULL) {
  outdir <- if (is.null(subdir)) ensure_dir(dir_fig) else ensure_dir(dir_fig, subdir)

  files <- vapply(FIG_FORMATS, function(fmt) {
    f <- file.path(outdir, paste0(name, ".", fmt))

    if (fmt == "tiff") {
      tiff(f, width = width, height = height, units = "in", res = 300, compression = "lzw")
    } else {
      pdf(f, width = width, height = height)
    }
    on.exit(dev.off(), add = TRUE)

    if (inherits(p, "pheatmap")) grid.draw(p$gtable) else print(p)
    f
  }, character(1))

  invisible(unname(files))
}
