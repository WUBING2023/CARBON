# Run from an ASCII-only working directory containing this source bundle.
lib <- normalizePath(".", winslash="/")
lib <- file.path(lib,"Rlibrary")
dir.create(lib,showWarnings=FALSE)
.libPaths(c(lib,.libPaths()))
for (name in c("GlobalOptions","shape","circlize","ggseqlogo")) {
  if (!requireNamespace(name,quietly=TRUE)) {
    tarball <- list.files("package_sources",pattern=paste0("^",name,"_[0-9].*\\.tar\\.gz$"),full.names=TRUE)
    stopifnot(length(tarball)==1L)
    install.packages(tarball,repos=NULL,type="source",lib=lib)
  }
}
stopifnot(requireNamespace("circlize",quietly=TRUE),requireNamespace("ggseqlogo",quietly=TRUE))
writeLines(capture.output(sessionInfo()),"R_SESSION_REPLAY.txt")
