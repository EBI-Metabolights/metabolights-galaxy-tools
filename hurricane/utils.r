loadAndDisplayPackages <- function(pkgs) {
    for (pkg in pkgs) suppressPackageStartupMessages(stopifnot(library(pkg, quietly = TRUE, logical.return = TRUE, character.only = TRUE)))

    sessioninfo <- sessionInfo()
    cat(sessioninfo$R.version$version.string, "\n")
    cat("Main packages:\n")
    for (pkg in names(sessioninfo$otherPkgs)) {
      cat(paste(pkg, packageVersion(pkg)), "\t")
    }
    cat("\n")
    cat("Other loaded packages:\n")
    for (pkg in names(sessioninfo$loadedOnly)) {
      cat(paste(pkg, packageVersion(pkg)), "\t")
    }
    cat("\n")
}

parseCommandArgs <- function(...) {
    args <- batch::parseCommandArgs(...)
    for (key in names(args)) {
        if (args[key] %in% c("TRUE", "FALSE"))
            args[key] <- as.logical(args[key])
    }
    return(args)
}
