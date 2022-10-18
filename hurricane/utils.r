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
        if (args[key] %in% c("TRUE", "FALSE")) {
            args[key] <- as.logical(args[key])
        }
    }
    return(args)
}

createYamlObj <- function(args) {
    params_string <- yaml::as.yaml(
        list(
            general_pars = list(
                peaks_location = args[["peaks_loc"]],
                spec_location = args[["spec_loc"]],
                output_dir = args[["output_dir"]]
            ),
            sd_pars = list(
                cutoff = args[["cutoff"]]
            ),
            am_pars = list(
                rank_limit = args[["rank_limit"]],
                dist_thresh = args[["dist_thresh"]],
                matchMethod = args[["match_method"]],
                refdb_file = args[["refdb_file"]]
            ),
            vis_pars = list(
                matlab_root = ""
            )
        )
    )
    params_obj <- yaml::yaml.load(params_string)
    return(params_obj)
}
