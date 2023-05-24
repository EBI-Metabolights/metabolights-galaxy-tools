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

createTestYamlObj <- function(args) {
  provided_config <- yaml::yaml.load_file(args[["params_file"]])
  provided_config$galaxy$gissmo_location <- args[["gissmo_ref"]]
  provided_config$files$spectral.matrix <- args[["matrix"]]
  return(provided_config)
}

create_yaml_obj_v2 <- function(args) {
  config <- list(
    dirs = list(
      temp = args[["dirs.temp"]],
      lib = args[["dirs.lib"]]
    ),
    study = list(
      id = args[["study.id"]],
      spectrometer.frequency = args[["study.spectrometer.frequency"]]
    ),
    files = list(
      spectral.matrix = args[["files.spectral.matrix"]]
    ),
    corrpockets = list(
      half.window = args[["corrpockets.half.window"]],
      noise.percentile = args[["corrpockets.noise.percentile"]],
      only.region.between = args[["corrpockets.only.region.between"]],
      rcutoff = args[["corrpockets.rcutoff"]]
    ),
    storm = list(
      correlation.r.cutoff = args[["storm.correlation.r.cutoff"]],
      q = args[["storm.q"]],
      b = args[["storm.b"]],
      number.of.plots = args[["storm.number.of.plots"]]
    ),
    tina = list(
      bounds = args[["tina.bounds"]],
      min.subset = args[["tina.min.subset"]],
      prom.ratio = args[["tina.prom.ratio"]],
      max.eps = args[["tina.max.eps"]],
      minPts = args[["tina.minPts"]],
      eps.stepsize = args[["tina.eps.stepsize"]],
      max.plots = args[["tina.max.plots"]],
      nfeats = args[["tina.nfeats"]]
    ),
    matching = list(
      ref.sig.SD.cutoff = args[["matching.ref.sig.SD.cutoff"]],
      max.hits = args[["matching.max.hits"]],
      r.thresh = args[["matching.r.thresh"]],
      p.thresh = args[["matching.p.thresh"]],
      filtering = list(
        res.area.threshold = args[["matching.filtering.res.area.threshold"]],
        ppm.tol = args[["matching.filtering.ppm.tol"]]
      )
    ),
    par = list(
      ncores = args[["par.ncores"]],
      type = args[["par.type"]]
    ),
    galaxy = list(
      enabled = args[["galaxy.enabled"]]
    ),
    debug = list(
      enabled = args[["debug.enabled"]],
      throttle_matches = args[["debug.throttle_matches"]]
    )
  )

  params_string <- yaml::as.yaml(config)
  params_obj <- yaml::yaml.load(params_string)

  return(params_obj)
}
