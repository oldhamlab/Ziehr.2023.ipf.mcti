# cosmos.R

set_carnival_options <- function(){
  path <- "analysis/carnival"
  carnival_options <- cosmosR::default_CARNIVAL_options(solver = "cplex")
  carnival_options$solverPath <- "/Applications/CPLEX_Studio2211/cplex/bin/x86-64_osx/cplex"
  carnival_options$outputFolder <- path
  carnival_options$workdir <- path
  carnival_options$threads <- 4
  carnival_options$clonelog <- -1
  carnival_options$poolIntensity <- 2
  carnival_options$mipGap <- 0.2
  carnival_options
}

get_cosmos_data <- function() {
  env <- rlang::new_environment()
  x <- data(package = "cosmosR")[["results"]][, "Item"]
  data(list = x, package = "cosmosR", envir = env)
  mget(x, envir = env)
}

preprocess <- function(
    route = c("forward", "reverse"),
    network,
    signaling,
    metabolites,
    transcripts,
    diff_exp_threshold = 1,
    maximum_network_depth = 6,
    carnival_options
){
  carnival_options$timelimit <- 3600

  network_members <- unique(c(network$source, network$target))
  signaling <- signaling[names(signaling) %in% network_members]
  metabolites <- metabolites[names(metabolites) %in% network_members]

  route <- match.arg(route)
  if (route == "forward") {
    fun <- cosmosR::preprocess_COSMOS_signaling_to_metabolism
    remove_nodes <- TRUE
  } else if (route == "reverse") {
    fun <- cosmosR::preprocess_COSMOS_metabolism_to_signaling
    remove_nodes <- TRUE
  }

  args <-
    list(
      meta_network = network,
      signaling_data = signaling,
      metabolic_data = metabolites,
      diff_expression_data = transcripts,
      diff_exp_threshold = diff_exp_threshold,
      maximum_network_depth = maximum_network_depth,
      remove_unexpressed_nodes = remove_nodes,
      CARNIVAL_options = carnival_options
    )

  do.call(fun, args)
}

run_cosmos <- function(
    route = c("forward", "reverse"),
    data,
    carnival_options,
    timelimit = 14400
){
  carnival_options$timelimit <- timelimit

  if (route == "forward") {
    fun <- cosmosR::run_COSMOS_signaling_to_metabolism
  } else if (route == "reverse") {
    fun <- cosmosR::run_COSMOS_metabolism_to_signaling
  }

  args <- list(data = data, CARNIVAL_options = carnival_options)

  do.call(fun, args)
}

combine_cosmos <- function(forward, reverse) {
  full_sif <- as.data.frame(rbind(forward[[1]], reverse[[1]]))
  # full_sif <- full_sif[full_sif$Weight > 0,]
  full_attributes <- as.data.frame(rbind(forward[[2]], reverse[[2]]))

  full_sif <- unique(full_sif)
  full_attributes <- unique(full_attributes)
  list(sif = full_sif, attributes = full_attributes)
}

# data --------------------------------------------------------------------

format_tf <- function(df) {
  df |>
    dplyr::filter(adj.P.Val <= 0.05) |>
    dplyr::select(tf, activity = t) |>
    tibble::deframe()
}

format_metab <- function(df) {
  out <-
    df |>
    dplyr::filter(adj.P.Val <= 0.05) |>
    dplyr::select(hmdb, activity = t) |>
    tibble::deframe() |>
    cosmosR::prepare_metab_inputs(compartment_codes = c("c", "m"))
  out[["Metab__HMDB0000190_e"]] <- 10
  out
}

format_deg <- function(df) {
  df |>
    dplyr::filter(!is.na(symbol) & symbol != "") |>
    dplyr::filter(abs(stat) == max(abs(stat)), .by = "symbol") |>
    dplyr::filter(padj < 0.05) |>
    dplyr::select(symbol, stat) |>
    tibble::deframe()
}
