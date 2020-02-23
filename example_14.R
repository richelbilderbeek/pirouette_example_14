# Works under Linux and MacOS only

library(pirouette)
suppressMessages(library(ggplot2))

root_folder <- getwd()
example_no <- 14
rng_seed <- 314
example_folder <- file.path(root_folder, paste0("example_", example_no, "_", rng_seed))
dir.create(example_folder, showWarnings = FALSE, recursive = TRUE)
setwd(example_folder)
set.seed(rng_seed)
testit::assert(is_beast2_installed())

phylogeny  <- ape::read.tree(
  text = "(((A:8, B:8):1, C:9):1, ((D:8, E:8):1, F:9):1);"
)

alignment_params <- create_alignment_params(
  sim_tral_fun = get_sim_tral_with_std_nsm_fun(
    mutation_rate = 0.1
  ),
  root_sequence = create_blocked_dna(length = 1000),
  rng_seed = rng_seed
)

# JC69, strict, Yule
generative_experiment <- create_gen_experiment()
check_experiment(generative_experiment)

# All non-birth-death tree priors
candidate_experiments <- create_all_bd_experiments(
  exclude_model = generative_experiment$inference_model
)
check_experiments(candidate_experiments)

experiments <- c(list(generative_experiment), candidate_experiments)

# Set the RNG seed
for (i in seq_along(experiments)) {
  experiments[[i]]$beast2_options$rng_seed <- rng_seed
}

check_experiments(experiments)

# Shorter on Travis
if (is_on_travis()) {
  for (i in seq_along(experiments)) {
    experiments[[i]]$inference_model$mcmc$chain_length <- 3000
    experiments[[i]]$inference_model$mcmc$store_every <- 1000
    experiments[[i]]$est_evidence_mcmc$chain_length <- 3000
    experiments[[i]]$est_evidence_mcmc$store_every <- 1000
    experiments[[i]]$est_evidence_mcmc$epsilon <- 100.0
  }
}

# Log-transformed nLTT function
lt_nltt <- function(tree, trees) {
  nLTT::nltts_diff(
    tree = tree,
    trees = trees,
    log_transform = TRUE
  )
}

pir_params <- create_pir_params(
  alignment_params = alignment_params,
  experiments = experiments,
  twinning_params = create_twinning_params(
    rng_seed_twin_tree = rng_seed,
    rng_seed_twin_alignment = rng_seed
  ),
  error_measure_params = create_error_measure_params(
    error_fun = lt_nltt
  )
)

rm_pir_param_files(pir_params)

errors <- pir_run(
  phylogeny,
  pir_params = pir_params
)

utils::write.csv(
  x = errors,
  file = file.path(example_folder, "errors.csv"),
  row.names = FALSE
)

pir_plot(errors) +
  ggsave(file.path(example_folder, "errors.png"))

pir_to_pics(
  phylogeny = phylogeny,
  pir_params = pir_params,
  folder = example_folder
)

pir_to_tables(
  pir_params = pir_params,
  folder = example_folder
)
