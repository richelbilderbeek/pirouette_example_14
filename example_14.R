#
# Difference from standard: 
# - Uses LTN as error function
#
# Works under Linux and MacOS only

library(pirouette)
suppressMessages(library(ggplot2))

################################################################################
# Constants
################################################################################
is_testing <- is_on_travis()
example_no <- 14
rng_seed <- 314
folder_name <- paste0("example_", example_no, "_", rng_seed)

set.seed(rng_seed)
phylogeny  <- ape::read.tree(
  text = "(((A:8, B:8):1, C:9):1, ((D:8, E:8):1, F:9):1);"
)

pir_params <- create_std_pir_params(
  folder_name = folder_name
)

# Log-transformed nLTT function
lt_nltt <- function(tree, trees) {
  nLTT::nltts_diff(
    tree = tree,
    trees = trees,
    log_transform = TRUE
  )
}

pir_params$error_measure_params$error_fun <- lt_nltt

if (is_testing) {
  pir_params <- shorten_pir_params(pir_params)
}

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
