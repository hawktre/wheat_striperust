# Source-detection analyses

This folder contains the primary analysis. It fits one latent source set for
each block × treatment experiment across the complete study. Each transition
has its own parameter vector, shared by all candidate source sets.

The model uses a zero-inflated beta response. `beta`, `delta`, and `gamma` are
unconstrained. Positive `kappa` and `phi` are optimized as `log_kappa` and
`log_phi`. The zero probability is the empirical proportion of zeros within a
transition.

`source_detection_functions.R` contains the forward model, EM algorithm,
likelihood, gradient, posterior summaries, and supporting functions.

`run_main_analysis.R` is the pipeline entry point. It first fits the full
forward model using method-of-moments starts over a grid of initial kappa
values. Those transition-specific estimates initialize the EM model. The
pipeline uses all resolutions through 64 individual plants for one-source
inference and resolutions through 16 cells for multi-source inference.
At each resolution, the source-set size is the number of distinct cells
occupied by the known inoculated plants. Multiple inoculated plants in one
coarse cell therefore represent one true source cell at that resolution.

Rows of the transition/cumulative EM fit grid run in parallel with
`parallel::mclapply()`. By default the runner uses one fewer than the detected
number of physical cores. Set the `SOURCE_DETECTION_CORES` environment variable
to an explicit positive integer to control resource use, for example
`SOURCE_DETECTION_CORES=4 Rscript src/source_detection/run_main_analysis.R`.
Transition-specific M-step optimizations inside each worker remain sequential.

After the analysis finishes, run `plot_main_analysis.R` to create result and
diagnostic figures under `output/source_detection/plots`.
It also writes compact posterior and convergence diagnostic CSV files alongside
the primary results.

The primary output is `whole_study_posteriors.csv`. The same fitted EM model
also produces `transition_and_cumulative_posteriors.csv`: `transition` rows use
one transition at a time, while `cumulative` rows use all evidence through the
listed visit. These are diagnostic predictions; they do not refit the latent
source independently at each transition.

The primary analysis includes independent transition refits and cumulative
refits. Each transition refit uses one interval only, so its inference does not
borrow information from earlier or later intervals. The visit-2 cumulative fit uses only the
first transition, the visit-3 fit uses the first two transitions, and so on.
The final cumulative fit is the whole-study fit. Earlier fits therefore do not
borrow source responsibilities or fitted EM parameters from later visits.

`cumulative_fit_posteriors.csv` contains every candidate posterior at every
cumulative endpoint. `whole_study_posteriors.csv` is its final-visit subset.
`transition_fit_posteriors.csv` contains the independently fitted transition
results, and `all_primary_fit_posteriors.csv` combines all three fit scopes.
The older within-fit transition and cumulative decomposition is retained only
for the final fit in `transition_and_cumulative_posteriors.csv` and should be
treated as a retrospective diagnostic.

For the 4-, 8-, and 16-cell resolutions, the primary pipeline fits candidate
source-set sizes from one through four and selects the source count with the
smallest BIC within each experiment, resolution, and temporal fit. The
64-plant resolution remains restricted to one-source analyses. BIC uses the
marginal observed-data likelihood and counts the five optimized parameters
plus the empirically estimated zero-mass probability for each included
transition; candidate prior weights are fixed and are not counted as estimated
parameters. All source-count fits and the BIC-selected subsets are saved
separately.

Every fitted source count also receives a naive endpoint prediction: the
groups containing the largest maximum plant severities at the endpoint visit.
The model and naive predictions are evaluated with the same exact and
distance-weighted metrics. When predicted and true source counts differ,
precision, recall, F1, and a padded Hungarian assignment are retained so the
joint count-and-location procedure can be evaluated without dropping errors.

`source_accuracy.csv` reports exact and distance-weighted accuracy for each
cumulative refit, including the final whole-study fit. `RcppHungarian` pairs
predicted and true cells to minimize total distance. Each paired prediction is
scored as one minus its distance divided by the maximum distance from that
particular true cell to any candidate-cell centroid.

## Sensitivity analyses

The existing `backward_model_individual` analysis evaluates candidate-specific
parameter vectors. This is a supplementary sensitivity analysis because the
main model uses one shared vector per experiment-transition.

The existing `backward_model_sensitivity` analysis evaluates sensitivity to
starting kappa values. The main model instead uses the best finite converged
forward fit from its initial-kappa grid as a warm start.

`run_initialization_sensitivity.R` is the current initialization sensitivity
pipeline. It fits only the whole-study scenario and bypasses the forward-model
optimization. For each value in the initial-kappa grid, transition-specific
method-of-moments estimates are passed directly to the same shared-parameter EM
algorithm used in the primary analysis. Every start is retained, and results
are written to
`output/source_detection_initialization_sensitivity`. When the
primary outputs are available, `initialization_comparison.csv` also reports
agreement with the warm-start source prediction and the observed-likelihood
difference.

After that pipeline completes, `plot_initialization_sensitivity.R` creates
convergence, likelihood, prediction-stability, accuracy, and optimized-kappa
diagnostics. It also writes compact stability summaries alongside the model
outputs.

The Stan files are exploratory attempts and are not part of the primary or
supplementary pipelines.
