# Current whole-epidemic simulation design

The current simulation propagates disease recursively from the observed first
visit. At every later transition, all infected plants may contribute to
secondary transmission. The generating parameter vector, wind matrix, and
empirical zero probability are specific to the block, treatment, and
transition. Consequently, block is a generating scenario rather than a Monte
Carlo replicate.

Simulation is exact under the fitted conditional ZIB model. Because that
likelihood does not constrain disease severity to be monotone, a simulated
plant can decrease in severity or return to zero. The simulation records the
frequency and magnitude of these events by transition so this limitation can
be quantified rather than hidden by modifying the generating distribution.

One computational unit is one simulation replicate, block, and treatment. The
forward model is refitted to the simulated epidemic. Its transition-specific
estimates are then used to initialize source-detection fits with the known
number of occupied source groups. As in the primary analysis, one-source
experiments use resolutions 4, 8h, 8v, 16, and 64; multi-source experiments
use resolutions 4, 8h, 8v, and 16.

For every resolution, the backward model is fit independently to each
transition and cumulatively through every endpoint. The final cumulative fit
is the whole-study fit.

The saved outputs track:

- epidemic summaries by visit;
- counts, proportions, and magnitudes of non-monotone disease trajectories;
- forward parameter truth, estimates, raw-scale bias, and log-scale bias for
  kappa and phi;
- forward optimizer convergence, evaluations, selected initial kappa,
  numerical-Hessian rank, eigenvalues, condition number, invertibility,
  standard errors, Wald coverage, and elapsed time;
- backward exact, component, distance-weighted, and naive prediction metrics;
- posterior probability and rank of the true source set, entropy, and the
  effective number of candidates;
- backward EM status, iterations, monotonicity, M-step status, final changes,
  errors, and elapsed time; and
- backward parameter estimates for secondary diagnostic use.

`run_simulation_study_local.R` is intentionally sequential and defaults to one
simulation across all four blocks and all three treatments. It checkpoints
after every block-treatment scenario. An HPC array runner will be added after
this local design passes a small smoke test.
