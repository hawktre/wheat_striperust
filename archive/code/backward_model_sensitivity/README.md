# Initialization sensitivity analysis

This directory contains the supplementary analysis of sensitivity to starting
values, particularly the initial value of the nonlinear distance parameter
`kappa`.

The primary analysis is in `Code/R/source_detection`. Its forward model uses
method-of-moments estimates over a grid of initial kappa values, retains the
best finite converged fit for each transition, and passes those estimates to
the study-level EM model as warm starts.

The scripts here predate the primary study-level implementation and remain
separate so alternative-start results can be reported in the supplement.
