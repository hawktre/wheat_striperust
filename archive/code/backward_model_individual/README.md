# Candidate-specific parameter sensitivity analysis

This directory contains the supplementary analysis in which each candidate
source set receives its own parameter vector during the M-step. It is retained
to evaluate sensitivity to that modeling choice.

The primary analysis is in `Code/R/source_detection`. There, all candidates
share one parameter vector within each block × treatment × transition, and one
latent source set is inferred using evidence across the complete study.

The scripts here predate the primary study-level implementation and should not
be interpreted as its production runner without reconciliation with the main
likelihood and posterior definitions.
