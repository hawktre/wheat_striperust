# Readable source-detection reference

This directory contains a deliberately slow, sequential restatement of the
primary source-detection analysis. It exists for reading and quality control;
it is not the production analysis.

The reference implementation uses explicit loops for plants, candidates,
transitions, EM iterations, source counts, and analysis scenarios. It also
uses the finite-difference derivatives supplied by `optim()` instead of the
analytic gradient used by the production implementation. This makes it useful
for checking the likelihood and EM objective through an independent code
path.

The scripts contain the same substantive choices as the primary pipeline:

- empirical transition-specific zero probabilities;
- unconstrained beta, delta, and gamma;
- positive kappa and phi represented on the log scale;
- method-of-moments forward initialization over a starting-kappa grid;
- one persistent candidate source set per fitted window;
- one shared parameter vector per transition;
- fixed uniform priors over candidate sets;
- independent-transition, cumulative, and whole-study fits;
- known source counts at each spatial resolution;
- endpoint maximum-severity predictions; and
- Hungarian distance-weighted accuracy, including unequal source counts.

No production script sources these files. If the readable runner is ever
executed, it writes to `output/source_detection_readable_reference`
so it cannot overwrite the primary analysis.
