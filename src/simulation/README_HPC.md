# HPC protocol for the whole-epidemic simulation

One Slurm array task runs one complete Monte Carlo replicate across all four
blocks and all three treatments. Eight block-treatment scenarios run in
parallel within a task. Forward transitions and backward fits within a
scenario remain sequential, and threaded math libraries are restricted to one
thread to prevent oversubscription.

The default Slurm request is 8 CPU cores, 20 GB RAM, and 2 hours per task. The
submission helper defaults to 12 concurrent tasks, for a maximum request of 96
cores and 240 GB. These values deliberately remain well below the `share`
partition limits.

Submit a short pilot with:

```bash
bash src/simulation/submit_whole_epidemic_simulations.sh 1 10 4
```

After checking runtimes and memory use, submit simulations 11 through 1000:

```bash
bash src/simulation/submit_whole_epidemic_simulations.sh 11 1000 12
```

A task writes exactly one atomic RDS file. A result enters `complete/` only if
all 12 block-treatment scenarios, 48 forward transitions, and 416 backward
fits are returned and all required optimizers converge. Otherwise, the entire
replicate is written to `failed/` and the Slurm task exits unsuccessfully. No
partial replicate is included in combined scientific results.

Hessian invertibility is recorded as a diagnostic but is not a completion
requirement. This lets the study quantify weak identifiability without silently
discarding those simulations.

Rerunning the same simulation ID reproduces its deterministic seed. A
persistently failing replicate can be replaced with a new unused simulation
ID above the planned range. Failure rates and reasons should still be retained
and reported.

Combine 1,000 complete results only after all replacements are available:

```bash
Rscript --vanilla src/simulation/combine_hpc_simulations.R 1000
```

The combine step requires exactly the requested number of complete replicate
files. It also writes `all_attempt_status.csv`, which retains the status and
failure reasons for excluded attempts.
