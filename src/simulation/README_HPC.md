# HPC protocol for the whole-epidemic simulation

The study uses 100 Slurm array jobs with 10 Monte Carlo simulations per job,
for 1,000 simulations total. Each job maintains a pool of up to ten
single-core simulations. When one simulation finishes, its free core
immediately starts the next simulation assigned to that job. The twelve
block-treatment scenarios within one simulation run sequentially; forward
transitions and backward fits within a scenario also remain sequential.
Threaded math libraries are restricted to one thread.

The default Slurm request is 10 CPU cores, 5 GB RAM, and 2 hours per job.
Five local replicates required an estimated 3.3--5.7 minutes with eight workers,
and a production-style local run required 5.6 minutes when parallelism was
within a simulation. The cross-simulation worker pool avoids waiting at every
simulation boundary; the two-hour request allows for slower nodes, long EM
fits, and Monte Carlo runtime variation. A full 10-simulation HPC pilot used
2.08 GB of resident memory, so 5 GB retains more than twofold headroom.

The submission helper submits all 100 jobs in one array. Its third argument
throttles how many jobs may run concurrently. With 10 cores and 5 GB per job,
the partition limits would permit at most 102 simultaneous jobs by CPU and 150
by memory. The helper defaults to 60, which requests at most 600 cores and
300 GB. Reduce this throttle when other jobs are using the same per-user
allocation.

A complete 100-job array reserves at most 83.3 CPU-days and 41.7 GB-days when
the full two-hour requests are charged, below the documented 1,024 CPU-day and
3,072 GB-day limits.

The cluster's effective limits may also depend on the user's account and QOS.
They can be rechecked on the login node with:

```bash
scontrol show partition share
scontrol show config | grep MaxArraySize
sacctmgr show assoc where user="$USER" format=User,Account,Partition,QOS,MaxJobs,MaxSubmit
sacctmgr show qos format=Name,MaxJobsPU,MaxSubmitJobsPU,MaxTRESPU
```

Submit one pilot job containing only simulations 1--2 with:

```bash
bash src/simulation/submit_whole_epidemic_simulations.sh 1 1 1 2
```

The fourth argument temporarily changes the simulations-per-job value from 10
to 2. After checking the pilot log and results, submit the complete production
array with the normal 10-simulation batches:

```bash
bash src/simulation/submit_whole_epidemic_simulations.sh 1 100 60 10
```

The first production job will detect and skip simulations 1 and 2 if their
pilot results are already complete.

Each simulation writes exactly one atomic RDS file. A result enters `complete/` only if
all 12 block-treatment scenarios, 48 forward transitions, and 416 backward
fits are returned, every forward optimizer converges, every backward EM fit
converges, and every backward fit has a nondecreasing observed likelihood.
Otherwise, the entire replicate is written to `failed/` and the Slurm task
exits unsuccessfully. No partial replicate is included in combined scientific
results.

An intermediate BFGS call that reports a nonzero code is retained in
`any_m_step_failure` as a diagnostic, but does not invalidate a fit when later
M-steps and the overall EM algorithm converge monotonically. Requiring no
transient BFGS codes would reject all five completed local pilot replicates.

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

After the pilot array finishes, inspect resource use with (replace the job ID):

```bash
sacct -j JOB_ID --units=G \
  --format=JobID,State,Elapsed,AllocCPUS,ReqMem,MaxRSS,ExitCode
```

The combine step requires exactly the requested number of complete replicate
files. It also writes `all_attempt_status.csv`, which retains the status and
failure reasons for excluded attempts.
