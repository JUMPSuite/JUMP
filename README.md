## Introduction
JUMP (JUmbo Mass spectrometry-based Proteomics tool) is a software suite for 
processing mass spectrometry-based proteomics data including peptide/protein identification,
filtering of identified peptide-spectrum matches (PSMs), quantification of peptides/proteins and so forth.

## Basic installation
A basic install is sufficient for multicore laptops, desktops,
workstations and servers.  For a basic install, use the `bootstrap.sh`
script provided in the repo.  Then, 

1. Place the JUMP distribution source in the desired location (call
this `<path to JUMP>`)
1. Change your working directory to the top level of the JUMP install
and run ` bashbootstrap.sh`

### Obtaining JUMP source
You can obtain the latest version of JUMP from git; simple clone the
git repository:

```
    git clone https://github.com/JUMPSuite/JUMP.git JUMP
```

in the directory _where you would like JUMP to be installed_ (call this directory <path to JUMP>).  Note
that JUMP does not support out-of-place installs; the JUMP git
repository _is_ the entire installation.  History of JUMP releases is
provided by git tags.

### Bootstrapping
To get dependencies installed, we recommend Conda, and we have
provided a bootstrapping script `bootstrap.sh` that downloads all
dependencies and installs them alongside JUMP.  Execute

```
    ./bootstrap.sh
```

and it will create a new directory `conda` in the current working
directory.  The directory `conda` will contain the conda
environment to be used by JUMP.  JUMP will be set up to use the PERL
and python interpreters in that environment.

Once `bootstrap.sh` is finished, the JUMP suite will be configured.
JUMP components then are available via the executable script `<path to
JUMP>/JUMP/bin/jump`.  For convenience, you may simply append `<path
to JUMP>/JUMP/bin` to your `PATH` environment variable as
`PATH=$PATH:<path to JUMP>/JUMP/bin`.  You will not need to activate
the conda environment `bootstrap.sh` creates. You are now ready to run
JUMP!

*You may use existing installations of PERL, python and R for JUMP,
 but this is not recommended.  PERL 6 is not supported*

### How to run JUMP
Please see [the link](manual.md) for JUMP manual.
Sample datasets (mzXML files for TMT-11plex and 16plex), parameter files and 
database files can be downloaded at https://www.stjuderesearch.org/site/lab/peng

### Installation on a HPC cluster
JUMP can be configured to interact with an HPC job manager on a
compute cluster.  Extra command-line arguments are necessary to
bootstrap.sh.  An example for an LSF cluster with no "service nodes" is

```
    ./bootstrap.sh CLUSTER=1 \
        DEFAULT_BATCH_CMD='bsub -K -P peng_grp -q standard' \
        COMPUTE_ON_LOGIN_NODE=0 \
	BATCH_DISPATCH_LAG=.01 \
	MAX_DISPATCH_WORKER_PROCS=128

```
The extra arguments to `bootstrap.sh` explained:

| Argument | Meaning |
| :--- | :--- |
| `CLUSTER` | 1 for cluster installation, 0 for local installation |
| `DEFAULT_BATCH_CMD` | The command used to submit scripts *from the command line* (see below for more information) |
| `COMPUTE_ON_LOGIN_NODE` | Set to 1 if your cluster has "service nodes" from which jobs are launched, 0 otherwise | 
| `BATCH_DISPATCH_LAG` | Lag between batch job submissions to avoid overwhelming the job manager | 
| `MAX_DISPATCH_WORKER_PROCS` | Throttling factor to prevent JUMP from spawning more than `MAX_DISPATCH_WORKER_PROCS` search jobs; useful for ensuring fairness in job priority or to prevent overuse of resources on service nodes |

#### Batch command configuration
When `CLUSTER=1`, JUMP will use your job manager (e.g. Slurm, SGE,
LSF) to submit and monitor jobs on your compute cluster.  Your job
manager must allow jobs defined as throwaway command line parameters.
For example, you must be able to do

```
$ bsub -P peng_grp -q standard -K "perl myscript1.pl && perl myscript2.pl"
```  

The command strings you create *must* create job submissions that
block until the batch job completes.  Use the following flags, based
on your job manager:

| Job Manager | Flag |
| :--- | :--- |
| SGE | `-synch` | 
| LSF | `-K` | 
| Slurm | `-W` |

#### JUMP, batch commands, and "tool types" The various stages of
JUMP's pipeline have differing compute requirements; some require more
memory, some may have longer wall clock times.  You can create
customized job submission strings for each of these "tool types" to
allow for more memory or longer run times. If you do not create a
customized job submission command for a tool type, JUMP will just use
the `DEFAULT_BATCH_CMD`.  The available "tool types" are:

| Tool type | Explanation |
| :--- | :--- |
| `JUMP_DATABASE_BATCH_CMD` | Components used for database construction.  Memory and wall clock time intensive. |
| `JUMP_SEARCH_BATCH_CMD` | Components used to run the job management and aggregation front-end of JUMP search.  Memory intensive when `COMPUTE_ON_LOGIN_NODE=0`, these jobs must persist until the entire mzXML input is searched. |
| `JUMP_QUANTIFICATION_BATCH_CMD` | Components used for quantification.  Memory intensive. |
| `JUMP_FILTER_BATCH_CMD` | Components used for filter.  Very memory intensive. |
| `JUMP_LOCALIZATION_BATCH_CMD` | Components used for site localization.  Very memory intensive. |
| `RUNSEARCH_SHELL_BATCH_CMD` | Sub-components used for searching individual spectra.  Neither compute, memory or time intensive. | 

