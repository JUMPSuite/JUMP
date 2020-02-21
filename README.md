## Basic installation
There are three steps to a basic install:

1. Place the JUMP distribution source in the desired location (call
1. Obtain a PERL and python executable and install all the requisite
dependencies
this `<path to JUMP>`)
1. Change your working directory to the top level of the JUMP install
and run `Makefile.PL`

### Obtaining JUMP source
You can obtain the latest version of JUMP from git; simple clone the
git repository:

```
    git clone https://github.com/JUMPSuite/JUMP_v1.13.1.git JUMPsuite
```

in the directory _where you would like JUMP to be installed_ (call
this directory <path to JUMP>).  Note
that JUMP does not support out-of-place installs; the JUMP git
repository _is_ the entire installation.  History of JUMP releases is
provided by git tags.

### Installing PERL, python and dependencies
To get dependencies installed, we recommend Conda, and we have
provided a bootstrapping script `bootstrap_conda.sh`; just execute

```
    sh bootstrap_conda.sh
```

and it will create a new directory `conda/jump` in the current working
directory.  The directory `conda/jump` will contain the conda
environment to be used by JUMP.  One can then activate that
environment with

```
    conda activate <path to directory where you executed boostrap_conda.sh>/conda/jump
```

Once the PERL and python dependencies for JUMP have been installed,
one is ready to perform the configuration and install of JUMP.

### Configuring JUMP
To complete installation of JUMP, one must execute `Makefile.PL`; that
script is in the top-level of the JUMP distribution.  Therefore:

1. Activate the conda environment with
`conda activate <path to directory where you executed boostrap_conda.sh>/conda/jump`
1. `cd <path to JUMP>/JUMPsuite`
1. `perl Makefile.PL MAX_SEARCH_WORKER_PROCS=<num cores>`

For the last step, set `<num cores>` to be the maximum number of cores you
would like JUMP to use.  For most purposes, you can set `<num cores>`
to be the total number of cores on your machine; if your workstation
has 20 cores, then execute `perl Makefile.PL MAX_SEARCH_WORKER_PROCS=20`

The JUMP suite will then be configured.  JUMP
components then are available via the executable script `<path to
JUMP>/JUMP/bin/jump`.  For convenience, you may simply append `<path to
JUMP>/JUMP/bin` to your `PATH` environment variable as `PATH=$PATH:<path to
JUMP>/JUMP/bin`.  You are now ready to run JUMP!
