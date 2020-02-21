## Basic installation
There are two steps to a basic install:

1. Obtain a PERL and python executable and install all the requisite
dependencies
1. Place the JUMP distribution source in the desired location (call
this `<path to JUMP>`) and run `Makefile.PL`

To get dependencies installed, we recommend Conda, and we have
provided a bootstrapping script `bootstrap_conda.sh`; just execute

```
    sh bootstrap_conda.sh
```

and it will create a new conda environment in the directory
`conda/jump` in the current working directory.

Once the PERL and python dependencies for JUMP have been installed,
one is ready to perform a basic install of JUMP.

1. Ensure that the PERL and python executables you are using are the
ones in which all the JUMP dependencies have been installed
1. Execute `perl Makefile.PL <options>`

The JUMP suite will then be configured.  Once, complete, JUMP
components are available via the executable script `<path to
JUMP>/JUMP/bin/jump`.  For convenience, you may simply append `<path to
JUMP>/JUMP/bin` to your `PATH` environment variable as `PATH=$PATH:<path to
JUMP>/JUMP/bin`.  You are now ready to run JUMP!

## Configuring JUMP to use multiple cores
JUMP contains configuration parameters that control how many cores
JUMP will use.  This parameter can be set at installation time or
afterwards.

To configure at installation time, do

```
    perl Makefile.PL MAX_SEARCH_WORKER_PROCS=<num cores>
```

This parameter can be overridden at run time with JUMP's user-level
configuration system.  User-level parameters will always override
install-time parameter settings.  JUMP will search for a file in
`$HOME/.jump/config`, and read whitespace-separated key-value pairs
from that file.  

To set the `max_search_worker_procs`
parameter at run time, add

```
    max_search_worker_procs <num cores>
```

to the file `$HOME/.jump/config`, creating the directory and file if
they do not already exist.  For example, if 20 cores are available,
then add

```
    max_search_worker_procs 20
```

to `$HOME/.jump/config`.