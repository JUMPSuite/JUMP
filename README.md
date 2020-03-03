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
    git clone https://github.com/JUMPSuite/JUMP_v1.13.1.git JUMPsuite
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
    bash bootstrap.sh
```

and it will create a new directory `conda/jump` in the current working
directory.  The directory `conda/jump` will contain the conda
environment to be used by JUMP.  JUMP will be set up to use the PERL
and python interpreters in that environment.

Once `bootstrap.sh` is finished, the JUMP suite will be configured.  JUMP
components then are available via the executable script `<path to
JUMP>/JUMP/bin/jump`.  For convenience, you may simply append `<path to
JUMP>/JUMP/bin` to your `PATH` environment variable as `PATH=$PATH:<path to
JUMP>/JUMP/bin`.  You are now ready to run JUMP!
