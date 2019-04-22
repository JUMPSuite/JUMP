.. JUMP documentation master file, created by
   sphinx-quickstart on Wed Feb 14 15:44:52 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


.. toctree::
   :maxdepth 2:

.. author :: St. Jude CRH

**************************************************
Top Level 
**************************************************


Introduction
**************************************************
With the version 1.13.1 release, several of the key JUMP pipeline
tools have been migrated to St. Jude's central HPC resources.  This
allows for faster processing of shotgun mass spectrometry proteomics
data, with tools to better support version management and, through
that, reproducibility of results.

The purpose of this document is to discuss differences in use of JUMP
between the Structural Biology compute resources and the HPC
implementation.  These differences are discussed both in terms of
basic JUMP uses, and expert-level and development JUMP uses.  We
discuss version and module management; job dispatching, monitoring and
cleanup; and file storage.

Getting started: Basic JUMP uses
**************************************************
This section discusses basic uses of the JUMP pipeline tools; by the
end of this section, the reader should have a working knowledge of how
to run JUMP analyses on the St. Jude HPC system, and use the
``module`` tools to manage software versions.

The ``module`` tools
--------------------------------------------------
The JUMP installation on the St Jude HPC uses the ``module`` tools to
manage versions and library dependencies.  For example, you can switch
between different versions of JUMP to compare differences or reproduce
results obtained in the past.

To see what versions of jump are available,
use the ``module avail`` command::

  module avail jump

You will see numbered versions, and may see an alpha version, which is
for pre-release testing.

.. attention:: Do not use the ``jump/alpha`` release for production
	       runs. It may contain bugs.
  
To load a version of JUMP, use the ``module load`` command.  For
example, to load the bleeding-edge alpha version, use::
  
  module load jump/alpha

There is also a ``module unload`` command that simply unloads a
module.
  
To see what version of JUMP you currently have loaded,
use::
  
  module list

To switch between versions of JUMP, use the ``module swap`` command;
for example, to switch from the alpha version to version 1.13.1,
use::
  
  module swap jump/alpha jump/1.13.1

The above command is equivalent to unloading the alpha version and
loading version 1.13.1::
  
  module unload jump/alpha
  module load jump/1.13.1

It is a good idea to create production scripts that incorporate module
commands.  See `Software Carpentry` for more details.

The ``module`` command manipulates your environment variables to load
all the requisite dependencies of JUMP.  If you wish to use ``module``
to manage PERL or R, refer to `Expert-level use of the JUMP search tool` for more
detailed information on how these will interact with JUMP.

The top-level JUMP tool
--------------------------------------------------
JUMP on the HPC has the same top-level ``jump`` command that exists on
the proteomics compute cluster.  To run JUMP search, simply use::

  jump -s -p jump.params data.mzXML

where ``jump.params`` is your parameter file and ``data.mzXML`` is
your MZXML file.

Likewise,
for JUMP filtering, use::

  jump -f jump.params

For
JUMP quantification, use::

  jump -q jump.params

For
JUMP localization, use::

  jump -l jump.params

For
JUMP database, use::

  jump -d jump.params  

Automatic log file generation
++++++++++++++++++++++++++++++++++++++++++++++++++
The individual tools will automatically create log files that contain
their standard output and standard error.  The log files are created
in the directory in which the jump tool is invoked.  They are given
names based on the names of the input files that are passed to JUMP.

JUMP search creates log file names based on the mzXML or RAW input
file names; the input file names are concatenated together, ".out"
is appended to that concatenation, and the local date and time are
prepended to that.  For example, if one had invoked::

  jump -s jump.params data1.raw data2.mxXML

on 4/2/2016 at 13:45:47, then the log file will be named ``4-2-2016_13:45:47_data1-data2.out``.

The remaining JUMP tools will create a log file based on the name of
the parameter file.  If ``jump.params`` is passed in, then the log
file will be named ``4-2-2016_13:45:47_jump.params.out``.

.. attention:: The log files produced by JUMP may be large.  The user
	       should excercise care and delete those that are not
	       necessary in order to avoid storage chargebacks.

The individual JUMP tools
--------------------------------------------------
For each option of the top-level JUMP command, there is an individual
jump command.  ``jump_sj.pl`` is the same as ``jump -s``,
``jump_f.pl`` is the same as ``jump -f``, ``jump_q.pl`` is the same
as ``jump -q``, ``jump_d.pl`` is the same as ``jump -d``, and
``jump_l.pl`` is the same as ``jump -l.``


Expert-level JUMP usage
**************************************************
Some of the JUMP tools allow for customization beyond the basic use
cases described in `Getting started: Basic JUMP uses`.  These extra
options are for controlling how JUMP allocates and uses HPC resources;
they do not change the *results* that JUMP produces.  In this section,
we will only provide examples using the individual tools, as the
top-level tool will accept the same command-line options.

The ``module`` tool and environment configuration
--------------------------------------------------
JUMP uses the ``module`` tool to manage versions and dependencies.  The
`module` tools manipulate environment variables to do this.  Version
1.13.1 of JUMP no longer uses ``FindBin`` to determine the location of
JUMP's Spiders modules to provide more flexibility in how JUMP may be
installed.  Instead, the module tools will set the environment variables
``JUMP_SJ_LIB``, ``JUMP_F_LIB``, and ``JUMP_Q_LIB`` to point to the
directories containing the Spiders modules for JUMP search, filter and
quantification, respectively.  If those variables are already defined
in one's environment prior to invoking ``module load jump/<version>``,
then the old values are saved and restored upon unloading the jump module.

The JUMP module file also loads the correct versions of PERL and R.
It has been witnessed that JUMP is sensitive to newer versions of
PERL, specifically, PERL 5.20.1 appears to produce non-deterministic
results for JUMP search.  Care has been taken to flag conflicts so
that JUMP modules cannot be loaded if an inappropriate version of PERL
has been loaded first.  However, one can always load JUMP *first* and
then load an inappropriate version of PERL.

.. attention:: Take care *not* to load ``perl/5.20.1`` when using
	       JUMP.  You may get
	       non-deterministic---i.e. non-reproducible---results.

Expert-level use of the JUMP search tool
--------------------------------------------------
JUMP search is one of the most computationally expensive tool, and
therefore has expert-level options for controlling how it runs on
the HPC cluster.  The new feature introduced for JUMP 1.13.1 on the
HPC lets it dispatch jobs on the HPC in several ways to cope with
resource use limitations.

JUMP search and job dispatch
++++++++++++++++++++++++++++++++++++++++++++++++++
The HPC login nodes have resource usage limitations that prohibit
long-running jobs or jobs that have large memory usage.  Therefore, a
large-scale JUMP search analysis may not run successfully on an HPC
login node.  To cope with that, ``jump_sj.pl`` has the ``--dispatch``
flag.  JUMP search can be run interactively on a compute node, or run
non-interactively on a compute node, or run interactively on the login
node.

The default behavior of JUMP search is to run interactively on a
compute node.  Simply using::

  jump_sj.pl -p jump.params data.mzXML

will work, but one may explicitly specify interactive compute node
dispatch as::

  jump_sj.pl -p jump.params data.mzXML --dispatch=batch-interactive

Note that the ``--dispatch=batch-interactive`` flag may be put at any
place in the argument list; it does not need to be the last argument.

In some cases, one might want to redirect output from JUMP search to a
file.  For example, one might want this when one's connection to the
HPC may be interrupted, such as when one is logging in via a VPN.  To
get a non-interactive batch dispatch, use::

  jump_sj.pl -p jump.params data.mzXML --dispatch=batch
  
The output will be put in a file named ``data.out``.

Finally, one can run JUMP search directly on the login node.  This is
useful for debugging JUMP, for example.  To get login node dispatch,
use::

    jump_sj.pl -p jump.params data.mzXML --dispatch=localhost
    
.. attention:: The login nodes prohibit CPU and memory intensive
	       jobs.  Use the login node mode *only* when the data set
	       is small.  Long-running or memory-intense jobs slow
	       down the login node for all users and may be
	       spontaneously terminated.
	       
 
JUMP search and HPC queues
++++++++++++++++++++++++++++++++++++++++++++++++++
By default, JUMP will dispatch to the ``heavy_io`` queue, which is
optimzed for input/output operations.  During periods of high load,
this may cause a nontivial wait time for dispatch.  One can specify
an alternate queue; for example, using::

  jump_sj.pl -p jump.params data.mzXML --queue=normal

will dispatch to the ``normal`` queue.  The HPC cluster has several
queues, which one may view with ``bqueues``.  Notice that some queues
require special access.  We recommend using either the ``normal`` or
``heavy_io`` queues.


Expert-level use of the JUMP localization and filter tools
--------------------------------------------------
JUMP filter and localization both may require several gigabytes of
memory to run. To manage resource usage, both JUMP filter and
localization present advanced options to control the resources they
request.  The first option is the ``--queue`` option, analagous to
JUMP search's ``--queue`` option, and the second is the ``--mem``
option.  By default, JUMP filter and localization dispatch to the
``normal`` queue, and request 200G of memory.

The ``mem`` option requires use of the ``queue`` option, and it is
incumbent on the user to make sure that the queue selected for
dispatch can also support the memory request.  For example, if one
wants to request 2TB of memory, one must use the ``large_mem`` queue
as::

  jump_f.pl jump.params --queue=large_mem --mem=2000000

For JUMP localization, the invocation is the same::

  jump_l.pl jump.params --queue=large_mem --mem=2000000

  
Software Carpentry
**************************************************
The ``module`` tools combined with version control [#f1]_ give some basic
control over versioning that can be useful for reproducing results
obtained in the past.  This can be useful, for example, when one has
submitted a paper and revisions are requested.  The principles we
discuss fall under the rubric of "software carpentry," and the
interested reader can find more at on the `software carpentry website
<https://software-carpentry.org/>`.

To put some files in a directory ``~/analysis`` under version control,
use::

  cd ~/analysis
  git init

We presume that the directory ``~/analyis`` already has the parameter
and MZXML in it.  Then create a shell script that just loads the right
JUMP module and executes a JUMP workflow.  The shell script---we name
it ``run-jump.sh``---might look
like

.. code-block:: shell
		
   #!/bin/bash
   
   module load jump/1.13.1
   jump_sj.pl -p jump.params data.mzXML

Once the parameters and data produce the result you want, do::

  git add run-jump.sh jump.params data.mzXML
  git commit -m 'Meaningful message that will jog your memory'
  git tag 'run-for-fancy-journal'

to get a snapshot of the parameters and data. Then you can reproduce
these results at any later time with::

  git checkout 'run-for-fancy-journal' -b recovery

which recovers the snapshots of all your files and puts them in a
branch called "recovery."  Then invoking
``run-jump.sh`` will reproduce your old results.  You can recover your
latest version of parameters and data with::

  git checkout master

For more powerful use of GIT, contact the RIS team.

.. rubric:: Footnotes

.. [#f1] For version control, we use GIT in our examples, based on its
	 ease of setup and use.

