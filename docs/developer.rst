=====================
JUMP Developer Manual
=====================

The JUMP configuration system
-----------------------------
JUMP has a configuration system for storing settings that control
*how* JUMP runs its jobs (not *what* the jobs do).  The configuration
system is the foundation on which JUMP's job manager types build.  

JUMP's configuration system allows for global (i.e. site-level)
settings that multiple users share.  Users can also override the
site-level settings with a file in their home directory:
``~/.jump/config``.  This file consists of two-column, whitespace-separated text formatted as

``setting_name setting_value``

An example configuration file is below:

.. code-block:: 
    :linenos:

    default_batch_cmd bsub -P hprc -K -R "rusage[mem=8192]"
    max_search_worker_procs 256

The above configuration file sets the ``default_batch_cmd`` setting to
``bsub -P hprc -K -R "rusage[mem=8192]"`` and
``max_search_worker_procs`` to ``256``.

JUMP's configuration system consists of these parts:
  * ``Spiders::Config.pm``, which loads and stores settings, and resolves
     conflicts between site- and user-level settings
  *  The command-line interface ``jump -config`` which can print and
     store settings to the terminal
  *  ``JUMPutils.config``, which  provides a python interface into
     the JUMP config system.
  
The JUMP job manager system
-----------------------------
JUMP uses fork-join parallelism to speed up search and database 