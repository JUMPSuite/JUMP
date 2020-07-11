import JUMPsl.spectral_data as sd, JUMPutil.job_manager as jm, numpy as np, sys, build_mass_window

lib = sd.CSRSpectralDataReader( sys.argv[2] )
delta = 1e-2
fmasses = np.arange(lib.minMZ(),lib.maxMZ()+delta,delta)

flib = sd.CSRFilterDataWriter( sys.argv[1] )
flib.initialize_masses( fmasses, build_mass_window.nsamples )
flib.close()

blksz = 16

with jm.JobManager() as job_manager:
    cmds = []
    for m in fmasses:
        cmds.append( 'python build_mass_window.py {} {} {} {}'.format(sys.argv[1],m,
                                                                      sys.argv[2],sys.argv[3]))
        if not len(cmds) < blksz:
            job_manager.add_job('JUMP_SPECTRAL_DATABASE', str.join(' ; ',cmds))
            cmds = []

    if len(cmds) > 0:
        job_manager.add_job('JUMP_SPECTRAL_DATABASE', str.join(' ; ',cmds))


