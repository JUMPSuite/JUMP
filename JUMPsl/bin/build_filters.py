import JUMPsl.spectral_data as sd, JUMPutil.job_manager as jm, numpy as np, sys

lib = sd.CSRSpectralDataReader( sys.argv[1] )
fmasses = np.linspace(lib.minMZ(),lib.maxMZ(),5)

flib = sd.CSRFilterDataWriter( sys.argv[2] )
flib.initialize_masses( fmasses, 61 )
flib.close()

with jm.JobManager() as job_manager:
    for m in fmasses:
        job_manager.add_job( 'JUMP_SPECTRAL_DATABASE', 'python build_mass_window.py {} {}'.format(sys.argv[2],m) )



