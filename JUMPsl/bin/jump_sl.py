import JUMPsl.binned_search as bs, JUMPsl.spectral_data as sd, JUMPutil.params as pa, JUMPutil.config as cf, JUMPutil.job_manager as jm, pickle

blksz = 96

def get_output_extension( output_format ):
    if 'pickle' == output_format:
        return 'ppd'

def decompose_domain( args ):
    q = sd.DTASReader( args.input )
    masses = list(sorted(q.spectra(),key=lambda x: q.prec_mass(x)))
    with jm.JobManager() as manager:
        for i in range(0,len(masses),blksz):
            fmt = str.join(',',[str(m) for m in masses[i:i+blksz]])
            manager.add_job( 'JUMP_SPECTRAL_SEARCH',
                             'python jump_sl.py --runsearch-shell=True --idxs={} --output-format={} --output={} {} {}'.format(
                    fmt,
                    args.output_format,
                    args.input + '.rank{}.ppd'.format(i),
                    args.params_file,
                    args.input) )
        

def runsearch_shell( args ):
    params = pa.Params( args.params_file )
    q = sd.DTASReader( args.input )

    rdr = sd.CSRSpectralDataReader( params.spectral_library )
    if 'binned' == params.spectral_library_method:
        s = bs.BlockEagerBinnedSearch( rdr, q.maxMZ(), params.bin_size, params.mass_window, params.n_hits, search_subset=args.idxs )
    
    if 'pickle' == args.output_format:
        pickle.dump( list(s(q)), open(args.output,'wb') )
    else:
        raise TypeError('do not know how to product output of type {}'.format(args.output))

if __name__ == '__main__':
    import sys, argparse

    parser = argparse.ArgumentParser()
    parser.add_argument('--output-format',default='pickle')
    parser.add_argument('--output',default='')
    parser.add_argument('--idxs',default=None)
    parser.add_argument('--runsearch-shell',default=True,type=bool)
    parser.add_argument('params_file')
    parser.add_argument('input')

    args = parser.parse_args()
    if not None == args.idxs:
        args.idxs = [int(i) for i in args.idxs.split(',')]

    if '' == args.output:
        args.output = args.input.replace('dtas',get_output_extension(args.output_format))

    if args.runsearch_shell:
        runsearch_shell( args )
    else:
        decompose_domain( args )

