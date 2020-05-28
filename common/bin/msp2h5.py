import JUMPsl.spectral_data as sd, heapq, numpy as np, h5py, tempfile, os.path, re, bz2, JUMPutil.sorting as so


def parse_comment( comment_line ):
    return dict([(m.group(1).lower(),m.group(2)) for m in
                 re.finditer( '(\\S+)=(\\S+)', comment_line )])

def multiopen( fname ):
    if '.bz2' == fname[-4:]:
        return bz2.open(fname,'rt')
    else:
        return open(fname)

if __name__ == '__main__':
    import optparse, sys
    p = optparse.OptionParser()
    p.add_option( '-v', '--verbose', dest='verbose', action='store_true' )
    p.add_option( '-n', '--num-spectra', dest='num_spectra', type='int', default=0 )
    p.add_option( '-p', '--precursor-mass-key', dest='precursor_mass_key',
                  default='precursormz' )

    opts,args = p.parse_args()

    def get_precursor_mass( mdata ):
        if opts.precursor_mass_key.lower() in mdata:
            return float(mdata[opts.precursor_mass_key.lower()])
        else:
            return float(mdata['parent'])

    if os.path.exists(sys.argv[-1]):
        raise Exception('will not overwrite file at ' + sys.argv[-1])

    sorter = so.OutOfCoreMerger(1024)
    num_written = 0
    for txtfile in args[:-1]:
        if opts.verbose:
            if num_written > 0: print('')
            print( 'processing file {}...'.format(txtfile), end='', flush=True )
        f = multiopen(txtfile)
        for l in f:
            if re.search( '^Name:', l ):
                break

        t = l.split(':')
        mdata = dict([(t[0].strip().lower(),t[1].strip())])
        mdata['source'] = args[0]

        mz = []
        inten = []
        
        for l in f:
            if re.search( '^Name:', l ):
                mdata['prec_mass'] = get_precursor_mass(mdata)
                sorter.push_record( mdata['prec_mass'], {'mz':np.array(mz),
                                                              'inten':np.array(inten)}, mdata )
                num_written += 1
                
                if opts.num_spectra > 0 and num_written == opts.num_spectra:
                    break

                t = l.split(':')
                mdata = dict([(t[0].strip().lower(),t[1].strip())])
                mdata['source'] = args[0]
            
                mz = []
                inten = []
                if num_written % 1000 == 0 and opts.verbose:
                    print( '{}...'.format(num_written), end='', flush=True )
                    
            elif re.search('^Comment:',l):
                mdata.update(parse_comment(l))
            elif re.search('^\\w+:',l):
                t = l.split(':')
                mdata[t[0].strip().lower()] = t[1].strip()
            elif re.search('^\\d+',l):
                t = re.split('\\s+',l)
                mz.append( float(t[0]) )
                inten.append( float(t[1]) )
            
        if not opts.num_spectra > 0 or not num_written == opts.num_spectra:
            mdata['prec_mass'] = get_precursor_mass(mdata)
            sorter.push_record( mdata['prec_mass'], {'mz':np.array(mz),
                                                          'inten':np.array(inten)}, mdata )
            num_written += 1

if opts.verbose: 
    print( 'done.' )
    print( 'Sorting masses...', end='', flush=True )

with sd.CSRSpectralDataWriter(sys.argv[-1],verbose=opts.verbose) as w:
    for i,(mass,arrd,metad) in enumerate(sorter):
        if i % 1000 == 0 and opts.verbose:
            print( '{}...'.format(i), end='', flush=True )
        w.write_record( np.hstack((arrd['mz'].reshape((-1,1)),arrd['inten'].reshape((-1,1)))), metad )

if opts.verbose: 
    print( 'file {} completed.'.format(sys.argv[-1]) )

    
