import JUMPsl.spectral_data as sd, heapq, numpy as np, h5py, tempfile, os.path, re, bz2

class OutOfCoreMerger:
    def __init__( self, bufsz=1048576 ):
        self.tdir = tempfile.TemporaryDirectory()
        self.h5file = h5py.File(os.path.join(self.tdir.name,'_cache.h5'),'w')
        self.data_buf = []
        self.bufsz = bufsz

    def push_record( self, idxval, array_data, kvdata ):
        if len(self.data_buf) >= self.bufsz:
            self.flush()
        self.data_buf.append( (idxval,array_data,kvdata) )
    
    def flush( self ):
        g = self.h5file.create_group(str(len(self.h5file)))
        self.data_buf.sort(key = lambda t: t[0])
        g.create_dataset( 'idx', data=[t[0] for t in self.data_buf] )
        
        keys = self.data_buf[0][1].keys()
        for k in keys:
            gs = g.create_group(k)
            gs.create_dataset( 'data', data=np.hstack([self.data_buf[i][1][k] for i in range(len(self.data_buf))]) )
            accum = 0
            offsets = []
            for i in range(len(self.data_buf)):
                offsets.append( accum )
                accum += self.data_buf[i][1][k].shape[0]

            offsets.append( accum )
            gs.create_dataset( 'offsets', data=offsets )

        for i in range(len(self.data_buf)):
            g.attrs[str(i)] = repr(self.data_buf[i][2])

        self.data_buf = []

    def idx_value( self, i, j ):
        return self.h5file[str(i)]['idx'][j]

    def nitems( self, i ):
        return self.h5file[str(i)]['idx'].shape[0]

    def read( self, i, j ):
        arr_data = {}
        for k in self.h5file[str(i)]:
            if k != 'idx':
                start = self.h5file[str(i)][k]['offsets'][j]
                end = self.h5file[str(i)][k]['offsets'][j+1]
                arr_data[k] = np.array(self.h5file[str(i)][k]['data'][start:end])

        kvdata = eval(self.h5file[str(i)].attrs[str(j)])
        return (arr_data,kvdata)

    def __iter__( self ):
        self.flush()

        for i in range(len(self.h5file)):
            heapq.heappush(self.data_buf,(self.idx_value(i,0),i,0))

        return self

    def __next__( self ):
        if len(self.data_buf) == 0:
            raise StopIteration

        t = heapq.heappop(self.data_buf)
        j = t[-1]
        if j + 1 < self.nitems(t[1]):
            heapq.heappush(self.data_buf,(self.idx_value(t[1],j+1),t[1],j+1))

        arr_data,kvdata = self.read(t[1],j)
        return (t[0],arr_data,kvdata)

def parse_comment( comment_line ):
    return dict([(m.group(1).lower(),m.group(2)) for m in
                 re.finditer( '(\\S+)=(\\S+)', comment_line )])

def multiopen( fname ):
    if '.bz2' == fname[-4:]:
        return bz2.open(fname,'rt')
    else:
        return open(fname)

if __name__ == '__main__':
    import argparse, sys
    p = argparse.ArgumentParser(description='import a library from MSP or SPTXT (SpectraST\'s text format)')
    p.add_argument( '-v', '--verbose', dest='verbose', action='store_true' )
    p.add_argument( '-n', '--num-spectra', dest='num_spectra', type=int, default=0, help='only process the first N spectra', metavar='N' )
    p.add_argument( '-p', '--precursor-mass-key', dest='precursor_mass_key',
                    help='key for precursor mass in spectra metadata.  Not case sensitive',
                    default='precursormz' )
    p.add_argument('import_files', help='MSP or SPTXT files to be imported', nargs='+')
    p.add_argument('output', help='name of output file')

    opts = p.parse_args()

    def get_precursor_mass( mdata ):
        if opts.precursor_mass_key.lower() in mdata:
            return float(mdata[opts.precursor_mass_key.lower()])
        else:
            return float(mdata['parent'])

    if os.path.exists(sys.argv[-1]):
        raise Exception('will not overwrite file at ' + sys.argv[-1])

    sorter = OutOfCoreMerger(1024)
    num_written = 0
    for txtfile in opts.import_files:
        if opts.verbose:
            if num_written > 0: print('')
            print( 'processing file {}...'.format(txtfile), end='', flush=True )
        f = multiopen(txtfile)
        for l in f:
            if re.search( '^Name:', l ):
                break

        t = l.split(':')
        mdata = dict([(t[0].strip().lower(),t[1].strip())])
        mdata['source'] = txtfile

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
                mdata['source'] = txtfile
            
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

with sd.CSRSpectralDataWriter(opts.output) as w:
    for i,(mass,arrd,metad) in enumerate(sorter):
        if i % 1000 == 0 and opts.verbose:
            print( '{}...'.format(i), end='', flush=True )
        w.write_record( np.hstack((arrd['mz'].reshape((-1,1)),arrd['inten'].reshape((-1,1)))), metad )

if opts.verbose: 
    print( 'done.' )

    
