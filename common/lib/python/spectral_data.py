import pickle, pickletools, re, struct, sys, numpy as np, numpy.random as nr, shelve, tempfile, os, os.path, shutil, glob, h5py, codecs, bisect

class SpectralData:
    """
    Abstract base class for objects that contain spectral data
    """
    def maxMZ( self ):
        raise TypeError( 'method not implemented' )

    def read_attrs( self, i ):
        raise TypeError( 'method not implemented' )

    def read_peptide( *args ):
        raise TypeError( 'method not implemented' )

    def window_by_precmass( self, minmass, maxmass, tol ):
        raise TypeError( 'method not implemented' )

    def spectra( self ):
        raise TypeError( 'method not implemented' )

class DTASReader(SpectralData):
    def __init__( self, dtas_file ):
        lines = open(dtas_file,'r').read().strip().split('\n')
        i = 0
        self.data = []
        while i < len(lines):
            name,prec_mass,charge = lines[i].strip().split()
            mz = np.array(lines[i+1].strip().split(),dtype=np.double)
            inten = np.array(lines[i+2].strip().split(),dtype=np.double)
            self.data.append( {'name':name,
                               'prec_mass':float(prec_mass),
                               'charge':int(charge),
                               'mz':mz,
                               'inten':inten} )
            i += 3

    def read_attrs( self, idx ):
        return self.data[idx]

    def maxMZ( self ):
        return max([d['mz'].max() for d in self.data])

    def minMZ( self ):
        return min([d['mz'].max() for d in self.data])

    def spectra( self ):
        return list(range(len(self.data)))

    def read_peptide( self, idx ):
        return np.hstack((self.data[idx]['mz'].reshape((-1,1)),
                          self.data[idx]['inten'].reshape((-1,1))))

def remap( s ):
    return s.replace('/','@')

def unremap( s ):
    return s.replace('@','/')

class CSRSpectralDataReader(SpectralData):
    def __init__( self, path ):
        self.h5file = h5py.File( path, 'r' )
        self.pmass = np.array(self.h5file['data/pmass'])
        self.offset = np.array(self.h5file['data/offset'])
        self.names = codecs.decode(np.array(self.h5file['data/names']).tostring(),
                                   'ascii').split('\n')
        self.maxMZ_ = self.h5file['data'].attrs['maxMZ']
        self.minMZ_ = self.h5file['data'].attrs['minMZ']

    def spectra( self ):
        return np.arange(self.pmass.shape[0])

    def window_by_precmass( self, minmass, maxmass, tol ):
        left_idx = bisect.bisect_left(self.pmass,minmass-tol)
        right_idx = bisect.bisect_right(self.pmass,maxmass+tol)
        return (np.arange(left_idx,right_idx),self.pmass[left_idx:right_idx])

    def maxMZ( self ): return self.maxMZ_

    def minMZ( self ): return self.minMZ_

    def read_attrs( self, i ):
        return dict(self.h5file[str(i)].attrs)

    def block_read_peptides( self, start_idx, end_idx ):
        start = self.offset[start_idx]
        if end_idx < len(self.offset) - 1:
            end = self.offset[end_idx+1]
        else:
            end = self.h5file['data/mz'].shape[0]

        mz_inten = np.empty((int(end-start),2))
        self.h5file['data/mz'].read_direct(mz_inten,np.s_[start:end],np.s_[:,0])
        self.h5file['data/inten'].read_direct(mz_inten,np.s_[start:end],np.s_[:,1])
        idxl = np.empty(mz_inten.shape[0],dtype=np.uint64)

        k = 0
        for i,j in zip(self.offset[start_idx:end_idx-1],self.offset[start_idx+1:end_idx]):
            idxl[i-start:j-start] = k
            k += 1
    
        idxl[self.offset[end_idx]-start:] = k
        return (mz_inten,idxl)

    def read_peptide( self, i ):
        start = self.offset[i]
        if i < len(self.offset) - 1:
            end = self.offset[i+1]
        else:
            end = self.h5file['data/mz'].shape[0]

        return np.hstack((np.array(self.h5file['data/mz'][start:end]).reshape((-1,1)),
                          np.array(self.h5file['data/inten'][start:end]).reshape((-1,1))))

class CSRSpectralDataWriter:
    def __init__( self, path ):
        self.path = path
        self.h5file = h5py.File( self.path, 'w' )

        self.mz = tempfile.TemporaryFile()
        self.inten = tempfile.TemporaryFile()
        self.pmass = tempfile.TemporaryFile()
        self.offset = tempfile.TemporaryFile()
        self.names = []

        self.i = 0
        self.end_ptr = 0
        
    def write_record( self, peaks, mdata ):
        assert peaks.shape[1] == 2
        self.mz.write( struct.pack('d'*peaks.shape[0],*(peaks[:,0])) )
        self.inten.write( struct.pack('d'*peaks.shape[0],*(peaks[:,1])) )
        self.pmass.write( struct.pack('d',mdata['prec_mass']) )
        self.offset.write( struct.pack('L',self.end_ptr ) )
        self.names.append( mdata['name'] )
        self.end_ptr += peaks.shape[0]

        h5data = self.h5file.create_group( 'spectra/{}'.format(self.i) )
        self.i += 1
        for k,v in mdata.items():
            h5data[k] = v

    def __enter__( self ):
        return self

    def __exit__( self, type, value, tb ):
        self.h5file.create_group('data')
        self.h5file.create_dataset( 'data/mz', dtype=np.double, shape=(self.end_ptr,) )
        self.h5file.create_dataset( 'data/inten', dtype=np.double, shape=(self.end_ptr,) )
        self.h5file.create_dataset( 'data/pmass', dtype=np.double, shape=(self.i,) )
        self.h5file.create_dataset( 'data/offset', dtype=np.uint64, shape=(self.i,) )
        
        self.mz.flush()
        self.inten.flush()
        self.mz.seek(0)
        self.inten.seek(0)
        blksz = 16777216
        i = 0
        while True:
            mz = self.mz.read(blksz*struct.calcsize('d'))
            inten = self.inten.read(blksz*struct.calcsize('d'))
            self.h5file['data/mz'][i*blksz:(i+1)*blksz] = struct.unpack('d'*(len(mz)//struct.calcsize('d')),
                                                                        mz)
            self.h5file['data/inten'][i*blksz:(i+1)*blksz] = struct.unpack('d'*(len(inten)//struct.calcsize('d')),
                                                                           inten)
            i += 1
            if len(mz)/struct.calcsize('d') < blksz:
                break
            
        self.pmass.flush()
        self.offset.flush()
        self.pmass.seek(0)
        self.offset.seek(0)
        i = 0
        while True:
            pmass = self.pmass.read(blksz*struct.calcsize('d'))
            offset = self.offset.read(blksz*struct.calcsize('L'))
            self.h5file['data/pmass'][i*blksz:(i+1)*blksz] = struct.unpack('d'*(len(pmass)//struct.calcsize('d')),
                                                                           pmass)
            self.h5file['data/offset'][i*blksz:(i+1)*blksz] = struct.unpack('L'*(len(offset)//struct.calcsize('L')),
                                                                            offset)
            i += 1
            if len(pmass)/struct.calcsize('d') < blksz:
                break

        self.h5file['data'].create_dataset( 'names', data=np.string_(str.join('\n',self.names)) )

        self.h5file.flush()
        self.h5file.close()

class H5SpectralDataWriter:
    def __init__( self, path ):
        self.path = path
        self.tempPath = tempfile.TemporaryDirectory()
        self.h5file = h5py.File(os.path.join(self.tempPath.name,
                                             os.path.basename(self.path)),'w',
                                libver='latest')
        self.h5file.create_group( 'index' )
        self._num_written = 0

    def __len__( self ):
        return self.num_written()

    def num_written( self ):
        return self._num_written

    def write_record( self, peaks, mdata ):
        key = remap(mdata['Name'])
        if key not in self.h5file:
            self.h5file.create_group( key )

        idx = str(len(self.h5file[key]))
        self.h5file[key].create_dataset( idx, peaks.shape, peaks.dtype )
        self.h5file[key][idx][...] = peaks
        for k,v in mdata.items():
            self.h5file[key][idx].attrs[k] = v

        self.h5file[key][idx].attrs['maxMZ'] = str(peaks[:,0].max())
        jdx = str(len(self.h5file['index']))
        self.h5file[key][idx].attrs['index'] = jdx
        self.h5file['index'][jdx] = h5py.SoftLink('/{}/{}'.format(key,idx))
        self._num_written += 1        

    def __enter__( self ):
        return self
        
    def __exit__( self, type, value, tb ):
        self.h5file.attrs['maxMZ'] = max([float(self.h5file['index'][str(i)].attrs['maxMZ'])
                                          for i in range(len(self.h5file['index']))])
        peptides = [k for k in self.h5file if k != 'index']
        masses = [float(self.h5file['index'][str(i)].attrs['PrecursorMZ'])
                  for i in range(len(self.h5file['index']))]
        
        maxNameLen = max([len(k) for k in peptides])
        self.h5file.create_dataset( 'lenMapKeys', (len(peptides),), dtype='S'+str(maxNameLen) )
        self.h5file.create_dataset( 'lenMapValues', (len(peptides),), dtype=np.uint )
        medianMass = np.empty(len(peptides))
        
        for i,p in enumerate(peptides):
            precMasses = np.array([float(self.h5file[p][str(i)].attrs['PrecursorMZ'])
                                   for i in range(len(self.h5file[p]))]) 
            medianMass[i] = np.median(precMasses)
        
        for i in range(len(peptides)):
            self.h5file['lenMapKeys'][i] = np.string_(peptides[i])
            self.h5file['lenMapValues'][i] = len(self.h5file[peptides[i]])
             
        self.h5file.create_dataset( 'MedPrecMassMap', (len(peptides),), dtype=np.double )
        self.h5file['MedPrecMassMap'][...] = medianMass[...]
        pm = dict([(v,i) for i,v in enumerate(sorted(set(peptides)))])
        self.h5file.create_dataset( 'indexMap', data=[pm[self.h5file['index'][str(i)].attrs['Name']]
                                                      for i in range(self.h5file['index'].shape[0])] )
        self.h5file.close()
        for f in glob.glob(os.path.join(self.tempPath.name,'*')):
            bf = os.path.basename(f)
            if os.path.exists(os.path.join(os.path.dirname(self.path),bf)):
                os.unlink(os.path.join(os.path.dirname(self.path),bf))
            
            shutil.move(f,os.path.join(os.path.dirname(self.path),bf))

class H5SpectralDataReader:
    def __init__( self, filePath ):
        self.h5file = h5py.File( filePath, 'r' )
        self.filePath = filePath
        self.psmMap = dict(zip([unremap(k.decode('ascii')) for k in np.array(self.h5file['lenMapKeys'])],
                               np.array(self.h5file['lenMapValues'])))
        self.massMap = dict(zip([unremap(k.decode('ascii')) for k in np.array(self.h5file['lenMapKeys'])],
                               np.array(self.h5file['MedPrecMassMap'])))

        self._len = len(self.h5file['index'])

    def spectra( self ):
        return np.array(self.h5file['index'])

    def typeMap( self ):
        return np.array(self.h5file['indexMap'])

    def peptides( self ):
        return set(self.psmMap.keys())

    def window_by_precmass( self, minmass, maxmass, tol ):
        idxs = []
        masses = []
        for p in self.peptides():
            mmass = self.medianPrecMass( p )
            if (mmass >= minmass - tol and mmass <= maxmass + tol):
                n = self.nPSMs(p)
                idxs += [(p,i) for i in range(n)]
                masses += [mmass]*n
        return (idxs,masses)

    def medianPrecMass( self, peptide ):
        return self.massMap[peptide]

    def precMasses( self, peptide ):
        return self.h5file[remap(peptide)]['precMasses']

    def __len__( self ):
        return self._len

    def maxMZ( self ):
        return float(self.h5file.attrs['maxMZ'])

    def nPSMs( self, peptide ):
        return self.psmMap[peptide]

    def get_index( self, peptide, pidx ):
        return int(self.read_attrs(peptide,pidx)['index'])

    def get_name( self, idx ):
        return self.read_attrs(idx)['Name']

    def read_peptide( *args ):
        self = args[0]
        if len(args) == 2 and isinstance(args[1],tuple):
            return np.array(self.h5file[remap(args[1][0])+'/'+str(args[1][1])])
        elif len(args) == 2:
            return np.array(self.h5file['index'][str(args[1])])
        elif len(args) == 3:
            return np.array(self.h5file[remap(args[1])][str(args[2])])
        else:
            raise KeyError('wrong number of args to read_peptide')

    def read_attrs( *args ):
        self = args[0]
        if len(args) == 2:
            return self.h5file['index'][str(args[1])].attrs
        elif len(args) == 3:
            return self.h5file[remap(args[1])][str(args[2])].attrs
        else:
            raise KeyError('wrong number of args to read_peptide')
        
class H5ConsensusReader(H5SpectralDataReader):
    def __init__( *args ):
        super(H5ConsensusReader,args[0]).__init__( *args[1:] )
        self = args[0]
        self.indexMap = list(self.peptides())
        self.rIndexMap = dict([(p,i) for i,p in enumerate(self.indexMap)])

    def nPSMs( self, peptide ):
        return 1

    def get_name( self, idx ):
        return self.indexMap[idx]

    def get_index( self, peptide, pidx ):
        assert pidx == 0
        return self.rIndexMap[peptide]

    def read_peptide( *args ):
        self = args[0]
        if len(args) == 2:
            pname = self.indexMap[args[1]]
        elif len(args) == 3:
            pname = args[1]
            assert args[2] == 0
        else:
            raise KeyError('wrong number of args to read_peptide')
        return np.vstack([super(H5ConsensusReader,self).read_peptide( pname, i )
                          for i in range(super(H5ConsensusReader,self).nPSMs(pname))])
        
    def read_attrs( *args ):
        self = args[0]
        if len(args) == 2:
            pname = self.indexMap[args[1]]
        elif len(args) == 3:
            pname = args[1]
            assert args[2] == 0
        else:
            raise KeyError('wrong number of args to read_peptide')
        return [super(H5ConsensusReader,self).read_attrs( pname, i )
                for i in range(super(H5ConsensusReader,self).nPSMs(pname))]

class SpectralDataReader:
    def __init__( *args, **kwargs ):
        self = args[0]
        if len(args) > 1:
            self.initialize( *(args[1:]), **kwargs )
        
    def initialize( self, path, keep_all_mdata=False ):
        self.path = path
        self.binfile = open(path + '.bin','rb')
        self.pep_map = [{'Name':t[0],'MW':t[1],'length':t[2],'offset':t[3]}
                        for t in zip(open(path+'.nam').read().strip().split('\n'),
                                     np.load(path + '.mw.npy'),
                                     np.load(path + '.len.npy'),
                                     np.load(path + '.off.npy'))]
        
        self.maxMW = max([float(d['MW']) for d in self.pep_map])
                
    def clone( self ):
        other = SpectralDataReader()
        other.path = self.path
        other.binfile = open(self.path + '.bin','rb')
        if isinstance(self.pep_map,dict):
            other.pep_map = dict([(k,dict(v)) for k,v in self.pep_map.items()])
        else:
            other.pep_map = [dict(d) for d in self.pep_map]
        other.mdata = shelve.open(self.path + '.db')
        other.maxMW = self.maxMW
        return other

    def __len__( self ):
        return len(self.pep_map)

    def shrink( self, pep_set ):
        self.pep_map = dict([(i,self.pep_map[i]) for i in pep_set])

    def myPeptides( self ):
        return [d['Name'] for d in self.pep_map]

    def read_record( self, i ):
        return read_peaks( self.binfile, self.pep_map[i]['offset'], self.pep_map[i]['length'] )

class DownsampledSpectralDataReader(SpectralDataReader):
    def __init__( self, path, pct, seed=0, pepSet=None ):
        super(DownsampledSpectralDataReader,self).__init__( path )
        nr.seed(seed)
        self.choice = nr.choice( len(self), int(len(self)*pct), replace=False )
        if None != pepSet:
            self.choice = [i for i in self.choice if self.pep_map[i]['Name'] in pepSet]
        
        self.pep_map = [self.pep_map[i] for i in self.choice]

class SpectralDataCollection:
    def __init__( self, paths ):
        self.readers = dict([(p,SpectralDataReader(p)) for p in paths])
        self.nindex = {}
        for p,r in self.readers.items():
            for i,item in enumerate(r.pep_map):
                if item['Name'] in self.nindex:
                    self.nindex[item['Name']].append( (p,i) )
                else:
                    self.nindex[item['Name']] = [(p,i)]
        
        self.names = list(self.nindex.keys())
        try:
            self.maxMW = max([r.maxMW for r in self.readers.values()])
        except ValueError:
            self.maxMW = None

    def clone( self ):
        other = SpectralDataCollection([])
        other.readers = dict([(k,v.clone()) for k,v in self.readers.items()])
        other.nindex = dict([(k,list(v)) for k,v in self.nindex.items()])
        other.names = list(self.names)
        other.maxMW = self.maxMW
        return other

    def shrink( self, pep_set ):
        for p,r in self.readers.items():
            pep_by_reader = set([t[1] for t in pep_set if t[0] == p])
            r.shrink( pep_by_reader )

        self.nindex = {}
        for p,r in self.readers.items():
            for i,item in r.pep_map.items():
                if item['Name'] in self.nindex:
                    self.nindex[item['Name']].append( (p,i) )
                else:
                    self.nindex[item['Name']] = [(p,i)]
        
        self.names = list(self.nindex.keys())
        self.maxMW = max([r.maxMW for r in self.readers.values()])
        
    def items( self ):
        items = []
        for reader,list_of_peptides in self.nindex.items():
            items += [(reader,(p,i)) for p,i in list_of_peptides]

        return sorted(items)

    def draw_sample( self, n, name, complement=False ):
        if complement:
            names = [self.names[i] for i in nr.choice(np.arange(len(self.names)),n+1,replace=False)
                     if self.names[i] != name]
            readerhits = [self.nindex[nom][nr.randint(0,len(self.nindex[nom]))] for nom in names[:n]]
        else:
            readerhits = [self.nindex[name][i] for i in
                          nr.choice(np.arange(len(self.nindex[name])),n,replace=False)]
        
        return [(name,t) for t in readerhits]

def write_peaks( outfile, peaks ):
    offset = outfile.tell()
    length = outfile.write(struct.pack( 'd'*len(peaks), *([t[0] for t in peaks])))
    length += outfile.write(struct.pack( 'd'*len(peaks), *([t[1] for t in peaks])))
    return (offset,length)
    
def read_peaks( infile, offset, length ):
    infile.seek(offset)
    bindata = infile.read(length)
    return zip(struct.unpack('d'*((length//struct.calcsize('d'))//2),bindata[:length//2]),
               struct.unpack('d'*((length//struct.calcsize('d'))//2),bindata[length//2:]))

