import pickle, pickletools, re, struct, sys, numpy as np, numpy.random as nr, shelve, tempfile, os, os.path, shutil, glob, h5py, codecs, bisect, JUMPutil.sorting as so

class MzIntenData:
    def __init__( self, mz, inten, idx ):
        self.mzarr = mz
        self.intenarr = inten
        self.idxarr = idx
        
    def mz( self, i=-1 ):
        return self._access_array( self.mzarr, i )

    def inten( self, i=-1 ):
        return self._access_array( self.intenarr, i )

    def _access_array( self, arr, i ):
        if i < 0:
            return arr
        elif i <= self.idx()[-1]:
            lb = bisect.bisect_left( self.idx(), i )
            ub = bisect.bisect_right( self.idx(), i )
            return arr[lb:ub]
        else:
            raise IndexError()

    def idx( self ):
        return self.idxarr

class SpectralData:
    """
    Abstract base class for objects that contain spectral data
    """
    def maxMZ( self ):
        raise TypeError( 'method not implemented' )

    def read_attrs( self, i ):
        raise TypeError( 'method not implemented' )

    def read_spectra( *args ):
        raise TypeError( 'method not implemented' )

    def window_by_precmass( self, minmass, maxmass, tol ):
        raise TypeError( 'method not implemented' )

    def spectra( self ):
        raise TypeError( 'method not implemented' )

    def prec_mass( self, idx ):
        raise TypeError( 'method not implemented' )

    def idx2name( self, *args ):
        raise TypeError( 'method not implemented' )

class LabeledSpectralData(SpectralData):
    def peptides( self ):
        raise TypeError( 'method not implemented' )

    def spectra_by_peptide( self, peptide ):
        raise TypeError( 'method not implemented' )

    def n_spectra( self, peptide ):
        raise TypeError( 'method not implemented' )

class DTASReader(SpectralData):
    def __init__( self, dtas_file ):
        lines = open(dtas_file,'r').read().strip().split('\n')
        i = 0
        self.mz = []
        self.inten = []
        self.data = []
        while i < len(lines):
            name,prec_mass,charge = lines[i].strip().split()
            mz = np.array(lines[i+1].strip().split(),dtype=np.double)
            inten = np.array(lines[i+2].strip().split(),dtype=np.double)
            self.mz.append( mz )
            self.inten.append( inten )
            self.data.append( {'name':name,
                               'prec_mass':float(prec_mass),
                               'charge':int(charge),
                               'mz':mz,
                               'inten':inten} )
            i += 3
        
        self.offset = np.empty(len(self.mz),dtype=np.uint64)
        accum = 0
        for i in range(len(self.offset)):
            self.offset[i] = accum
            accum += self.mz[i].shape[0]
        
        self.mz = np.hstack(self.mz)
        self.inten = np.hstack(self.inten)

    def read_attrs( self, idx ):
        return self.data[idx]

    def prec_mass( self, idx ):
        return self.data[idx]['prec_mass']

    def maxMZ( self ):
        return self.mz.max()

    def minMZ( self ):
        return self.mz.min()

    def idx2name( self, idx ):
        return self.data[idx]['name']

    def spectra( self ):
        return list(range(len(self.data)))

    def read_spectra( self, start_idx, end_idx ):
        start = self.offset[start_idx]
        idxl = np.empty(self.mz.shape[0],dtype=np.uint64)

        k = 0
        for i,j in zip(self.offset[start_idx:end_idx],self.offset[start_idx+1:end_idx+1]):
            idxl[i-start:j-start] = k
            k += 1
    
        idxl[self.offset[end_idx]-start:] = k
        return MzIntenData( self.mz, self.inten, idxl )

def dict2str( kvdata ):
    return str.join('\n',['{}\r{}'.format(k,v) for k,v in kvdata])

def str2dict( s ):
    return dict([l.split('\r') for l in s.strip().split('\n')])

def remap( s ):
    return s.replace('/','@')

def unremap( s ):
    return s.replace('@','/')

class CSRSpectralDataReader(LabeledSpectralData):
    def __init__( self, path, cache_names=True ):
        self.h5file = h5py.File( path, 'r' )
        self.pmass = np.array(self.h5file['data/pmass'])
        self.offset = np.array(self.h5file['data/offset'])
        self.maxMZ_ = self.h5file['data'].attrs['maxMZ']
        self.minMZ_ = self.h5file['data'].attrs['minMZ']
        if cache_names:
            self.peptide_names = codecs.decode(struct.pack('>' + ('B'*self.h5file['meta/names/array'].shape[0]),
                                                           *np.array(self.h5file['meta/names/array'])),
                                               'ascii').split('\n')
        else:
            self.peptide_names = ''
            
    def prec_mass( self, idx ):
        return self.pmass[idx]

    def spectra( self ):
        return np.arange(self.pmass.shape[0])

    def window_by_precmass( self, minmass, maxmass, tol ):
        left_idx = bisect.bisect_left(self.pmass,minmass-tol)
        right_idx = bisect.bisect_right(self.pmass,maxmass+tol)
        return (np.arange(left_idx,right_idx),self.pmass[left_idx:right_idx])

    def maxMZ( self ): return self.maxMZ_

    def minMZ( self ): return self.minMZ_

    def read_attrs( self, idx ):
        return str2dict(self.read_meta_data( idx, 'kvdata', 0 ))

    def idx2name( self, idx ):
        if len(self.peptide_names) > 0:
            return self.peptide_names[idx]
        else:
            return self.read_meta_data( idx, 'name', 1 )

    def read_meta_data( self, idx, key, backup ):
        if idx == self.h5file['meta/{}/offset'.format(key)].shape[0] - 1:
            idx = self.h5file['meta/{}/offset'.format(key)][idx]
            idx_p1 = self.h5file['meta/{}/offset'.format(key)].shape[0]
        else:
            idx,idx_p1 = [int(i) for i in self.h5file['meta/{}/offset'.format(key)][idx:idx+2]]
            
        byts = np.array(self.h5file['meta/{}/array'.format(key)][idx:idx_p1-backup])
        return codecs.decode(struct.pack('>' + ('B'*byts.shape[0]),*byts),'ascii')

    def read_spectra( self, start_idx, end_idx ):
        start = self.offset[start_idx]
        if end_idx < len(self.offset) - 1:
            end = self.offset[end_idx+1]
        else:
            end = self.h5file['data/mz'].shape[0]

        mz = np.empty((int(end-start)))
        inten = np.empty((int(end-start)))
        self.h5file['data/mz'].read_direct(mz,np.s_[start:end])
        self.h5file['data/inten'].read_direct(inten,np.s_[start:end])
        idxl = np.empty(mz.shape[0],dtype=np.uint64)

        k = 0
        for i,j in zip(self.offset[start_idx:end_idx],self.offset[start_idx+1:end_idx+1]):
            idxl[i-start:j-start] = k
            k += 1
    
        idxl[self.offset[end_idx]-start:] = k
        return MzIntenData(mz,inten,idxl)

    def n_spectra( self, peptide ):
        return len(self.h5file['peptides/{}'.format(remap(peptide))])

    def peptides( self ):
        if len(self.peptide_names) > 0:
            return set(self.peptide_names)
        else:
            return set(self.h5file['peptides'])

    def spectra_by_peptide( self, peptide ):
        return np.array(self.h5file['peptides/{}'.format(remap(peptide))])

class CSRSpectralDataWriter(SpectralData):
    def __init__( self, path, verbose=False ):
        self.path = path
        self.h5file = h5py.File( self.path, 'w' )
        self.h5file.create_group( 'peptides' )
        self.h5file.create_group( 'spectra' )

        self.mz = tempfile.TemporaryFile()
        self.inten = tempfile.TemporaryFile()
        self.pmass = tempfile.TemporaryFile()
        self.offset = tempfile.TemporaryFile()
        self.tnames = tempfile.TemporaryFile()
        self.tnames_offset = tempfile.TemporaryFile()
        self.mdata = tempfile.TemporaryFile()
        self.mdata_offset = tempfile.TemporaryFile()
        handle,self.mdataFile = tempfile.mkstemp()
        os.close(handle)
        self.sorter = so.OutOfCoreMerger(1024)
        self.spectra_block_sz = 1024
        
        self.mz_interval = [0,np.inf]

        self.i = 0
        self.end_ptr = 0
        self.tnames_end = 0
        self.mdata_end = 0
        self.verbose = verbose
        
    def write_record( self, peaks, mdata ):
        assert peaks.shape[1] == 2
        assert 'name' in mdata

        self.mz_interval[0] = min(self.mz_interval[0],peaks[:,0].min())
        self.mz_interval[1] = max(self.mz_interval[0],peaks[:,0].max())

        self.mz.write( struct.pack('d'*peaks.shape[0],*(peaks[:,0])) )
        self.inten.write( struct.pack('d'*peaks.shape[0],*(peaks[:,1])) )
        self.pmass.write( struct.pack('d',mdata['prec_mass']) )
        self.offset.write( struct.pack('L',self.end_ptr ) )
        self.end_ptr += peaks.shape[0]

        mdata['lidx'] = self.i
        self.sorter.push_record( np.string_(codecs.encode(mdata['name'],'ascii')), dict(), {'lidx':self.i} )

        self.tnames.write( np.string_(codecs.encode(mdata['name'].replace('\n','') + '\n','ascii')).tostring() )
        self.tnames_offset.write( struct.pack( 'L', self.tnames_end ) )
        self.tnames_end += len(mdata['name']) + 1

        str_mdata = codecs.encode(dict2str(mdata.items()),'ascii')
        self.mdata.write( np.string_(str_mdata).tostring() )
        self.mdata_offset.write( struct.pack( 'L', self.mdata_end ) )
        self.mdata_end += len(str_mdata)

        self.i += 1

    def __enter__( self ):
        return self

    def __exit__( self, type, value, tb ):
        if self.verbose:
            print( 'creating file hierarchy...', end='', flush=True )
            
        self.h5file.create_group('data')
        self.h5file.create_dataset( 'data/mz', dtype=np.double, shape=(self.end_ptr,) )
        self.h5file.create_dataset( 'data/inten', dtype=np.double, shape=(self.end_ptr,) )
        self.h5file.create_dataset( 'data/pmass', dtype=np.double, shape=(self.i,) )
        self.h5file.create_dataset( 'data/offset', dtype=np.uint64, shape=(self.i,) )
        self.h5file.create_group( 'meta' )
        self.h5file.create_group( 'meta/names' )
        self.h5file.create_dataset( 'meta/names/array', dtype=np.uint8, shape=(self.tnames_end,) )
        self.h5file.create_dataset( 'meta/names/offset', dtype=np.uint64, shape=(self.i,) )
        self.h5file.create_group( 'meta/kvdata' )
        self.h5file.create_dataset( 'meta/kvdata/array', dtype=np.uint8, shape=(self.mdata_end,) )
        self.h5file.create_dataset( 'meta/kvdata/offset', dtype=np.uint64, shape=(self.i,) )
        
        self.h5file['data'].attrs['minMZ'] = self.mz_interval[0]
        self.h5file['data'].attrs['maxMZ'] = self.mz_interval[1]
        if self.verbose:
            print('done.')
            print('mapping spectra to peptide groups...', end='', flush=True )
        
        cur_pep = ''
        spec_idxs = []
        for i,(name,empt1,md) in enumerate(self.sorter):
            if self.verbose and i % 1000 == 0:
                print( str(i) + '...', end='', flush=True )
            str_name = codecs.decode(name,'ascii')
            if str_name != cur_pep:
                if len(cur_pep) > 0:
                    self.h5file.create_dataset( 'peptides/{}'.format(remap(cur_pep)), data=spec_idxs )
                    
                spec_idxs = [md['lidx']]
                cur_pep = str_name
            else:
                spec_idxs.append( md['lidx'] )

        self.h5file.create_dataset( 'peptides/{}'.format(remap(cur_pep)), data=spec_idxs )
        if self.verbose:
            print( 'done.' )
            print( 'writing mz/intensity data...', end='', flush=True )
                
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
        #    
        self.tnames.flush()
        self.tnames.seek(0)
        i = 0
        while True:
            tnames = self.tnames.read(blksz*struct.calcsize('B'))
            self.h5file['meta/names/array'][i*blksz:(i+1)*blksz] = struct.unpack('>' + ('B'*(len(tnames)//struct.calcsize('B'))),
                                                                                 tnames)
            i += 1
            if len(tnames)/struct.calcsize('B') < blksz:
                break

        self.mdata.flush()
        self.mdata.seek(0)
        i = 0
        while True:
            mdata = self.mdata.read(blksz*struct.calcsize('B'))
            self.h5file['meta/kvdata/array'][i*blksz:(i+1)*blksz] = struct.unpack('>' + ('B'*(len(mdata)//struct.calcsize('B'))),
                                                                                  mdata)
            i += 1
            if len(mdata)/struct.calcsize('B') < blksz:
                break

        self.pmass.flush()
        self.offset.flush()
        self.tnames_offset.flush()
        self.mdata_offset.flush()
        self.pmass.seek(0)
        self.offset.seek(0)
        self.tnames_offset.seek(0)
        self.mdata_offset.seek(0)
        i = 0
        while True:
            pmass = self.pmass.read(blksz*struct.calcsize('d'))
            offset = self.offset.read(blksz*struct.calcsize('L'))
            tnames_offset = self.tnames_offset.read(blksz*struct.calcsize('L'))
            mdata_offset = self.mdata_offset.read(blksz*struct.calcsize('L'))
            self.h5file['data/pmass'][i*blksz:(i+1)*blksz] = struct.unpack('d'*(len(pmass)//struct.calcsize('d')),
                                                                           pmass)
            self.h5file['data/offset'][i*blksz:(i+1)*blksz] = struct.unpack('L'*(len(offset)//struct.calcsize('L')),
                                                                            offset)
            self.h5file['meta/names/offset'][i*blksz:(i+1)*blksz] = struct.unpack('L'*(len(tnames_offset)//struct.calcsize('L')),
                                                                                  tnames_offset)
            self.h5file['meta/kvdata/offset'][i*blksz:(i+1)*blksz] = struct.unpack('L'*(len(mdata_offset)//struct.calcsize('L')),
                                                                                   mdata_offset)

            i += 1
            if len(pmass)/struct.calcsize('d') < blksz:
                break

        if self.verbose:
            print( 'done.' )
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

