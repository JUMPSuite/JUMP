import pickle, pickletools, re, struct, sys, numpy as np, numpy.random as nr, shelve, tempfile, os, os.path, shutil, glob, h5py, codecs, bisect, JUMPutil.sorting as so

class MzIntenData:
    """
    Abtract base class for m/z,intensity data to be packaged and returned from :class:`SpectralData` instances
    """
    def __init__( self, mz, inten, idx ):
        self.mzarr = mz
        self.intenarr = inten
        self.idxarr = idx
        
    def mz( self, i=-1 ):
        """
        Get m/z data at index i, or all data if no i is provided
        """
        return self._access_array( self.mzarr, i )

    def inten( self, i=-1 ):
        """
        Get intensity data at index i, or all data if no i is provided
        """
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
        """
        Get indices for all data in the array
        """
        return self.idxarr

class SpectralData:
    """
    Abstract base class for objects that contain spectral data
    """
    def maxMZ( self ):
        """
        Get max m/z for all spectra in self
        """
        raise TypeError( 'method not implemented' )

    def read_attrs( self, i ):
        """
        Return a dictionary containing metadata for spectum at index i
        """
        raise TypeError( 'method not implemented' )

    def read_spectra( *args ):
        """
        Return a :class:`MzIntenData` instance for spectra indexed by args
        """
        raise TypeError( 'method not implemented' )

    def window_by_precmass( self, minmass, maxmass, tol ):
        """
        Get the indices in self that have precusor mass in the interval [minmass-tol,maxmass+tol]
        """
        raise TypeError( 'method not implemented' )

    def spectra( self ):
        """
        Return all indices of spectra in self
        """
        raise TypeError( 'method not implemented' )

    def prec_mass( self, idx ):
        """
        Get the precursor mass of spectrum at index i
        """
        raise TypeError( 'method not implemented' )

    def idx2name( self, *args ):
        """
        Get the spectrum name for indices in args
        """
        raise TypeError( 'method not implemented' )

class LabeledSpectralData(SpectralData):
    """
    Abstract class for SpectralData that has been identified and labeled
    """
    def peptides( self ):
        """
        Return all the peptides in self
        """
        raise TypeError( 'method not implemented' )

    def spectra_by_peptide( self, peptide ):
        """
        Get the indices of all spectra in self for peptide.
        """
        raise TypeError( 'method not implemented' )

    def n_spectra( self, peptide ):
        """
        Get the number of spectra in self.
        """
        raise TypeError( 'method not implemented' )

class DTASReader(SpectralData):
    """
    :class:`SpectralData` front-end for DTAS data
    """
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
        if end_idx < len(self.offset) - 1:
            end = self.offset[end_idx]
        else:
            end = self.mz.shape[0]

        start = self.offset[start_idx]
        mz = self.mz[start:end]
        inten = self.inten[start:end]
        idxl = np.empty(mz.shape[0],dtype=np.uint64)

        k = 0
        for i,j in zip(self.offset[start_idx:end_idx],self.offset[start_idx+1:end_idx+1]):
            idxl[i-start:j-start] = k
            k += 1
    
        if end_idx < len(self.offset) - 1:
            idxl[self.offset[end_idx]-start:] = k
        else:
            idxl[self.offset[end_idx-1]-start:] = k

        return MzIntenData( mz, inten, idxl )

def dict2str( kvdata ):
    return str.join('\n',['{}\r{}'.format(k,v) for k,v in kvdata])

def str2dict( s ):
    return dict([l.split('\r') for l in s.strip().split('\n')])

def remap( s ):
    return s.replace('/','@')

def unremap( s ):
    return s.replace('@','/')

class CSRSpectralDataReader(LabeledSpectralData):
    """
    :class:`LabeledSpectralData` instance for reading JUMPs own HDF5 file format
    """
    def __init__( self, path, cache_names=True ):
        self.h5file = h5py.File( path, 'r' )
        self.pmass = np.array(self.h5file['data/pmass'])
        self.offset = np.array(self.h5file['data/offset'])
        self.maxMZ_ = self.h5file['data'].attrs['maxMZ']
        self.minMZ_ = self.h5file['data'].attrs['minMZ']
        if cache_names:
            self.peptide_names = codecs.decode(struct.pack('>' + ('B'*self.h5file['meta/names/array'].shape[0]),
                                                           *np.array(self.h5file['meta/names/array'])),
                                               'ascii','replace').split('\n')
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
        return codecs.decode(struct.pack('>' + ('B'*byts.shape[0]),*byts),'ascii','replace')

    def read_spectra( self, start_idx, end_idx ):
        start = self.offset[start_idx]
        if end_idx < len(self.offset) - 1:
            end = self.offset[end_idx]
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
    
        if end_idx < len(self.offset):
            idxl[self.offset[end_idx]-start:] = k
        else:
            idxl[self.offset[end_idx-1]-start:] = k

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
    """
    Class for writing JUMP's own HDF5 library format
    """
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
        
        self.mz_interval = [np.inf,0]

        self.i = 0
        self.end_ptr = 0
        self.tnames_end = 0
        self.mdata_end = 0
        self.verbose = verbose
        
    def write_record( self, peaks, mdata ):
        assert peaks.shape[1] == 2
        assert 'name' in mdata

        self.mz_interval[0] = min(self.mz_interval[0],peaks[:,0].min())
        self.mz_interval[1] = max(self.mz_interval[1],peaks[:,0].max())

        self.mz.write( struct.pack('d'*peaks.shape[0],*(peaks[:,0])) )
        self.inten.write( struct.pack('d'*peaks.shape[0],*(peaks[:,1])) )
        self.pmass.write( struct.pack('d',mdata['prec_mass']) )
        self.offset.write( struct.pack('L',self.end_ptr ) )
        self.end_ptr += peaks.shape[0]

        mdata['lidx'] = self.i
        self.sorter.push_record( np.string_(codecs.encode(mdata['name'],'ascii','replace')), dict(), {'lidx':self.i} )

        self.tnames.write( np.string_(codecs.encode(mdata['name'].replace('\n','') + '\n','ascii','replace')).tostring() )
        self.tnames_offset.write( struct.pack( 'L', self.tnames_end ) )
        self.tnames_end += len(mdata['name']) + 1

        str_mdata = codecs.encode(dict2str(mdata.items()),'ascii','replace')
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
            str_name = codecs.decode(name,'ascii','replace')
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

class CSRFilterDataWriter:
    def __init__( self, h5fileName ):
        self.h5file = h5py.File( h5fileName, 'a' )
        if 'data/filters/masses' in self.h5file:
            self.mass_array = np.array(self.h5file['data/filters/masses'])

    def initialize_masses( self, masses, ninterp_pts ):
        if 'data' not in self.h5file: 
            self.h5file.create_group( 'data' )

        self.h5file.create_group( 'data/filters' )
        self.h5file['data/filters'].create_dataset( 'masses', data=masses )
        self.h5file['data/filters'].create_dataset( 'filters', dtype=np.double, shape=(len(masses),ninterp_pts) )

    def write_filter( self, mass, filter ):
        idx = bisect.bisect( self.mass_array, mass )
        assert self.mass_array[idx-1] == mass
        self.h5file['data/filters/filters'][idx-1,:] = filter.reshape((-1,))

    def close( self ):
        self.h5file.close()

class CSRFilterDataReader:
    def read_filter( self, mass, filter ):
        idx = bisect.bisect( self.mass_array, mass )
        if mass - self.mass_array[idx-1] < self.mass_array[idx] - mass:
            return np.array(self.h5file['data/filters/filters'][idx-1,:])
        else:
            return np.array(self.h5file['data/filters/filters'][idx,:])

