import JUMPsl.binning as bi, JUMPsl.spectral_data as sd, scipy.sparse as ss, numpy as np, bisect
def build_binned_mat( binner, collection, spectra, col_or_row='row' ):
    if np.all(spectra == np.arange(min(spectra),max(spectra)+1)):
        mzint = collection.read_spectra( min(spectra), max(spectra)+1 )
        c,data = binner.bin( mzint.mz(), mzint.inten() )
    
        Al = ss.coo_matrix( (data,(mzint.idx(),c)), shape=(len(spectra),binner.nbins) )
    else:
        mzint_data = [collection.read_spectra( s, s+1 ) for s in spectra]
        c,data = binner.bin( np.hstack([mzd.mz() for mzd in mzint_data]), 
                             np.hstack([mzd.inten() for mzd in mzint_data]) )

        Al = ss.coo_matrix( (data,(np.hstack([np.ones(mzd.mz().shape)*i for i,mzd in enumerate(mzint_data)]),
                                   c)), shape=(len(spectra),binner.nbins) )

    if 'row' == col_or_row:
        Al = Al.tocsr()
    else:
        Al = Al.tocsc()
    D = ss.dia_matrix( (np.power(Al.multiply(Al).sum(1).T + 
                                 (0 == Al.multiply(Al).sum(1).T),-.5),0), shape=(Al.shape[0],Al.shape[0]) )
        
    return (D,Al,[collection.idx2name(i) for i in spectra])

def get_minmax_median_precmass( search_spectra ):
    prec_masses = [search_spectra.prec_mass(p) for p in search_spectra.spectra()]
    return (min(prec_masses),max(prec_masses))

def precmass_window( minmass, maxmass, tol, spectra_collection ):
    wanted_spectra = []
    for p in spectra_collection.spectra():
        med_mass = spectra_collection.prec_mass(p)
        if med_mass >= minmass - tol and med_mass <= maxmass + tol:
            wanted_spectra.append( (p,med_mass) )

    return ([t[0] for t in wanted_spectra],[t[1] for t in wanted_spectra])

class BinnedSearch:
    """
    Basic binned search object
    """
    def __init__( self, library_collection, max_mz, binsz, tol, n,
                  blksz=2, search_subset=None, pmass_window=None ):
        self.binner = bi.Binner( binsz, max(library_collection.maxMZ(),
                                            max_mz), library_collection )

        if None != pmass_window:
            self.minmass = pmass_window[0]
            self.maxmass = pmass_window[1]
            self.spectra,self.med_masses = library_collection.window_by_precmass( self.minmass, 
                                                                                  self.maxmass, tol )
        else:
            self.spectra = library_collection.spectra()
            self.med_masses = np.array([library_collection.prec_mass(s) for s in self.spectra])
            self.minmass = self.med_masses.min()
            self.maxmass = self.med_masses.max()

        self.n = n
        self.tol = tol
        if len(self.spectra) == 0:
            self.empty = True
        else:
            self.empty = False
            self.med_masses = np.array(self.med_masses)
            self.Dlib,self.Alib,self.row_map = build_binned_mat( self.binner, library_collection, 
                                                                 self.spectra, col_or_row='row' )
            self.Dlib = self.Dlib.diagonal().reshape((-1,1))
            self.blksz = blksz

    def __call__( self, search_spectra, search_subset=None ):
        self.search_spectra = search_spectra
        if None == search_subset:
            self.search_subset = search_spectra.spectra()
            self.minmass = search_spectra.minMZ()
            self.maxmass = search_spectra.maxMZ()
        else:
            self.search_subset = search_subset
            massl = [search_spectra.prec_mass(p) for p in search_subset]
            self.minmass = min(massl)
            self.maxmass = max(massl)

        self.Dq,self.Aq,self.col_map = build_binned_mat( self.binner, search_spectra, 
                                                         self.search_subset,
                                                         col_or_row='row' )
        self.Dq = self.Dq.diagonal().reshape((-1,1))
        return self

    def __iter__( self ):
        self.idx = 0
        return self

    def __next__( self ):
        if self.idx >= len(self.search_subset):
            raise StopIteration

        # mask off spectra we were not supposed to be searching
        precmass = self.search_spectra.prec_mass(self.search_subset[self.idx])
        lb = bisect.bisect_left( self.med_masses, precmass - self.tol )
        ub = bisect.bisect_right( self.med_masses, precmass + self.tol )
        if lb == ub or self.empty:
            self.idx += 1
            return (self.search_spectra.idx2name(self.search_subset[self.idx-1]),[],[],[])
        else:
            ranking = np.array((self.Alib[lb:ub,:]*(self.Aq[self.idx,:].T.multiply(self.Dq[self.idx]))).multiply(self.Dlib[lb:ub]).todense())

            # only the last n items
            perm = ranking.argsort(0)[-self.n:,0].reshape((-1,))

            reverse = np.arange(min(self.n-1,perm.shape[0]-1),-1,-1)
            self.idx += 1
            return (self.search_spectra.idx2name(self.search_subset[self.idx-1]),
                    lb+perm[reverse],
                    ranking[perm[reverse]],
                    [self.row_map[i+lb] for i in perm[reverse].flat])

class BlockEagerBinnedSearch:
    """
    Binned search that lazy-evaluates query spectra.
    """
    def __init__( self, library_collection, max_mz, binsz, tol, n,
                  blksz=128, search_subset=None, pmass_window=None ):
        self.search = BinnedSearch( library_collection, max_mz, binsz, tol, n, 
                                    search_subset, pmass_window )
        self.blksz = blksz

    def __call__( self, search_spectra, search_subset=None ):
        self.search_spectra = search_spectra
        if None == search_subset:
            self.search_subset = search_spectra.spectra()
            self.minmass = search_spectra.minMZ()
            self.maxmass = search_spectra.maxMZ()
        else:
            self.search_subset = search_subset
            massl = [search_spectra.prec_mass(p) for p in search_subset]
            self.minmass = min(massl)
            self.maxmass = max(massl)
        self.search_subset.sort(key=lambda p: self.search_spectra.prec_mass(p))
        return self

    def __iter__( self ):
        self.next = iter(self.search(self.search_spectra,self.search_subset[:self.blksz])) 
        self.idx = 0
        return self

    def __next__( self ):
        if self.idx >= len(self.search_subset):
            raise StopIteration

        try:
            self.idx += 1
            return self.next.__next__()
        except StopIteration:
            self.next = iter(self.search(self.search_spectra,
                                         self.search_subset[self.idx:self.idx+self.blksz]))             
            self.idx += 1
            return self.next.__next__()

