import binning as bi, spectral_data as sd, scipy.sparse as ss, numpy as np, bisect

def build_binned_mat( binner, collection, spectra, col_or_row='row' ):
    mzint = collection.read_spectra( min(spectra), max(spectra) )
    c,data = binner.bin( mzint.mz(), mzint.inten() )
    
    Al = ss.coo_matrix( (data,(mzint.idx(),c)), shape=(len(spectra),binner.nbins) )
    if 'row' == col_or_row:
        Al = Al.tocsr()
    else:
        Al = Al.tocsc()
    D = ss.dia_matrix( (np.power(Al.multiply(Al).sum(1).T,-.5),0), shape=(Al.shape[0],Al.shape[0]) )
        
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
    def __init__( self, library_collection, search_spectra, binsz, tol, n,
                  blksz=256 ):
        self.binner = bi.Binner( binsz, max(library_collection.maxMZ(),
                                          search_spectra.maxMZ()), library_collection )
        self.minmass = search_spectra.minMZ()
        self.maxmass = search_spectra.maxMZ()
        self.spectra,self.med_masses = library_collection.window_by_precmass( self.minmass, 
                                                                              self.maxmass, tol )
        self.med_masses = np.array(self.med_masses)
        self.Dlib,self.Alib,self.row_map = build_binned_mat( self.binner, library_collection, 
                                                             self.spectra, col_or_row='row' )
        self.Dlib = self.Dlib.diagonal()

        self.n = n
        self.tol = tol

        self.search_spectra = search_spectra
        self.Dq,self.Aq,self.col_map = build_binned_mat( self.binner, search_spectra, 
                                                         search_spectra.spectra(),
                                                         col_or_row='row' )
        self.Dq = self.Dq.diagonal()
        self.blksz = blksz

    def __iter__( self ):
        self.idx = 0
        return self

    def __next__( self ):
        if self.idx >= len(self.search_spectra.spectra()):
            raise StopIteration

        # mask off spectra we were not supposed to be searching
        precmass = self.search_spectra.prec_mass(self.idx)
        lb = bisect.bisect_left( self.med_masses, precmass - self.tol )
        ub = bisect.bisect_right( self.med_masses, precmass + self.tol )
        if lb == ub:
            self.idx += 1
            return ([],[],[])
        else:
            ranking = np.array((self.Alib[lb:ub,:]*(self.Aq[self.idx,:].T.multiply(self.Dq[self.idx]))).multiply(self.Dlib[self.idx]).todense())

            # only the last n items
            perm = ranking.argsort(0)[-self.n:,0].reshape((-1,))

            reverse = np.arange(min(self.n-1,perm.shape[0]-1),-1,-1)
            self.idx += 1
            return (lb+perm[reverse],
                    ranking[perm[reverse]],
                    [self.row_map[i+lb] for i in perm[reverse].flat])
