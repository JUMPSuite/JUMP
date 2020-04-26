import binning as bi, spectral_data as sd, scipy.sparse as ss, numpy as np

def build_binned_mat( binner, collection, spectra, do_names=True ):
    Al = []
    dl = []
    for (cols,data),idx in zip(binner.bin_iter(spectra),spectra):
        Al.append(ss.coo_matrix( (data,(np.zeros(cols.shape),cols)), shape=(1,binner.nbins) ))
        if do_names:
            dl.append( collection.read_attrs( idx )['name'] )
        else:
            dl.append( idx )

    Al = ss.vstack(Al).tocsr()
    D = ss.dia_matrix( (np.power(Al.multiply(Al).sum(1).T,-.5),0), shape=(Al.shape[0],Al.shape[0]) )
        
    return D*Al,dl

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
    def __init__( self, library_collection, search_spectra, binsz, tol, n ):
        self.binner = bi.Binner( binsz, max(library_collection.maxMZ(),
                                          search_spectra.maxMZ()), library_collection )
        self.minmass = search_spectra.minMZ()
        self.maxmass = search_spectra.maxMZ()
        self.spectra,self.med_masses = library_collection.window_by_precmass( self.minmass, 
                                                                              self.maxmass, tol )
        self.med_masses = np.array(self.med_masses)
        self.libA,self.names = build_binned_mat(self.binner, library_collection, 
                                                self.spectra, do_names=False )
        self.libA = self.libA.tocsr()
        self.idx = 0
        self.search_spectra = search_spectra
        self.n = n
        self.tol = tol

    def __iter__( self ):
        return self

    def __next__( self ):
        if self.idx > self.n:
            raise StopIteration

        spectra_idx = ss.csc_matrix(build_binned_mat(self.binner, self.search_spectra, (self.idx,),
                                                     do_names=False)[0])
        ranking = np.array((self.libA*spectra_idx.T).todense())

        # mask off spectra we were not supposed to be searching
        medmass_idx = self.med_masses[self.idx]
        mask = (self.med_masses < self.minmass - self.tol) + (self.med_masses > self.maxmass + self.tol)
        ranking[np.where(mask)[0],:] = 0.

        # only the last n items
        perm = ranking.argsort(0)[-self.n:,0].reshape((-1,))

        self.idx += 1
        reverse = np.arange(self.n-1,-1,-1)
        return (perm[reverse],
                ranking[perm[reverse]],
                [self.names[i] for i in perm[reverse].flat])
