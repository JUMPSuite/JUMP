import numpy as np

class Binner:
    def __init__( self, binsz, maxmz, collection, tlower=180. ):
        self.binsz = binsz
        self.maxmz = maxmz
        self.nbins = int(np.ceil(self.maxmz/self.binsz))+1
        self.collection = collection
        self.tlower = tlower

    def iscontiguous( self, l ):
        l = np.array(l)
        return (l[1:] - (l[:-1] + 1)).sum() == 0

    def bin_iter( self, sid_list ):
        if hasattr(self.collection,'block_read_peptides') and self.iscontiguous( sid_list ):
            perm = np.argsort(sid_list)
            rperm = dict([(p,i) for i,p in enumerate(perm)])
            allmz = self.collection.block_read_peptides( sid_list[perm[0]], sid_list[perm[-1]] )
            
            return [self.bin_mzint(allmz[rperm[i]]) for i in range(len(allmz))]
                    
        else:
            return [self.bin_mzint(self(sid)) for sid in sid_list]

    def __call__( self, sid ):

        mzint = self.collection.read_peptide( sid )
        return self.bin_mzint( mzint )

    def bin_mzint( self, mzint ):
        mzint[:,1] **= .5
        mzint = mzint[(mzint[:,0] > self.tlower).nonzero()[0],:]
        
        return (np.ceil(mzint[:,0]/self.binsz).astype(np.int),mzint[:,1])


class CachedBinner:
    def __init__( self, maxmz, collection, pepidxs ):
        self.maxmz = maxmz
        self.collection = collection
        self.pepidxs = pepidxs
        self.mzint = [self.collection.read_peptide(pidx) for pidx in pepidxs]
        self.c = np.vstack([np.ones(mzint.shape[0])*i for i,mzint in enumerate(self.mzint)])
        self.mzint = np.vstack(self.mzint)

    def __call__( self, binsz ):
        mzint = np.array(self.mzint)
        mzint[:,1] **= .5
        mzint = mzint[(mzint[:,0] > self.tlower).nonzero()[0],:]

        nbins = int(np.ceil(self.maxmz/binsz))+1
        return (np.ceil(self.mzint[:,0]/binsz).astype(np.int),self.c,self.mzint[:,1])
        
