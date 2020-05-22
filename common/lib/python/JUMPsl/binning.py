import numpy as np

class Binner:
    def __init__( self, binsz, maxmz, collection, tlower=180. ):
        self.binsz = binsz
        self.maxmz = maxmz
        self.nbins = int(np.ceil(self.maxmz/self.binsz))+1
        self.tlower = tlower

    def bin( self, mz, inten ):
        inten **= .5
        inten[(mz <= self.tlower).nonzero()[0]] = 0
        
        return (np.ceil(mz/self.binsz).astype(np.int),inten)


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
        
