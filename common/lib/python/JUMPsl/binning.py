import numpy as np

class Binner:
    """
    Basic class for binning m/z,intensity data.  
    """
    def __init__( self, binsz, maxmz, collection, tlower=180. ):
        """
        Initialize self.  

        binsz: width of the bin in daltons

        maxmz: maximum m/z of the domain

        collection: an object that implements the ..autoclass JUMPsl.spectral_data.SpectralData interface
        """
        self.binsz = binsz
        self.maxmz = maxmz
        self.nbins = int(np.ceil(self.maxmz/self.binsz))+1
        self.tlower = tlower

    def bin( self, mz, inten ):
        """
        Bin the m/z, intensity dta
        """
        inten **= .5
        inten[(mz <= self.tlower).nonzero()[0]] = 0
        
        return (np.ceil(mz/self.binsz).astype(np.int),inten)


