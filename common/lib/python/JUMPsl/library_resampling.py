import numpy as np, numpy.fft as nft, numpy.random as nr, glob, sys, pickle, scipy.sparse as ss, numpy.linalg as nl, scipy.io as sio, os.path, JUMPsl.interval_map as im, torch as tc, warnings, time

try:
    ourDevice = tc.cuda.current_device()
except Exception:
    warnings.warn( 'CUDA not available; using CPU for math...this could be sloooooooow',
                   RuntimeWarning )
    ourDevice = tc.device('cpu')

def chebNodes( nFreq ):
    freq = np.zeros(nFreq//2 + 1)
    freq[1:] = 8*np.cos((2.*np.arange(1,nFreq//2 + 1) - 1)*(np.pi/(2*nFreq)))
    return freq

class uniformRandomFreq:
    def __init__( self ):
        self.memoize = {}

    def __call__( self, nFreq ):
        if nFreq not in self.memoize:
            self.memoize[nFreq] = np.zeros(nFreq//2 + 1)
            self.memoize[nFreq][1:] = nr.uniform(0,1,nFreq//2)

        return self.memoize[nFreq]

def choice( l, n, replace=True ):
    idxs = nr.choice(len(l),n,replace=replace)
    return [l[i] for i in idxs]

class fourierQuery:
    def __init__( self, collection, nsamples, searchPeptides ):
        self.nsamples = nsamples
        self.collection = collection
        self.searchPeptides = searchPeptides 
        self.Al,self.dl = buildLibMatFromH5Reader( nsamples, collection, searchPeptides )

    def __call__( self, queryPeptides ):
        if isinstance(queryPeptides,list):    
            Ar,dr = buildLibMatFromH5Reader( self.nsamples, self.collection, queryPeptides )
            Ar = pointMultiply(Ar,conjugate(self.Al))
            b = tc.DoubleTensor([int(self.dl[i] == dr[i]) for i in range(len(self.dl))])
            return Ar,self.dl,b
        else:
            Ar,dr = buildLibMatFromH5Reader( self.nsamples, self.collection, [queryPeptides] )
            Al = pointMultiply(self.Al,conjugate(Ar))
            b = tc.DoubleTensor([int(dr[0] == self.dl[i]) for i in range(len(self.dl))])
            return Al,self.dl,b

def buildLibMatFromH5Reader( nsamples, h5reader, peptides, blocksize=10 ):
    # A = np.empty((len(peptides),nFreq(nsamples)))
    # dl = []
    # for i,(key,idx) in enumerate(peptides):
    #     row = conjugate(toFourierDomain_(collection.readers[key].read_record(idx),nsamples))
    #     row /= hermNorm(row)
    #     A[i,:] = row
    #     dl.append( collection.readers[key].pep_map[idx]['Name'] )
    s = time.time()
    if read_contiguous:
        offset = min(peptides)
        mzint = h5reader.read_spectra(offset,max(peptides)+1)
        sl = [np.hstack((mzint.mz(pidx-offset).reshape((-1,1)),
                         mzint.inten(pidx-offset).reshape((-1,1)))) for pidx in peptides]
        dl = [h5reader.idx2name(pidx) for pidx in peptides]
    else:
        sl = []
        dl = []
        for pidx in peptides:
            mzint = h5reader.read_spectra(pidx,pidx+1)
            sl.append( np.hstack(mzint.mz().reshape((-1,1)),
                                 mzint.inten().reshape((-1,1))) )
            dl.append( h5reader.idx2name( pidx ) )

    e = time.time()
    print( 'spectra read time: {}'.format(e - s) )
    s = time.time()
    A = toFourierDomain_( sl, nsamples, blocksize )
    tc.div(tc.transpose(A,0,1),hermNorm(A),out=tc.transpose(A,0,1))
    e = time.time()
    print( 'spectra eval time: {}'.format(e - s) )

    return A,dl

def buildLibMatFromCollection( nsamples, collection, peptides ):
    # A = np.empty((len(peptides),nFreq(nsamples)))
    # dl = []
    # for i,(key,idx) in enumerate(peptides):
    #     row = conjugate(toFourierDomain_(collection.readers[key].read_record(idx),nsamples))
    #     row /= hermNorm(row)
    #     A[i,:] = row
    #     dl.append( collection.readers[key].pep_map[idx]['Name'] )
    sl = [np.array(list(collection.readers[k].read_record(idx))) for k,idx in peptides]
    dl = [collection.readers[k].pep_map[idx]['Name'] for k,idx in peptides]
    A = toFourierDomain_( sl, nsamples, 5 )
    for i in range(A.shape[0]):
        A[i,:] /= hermNorm(A[i,:]).item()
    
    return A,dl

def buildLibMat( nsamples, reader ):
    dataList = []        
    
    # things come transposed out of Epetra
    AnpView = A.ExtractView().T  
    for gd,ld in [(i,dimMap.LID(i)) for i in dimMap.MyGlobalElements()]:
#        pepName,(readerName,i) = dataList[gd]
        row = toFourierDomain(reader.read_record(keys[gd]),nsamples)
        row /= hermNorm(row)
        AnpView[ld,:] = row

    return dataList,A,dimMap

def buildAsymEdgeList( sampleRate, leftreader, rightreader, LHSpeptides,
                       nedges, seed=0, RHSpeptides=None ):
    nr.seed(seed)
#    leftPeps = leftCollection.items()
#    rightPeps = rightCollection.items()
    matchEdges = set()
    nonMatchEdges = set()
    if None == RHSpeptides:
        RHSpeptides = list(rightreader.peptides())
    # build up inter-peptide edges
    while len(nonMatchEdges) < nedges:
        lpep = LHSpeptides[nr.randint(0,len(LHSpeptides))]
        lidx = nr.randint(0,leftreader.n_spectra(lpep))
        source = leftreader.get_index(lpep,lidx)
        rpep = RHSpeptides[nr.randint(0,len(RHSpeptides))]
        ridx = nr.randint(0,leftreader.n_spectra(rpep))
        target = rightreader.get_index(rpep,ridx)
        if leftreader.idx2name(source) == rightreader.idx2name(target):
            matchEdges.add( (source,target) )
        else:
            nonMatchEdges.add( (source,target) )


    assert len(nonMatchEdges) > len(matchEdges)
    matchEdges = set([((s,t),1) for s,t in matchEdges])
    nonMatchEdges = set([((s,t),0) for s,t in nonMatchEdges])
#    del rightPeps
#    sharedPeps = list(set(leftCollection.nindex.keys()).intersection(set(rightCollection.nindex.keys())))
#    imap = im.NonUniformIntervalMap([rightCollection.nindex[t[0]] for t in leftPeps])
    # for eNo in range(len(imap)):
    #     lpid,rmpid,(readerId,readerIdx) = imap[eNo]
    #     assert leftPeps[lpid][0] == rightCollection.readers[readerId].pep_map[readerIdx]['Name']
    for pepName in LHSpeptides:
        n = leftreader.n_spectra(pepName)
        m = rightreader.n_spectra(pepName)
        if n*m <= nedges:
            for i in range(0,n):
                source = leftreader.get_index(pepName,i)
                for j in range(0,m):
                    target = rightreader.get_index(pepName,j)
                    matchEdges.add(((source,target),1))
        else:
            while len(matchEdges) < min(nedges,totIntraEdges):
                source = nr.randint(0,len(reader))
                target = nr.randint(0,len(reader))
                if reader.idx2name(source) == reader.idx2name(target):
                    matchEdges.add(((source,target),1))
        
    return list(matchEdges.union(nonMatchEdges))
            
def buildAsymEdgeListN( sampleRate, leftreader, rightreader, LHSpeptides, RHSpeptides,
                        nedges, seed=0 ):
    matchEdges = set()
    nonMatchEdges = set()
    nr.seed(seed)
    LHSMap = dict([(p,leftreader.spectra_by_peptide(p)) for p in LHSpeptides])
    RHSMap = dict([(p,rightreader.spectra_by_peptide(p)) for p in RHSpeptides])

    for pepName in LHSpeptides:
        lidx = LHSMap[pepName]
        ridx = RHSMap[pepName]
        assert len(lidx)*len(ridx) >= nedges
        tempEdges = set()
        while len(tempEdges) < nedges:
            source = lidx[nr.randint(0,len(lidx))]#leftreader.get_index(pepName,nr.randint(0,leftreader.n_spectra(pepName)))
            target = ridx[nr.randint(0,len(ridx))]#rightreader.get_index(pepName,nr.randint(0,rightreader.n_spectra(pepName)))
            tempEdges.add((source,target))
        matchEdges = matchEdges.union(tempEdges)
        
#    leftPeps = leftCollection.items()
#    rightPeps = rightCollection.items()
    # build up inter-peptide edges
    while len(nonMatchEdges) != len(matchEdges):
        lpep = LHSpeptides[nr.randint(0,len(LHSpeptides))]
        lidx = LHSMap[lpep]
        source = lidx[nr.randint(0,len(lidx))]
        rpep = RHSpeptides[nr.randint(0,len(RHSpeptides))]
        ridx = RHSMap[rpep]
        target = ridx[nr.randint(0,len(ridx))]
        if leftreader.idx2name(source) == rightreader.idx2name(target):
            matchEdges.add( (source,target) )
        else:
            nonMatchEdges.add( (source,target) )


    assert len(nonMatchEdges) == len(matchEdges)
    matchEdges = set([((s,t),1) for s,t in matchEdges])
    nonMatchEdges = set([((s,t),0) for s,t in nonMatchEdges])
#    del rightPeps
#    sharedPeps = list(set(leftCollection.nindex.keys()).intersection(set(rightCollection.nindex.keys())))
#    imap = im.NonUniformIntervalMap([rightCollection.nindex[t[0]] for t in leftPeps])
    # for eNo in range(len(imap)):
    #     lpid,rmpid,(readerId,readerIdx) = imap[eNo]
    #     assert leftPeps[lpid][0] == rightCollection.readers[readerId].pep_map[readerIdx]['Name']
    return list(matchEdges.union(nonMatchEdges))

def buildSymEdgeList( tm, idxs, nEdges ):
    el = set()
    while len(el) < nEdges:
        source,target = [idxs[i] for i in nr.randint(0,len(idxs),2)]
        val = int(tm[source] == tm[target])
        if ((source,target),val) not in el or ((target,source),val) not in el:
            el.add( ((source,target),val) )
    return list(el)
                       
def buildEdgeList( sampleRate, comm, leftCollection, rightCollection,
                   nclassSamples, nedges ):
    edgeList = []
    gid = 0
    pepEnum = list(rightCollection.nindex.keys())
    imap = im.NonUniformIntervalMap( [rightCollection.nindex[pnom] for pnom in pepEnum] )
    for nom in leftCollection.nindex.keys():
        instances = leftCollection.nindex[nom]
        if nclassSamples >= len(instances):
            sources = instances
        else:
            sources = choice(instances,nclassSamples,replace=False)
            
        if nedges >= len(rightCollection.nindex[nom]):
            for s in sources:
                for t in [v for v in rightCollection.nindex[nom] if v != s]:
                    if gid % comm.NumProc() == comm.MyPID():
                        edgeList.append( ((s,t),1) )
                    gid += 1

        else:
            for s in sources:
                for t in choice([v for v in rightCollection.nindex[nom] if v != s],
                                nedges,replace=False):
                    if gid % comm.NumProc() == comm.MyPID():
                        edgeList.append( ((s,t),1) )
                    gid += 1

        for s in sources:
            interClass = addInterClassEdgesRandom( s, nom,
                                                   min(len(leftCollection.nindex[nom])*
                                                       (len(rightCollection.nindex[nom])-1),
                                                       nedges),
                                                   imap, pepEnum )
            for i,edge in enumerate(interClass):
                if gid % comm.NumProc() == comm.MyPID():
                    edgeList.append( edge )
                gid += 1

    myElements = [i for i in range(gid) if i % comm.NumProc() == comm.MyPID()]
    bm = pe.BlockMap( -1, len(edgeList), 1, 0, comm )
    return (edgeList,bm)

def addInterClassEdgesRandom( s, nom, nedges, intervalMap, pepEnum ):
    edges = []
    for i in range(nedges):
        pepNo,offset,pep = intervalMap[nr.randint(0,len(intervalMap))]
        while pepEnum[pepNo] == nom: 
            pepNo,offset,pep = intervalMap[nr.randint(0,len(intervalMap))]

        edges.append( ((s,pep),0) )

    return edges

def addInterclassEdgesByScore( x, s, nom, names, nedges, mdata, libA, dimMap, comm ):
    Dr = ss.dia_matrix( (np.multiply(libA,np.conjugate(libA)).sum(1)**-.5,0),
                        (libA.shape[0],libA.shape[0]) )
    Dl = ss.dia_matrix( ([1.] + ([2.]*(x.shape[0]-1)),0),
                        (x.shape[0],x.shape[0]) )
    y = Dr*np.dot(libA,Dl*x)
    
    locScores = [(y[v].real,mdata[k][1:]) for k,v in dimMap.items() if mdata[k][0] != nom]
    locScores.sort(key=lambda t: t[0], reverse=True)
    globScores = sum(comm.allgather( locScores[:nedges] ),[])
    globScores.sort(key=lambda t: t[0], reverse=True)
    
    return [((s,t[1]),0) for t in globScores[:nedges]]

def createTrainingMat( edgeList, sampleRate, leftCollection, rightCollection,
                       bias=1. ):
    Al,dl = buildLibMatFromH5Reader( sampleRate, leftCollection, [s for ((s,t),v) in edgeList[:1]] )
    Ar,dr = buildLibMatFromH5Reader( sampleRate, rightCollection, [t for ((s,t),v) in edgeList] )
    # signals are hermitian, so side-ness does not matter and point multiply gets its shape
    # from the left argument
    Al = pointMultiply(Ar,conjugate(Al))
    b = tc.DoubleTensor([v for ((s,t),v) in edgeList]).reshape((-1,1)).to(ourDevice)
    idxs = b.nonzero()[:,0]
    Al[idxs,:] *= bias**.5
    b *= bias**.5

    return (Al,b)

def normalEq( A, b, bias ):
    G = tc.matmul(tc.transpose(A,0,1),A)
    maskNpGramMatrix(G)
    bNormal = tc.matmul(tc.transpose(A,0,1),b)#*(bias**2)
#    b *= bias
    
    return (G,bNormal)

def innerProductScalarMat( nsamples ):
    data = np.ones(nsamples)
    data[1:] *= 2
    data[imaginaryIndex(nsamples):] *= -1
    return data

def maskMVGramMatrix( A ):
    return maskNpGramMatrix(A.ExtractView())

def maskNpGramMatrix( A ):
    iidx = imaginaryIndex(A.shape[0])
    A[iidx:,:iidx] = 0
    A[:iidx,iidx:] = 0
    return A

fourierFrequencies = chebNodes

def nFreq( nBasisFunctions ):
    return fourierFrequencies(nBasisFunctions).shape[0]*2 - 1

def toFourierDomain_( spectralDataBlock, nBasisFunctions, blksz ):
    s = tc.DoubleTensor(2*np.pi*fourierFrequencies(nBasisFunctions)).to(ourDevice).reshape((1,1,-1))
    maxlen = max([a.shape[0] for a in spectralDataBlock])
    Dmz = np.zeros((maxlen,len(spectralDataBlock)))
    Dinten = np.zeros((maxlen,len(spectralDataBlock)))
    for i,a in enumerate(spectralDataBlock):
        Dmz[:a.shape[0],i] = a[:,0]
        Dinten[:a.shape[0],i] = a[:,1]

    A = tc.empty(len(spectralDataBlock),nFreq(nBasisFunctions)).to(tc.double).to(ourDevice)
    for i in range(0,A.shape[0]//blksz + int(A.shape[0] % blksz != 0)):
        mz = tc.DoubleTensor(Dmz[:,i*blksz:(i+1)*blksz].reshape((Dmz.shape[0],-1,1))).to(ourDevice)
        inten = tc.DoubleTensor(Dinten[:,i*blksz:(i+1)*blksz].reshape((Dinten.shape[0],-1,1))).to(ourDevice)**.5
        A[i*blksz:(i+1)*blksz,:imaginaryIndex(nBasisFunctions)] = (tc.cos(-mz*s)*inten).sum(0)
        A[i*blksz:(i+1)*blksz,imaginaryIndex(nBasisFunctions):] = (tc.sin(mz*s)*inten).sum(0)[:,1:]

    return A

def toFourierDomain( spectralData, nBasisFunctions ):
    s = 2*np.pi*fourierFrequencies(nBasisFunctions)*1j
    a = sum([inten*np.exp(-mz*s) for mz,inten in spectralData])
    return np.hstack((a.real,-a.imag[1:]))

def imaginaryIndex( lenX ):
    return lenX//2 + 1

def DCfreq( lenX ):
    return (0,)

def hermNorm( x ):
    x = x**2.
    x *= 2.
    if len(x.shape) > 1:
        x[:,DCfreq(x.shape[0])] /= 2.
        return tc.sqrt(x.sum(1))
    else:
        x[DCfreq(x.shape[0])] /= 2.
        return tc.sqrt(x.sum())

def hermNorm_( x ):
    dups = list(range(1,imaginaryIndex(x.shape[0]))) + list(range(imaginaryIndex(x.shape[0]),x.shape[0]))
    return np.sqrt((np.abs(x)**2).sum() + (np.abs(x[dups])**2).sum())

def conjugate( x ):
    conjX = tc.empty(x.shape).to(x.dtype).to(x.device)
    conjX[:] = x[:]
    if len(x.shape) > 1:
        conjX[:,imaginaryIndex(x.shape[1]):] *= -1
    else:
        conjX[imaginaryIndex(x.shape[0]):] *= -1
        
    return conjX

def hermitian2nonhermitianPerm( nBasisFunctions ):
    return (list(range((nBasisFunctions//2)+1)),list(reversed(range(1,(nBasisFunctions//2)+1))))

def realImagPerm( nBasisFunctions ):
    dcf = set(DCfreq( nBasisFunctions ))
    pr = np.array([i for i in range(imaginaryIndex( nBasisFunctions ))
                   if i not in dcf])
    pi = np.array([i for i in range(imaginaryIndex( nBasisFunctions ),nFreq( nBasisFunctions ))
                   if i not in dcf])
    return (pr,pi)

def pointMultiply( x, y ):
    if len(x.shape) < 2:
        x = x.reshape((1,-1))
    if len(y.shape) < 2:
        y = y.reshape((1,-1))

    z = tc.zeros(x.shape).to(tc.double).to(ourDevice)
    rperm,iperm = realImagPerm( x.shape[1] )
    z[:,DCfreq( x.shape[1] )] = tc.mul(x[:,DCfreq( x.shape[1] )],y[:,DCfreq( x.shape[1] )])
    z[:,rperm] = tc.mul(x[:,rperm],y[:,rperm]) - tc.mul(x[:,iperm],y[:,iperm])
    z[:,iperm] = tc.mul(x[:,rperm],y[:,iperm]) + tc.mul(x[:,iperm],y[:,rperm])
    return z

def pointMultiply_( x, y ):
    if len(x.shape) < 2:
        x = x.reshape((1,-1))
    if len(y.shape) < 2:
        y = y.reshape((1,-1))

    z = tc.zeros(x.shape).to(tc.double).to(ourDevice)
    rperm,iperm = realImagPerm( x.shape[1] )
    idxs0 = tc.arange(x.shape[0],device=ourDevice).repeat(len(DCfreq(x.shape[1])))
    idxs1 = tc.cat( [tc.ones(x.shape[0],dtype=tc.int64,device=ourDevice)*i for i in DCfreq(x.shape[1])] )
    tmp = tc.mul(x[idxs0,idxs1],y[idxs0,idxs1])
    z.index_put_( (idxs0,idxs1), tmp )
    del idxs1

    idxs0 = tc.arange(x.shape[0],device=ourDevice).repeat(len(iperm))
    idxs1r = tc.cat( [tc.ones(x.shape[0],dtype=tc.int64,device=ourDevice)*i for i in rperm] )
    idxs1i = tc.cat( [tc.ones(x.shape[0],dtype=tc.int64,device=ourDevice)*i for i in iperm] )
    tmp = tc.mul(x[idxs0,idxs1r],y[idxs0,idxs1r]) - tc.mul(x[idxs0,idxs1i],y[idxs0,idxs1i])
    z.index_put_( (idxs0,idxs1r), tmp )
    tmp = tc.mul(x[idxs0,idxs1r],y[idxs0,idxs1i]) + tc.mul(x[idxs0,idxs1i],y[idxs0,idxs1r])
    z.index_put_( (idxs0,idxs1i), tmp )    
    return z

def build_manifold_lib( nsamples, libSpectra, trainSpectra, seed, nedges,
                        savePath, doBias=False ):
    libCollection = sb.SpectralDataCollection( libSpectra )
    trainCollection = sb.SpectralDataCollection( trainSpectra )
    
    nr.seed(seed)

#     comm = mpi.COMM_WORLD
#     blkSz = len(mdata)//comm.size
#     print(blkSz)
#     if comm.rank < comm.size - 1:
#         dm = dict(zip(range(comm.rank*blkSz,(comm.rank+1)*blkSz),
#                       range(blkSz)))
#     else:
#         dm = dict(zip(range(comm.rank*blkSz,len(mdata)),
#                       range(blkSz + (len(mdata) % comm.size))))
    comm = pe.MpiComm()
    def mpiPrint( *args, **kwargs ):
        if 0 == comm.MyPID():
            print( *args, **kwargs )

    mpiPrint( 'building edge list...', end='', flush=True )
    el,elBm = buildAsymEdgeList( nsamples, comm, libCollection, trainCollection, nedges )
    mpiPrint( str.format( 'edge list has {} edges, min local size {}, max local size {}',
                          elBm.NumGlobalElements(), comm.MinAll(len(el)), comm.MaxAll(len(el)) ),
              flush=True )
    assert np.all([int(libCollection.idx2name(s) ==
                       trainCollection.idx2name(t) == val) for (s,t),val in el])

#    libCollection.shrink(set([s for ((s,t),v) in el]))
#    trainCollection.shrink(set([t for ((s,t),v) in el]))
    
    mpiPrint( 'building training mat...', end='', flush=True )
    A,b,G,bNormal = createTrainingMat( el, elBm, nsamples, comm, libCollection,
                                       trainCollection )
    mpiPrint( 'done.', flush=True )

    if 0 == comm.MyPID() and not os.path.exists(savePath):
        os.makedirs(savePath)

    comm.Barrier()
    dstore = umi.PartitionedShelf( comm, os.path.join(savePath,'dist') )
    dstore.put('A',A.ExtractView(),A.Map())
    dstore.put('edgeList',el,elBm)
    comm.Barrier()
    if 0 == comm.MyPID():
        sio.savemat( os.path.join(savePath,'normalEq.mat'),
                     {'G':G.ExtractView().T,'b':bNormal.ExtractView().T} )


def buildManifoldLibFromEdgeList( nsamples, sourceCollection, targetCollection, edgeList,
                                  blockMap, savePath ):

#     comm = mpi.COMM_WORLD
#     blkSz = len(mdata)//comm.size
#     print(blkSz)
#     if comm.rank < comm.size - 1:
#         dm = dict(zip(range(comm.rank*blkSz,(comm.rank+1)*blkSz),
#                       range(blkSz)))
#     else:
#         dm = dict(zip(range(comm.rank*blkSz,len(mdata)),
#                       range(blkSz + (len(mdata) % comm.size))))
    comm = blockMap.Comm()
    def mpiPrint( *args, **kwargs ):
        if 0 == comm.MyPID():
            print( *args, **kwargs )

    mpiPrint( 'building training mat...', end='', flush=True )
    A,b,G,bNormal = createTrainingMat( edgeList, blockMap, nsamples, comm, sourceCollection,
                                       targetCollection )
    mpiPrint( 'done.', flush=True )

    if 0 == comm.MyPID() and not os.path.exists(savePath):
        os.makedirs(savePath)

    comm.Barrier()
    dstore = umi.PartitionedShelf( comm, os.path.join(savePath,'dist') )
    dstore.put('A',A.ExtractView(),A.Map())
    dstore.put('b',A.ExtractView(),b.Map())
    dstore.put('edgeList',edgeList,blockMap)
    comm.Barrier()
    if 0 == comm.MyPID():
        sio.savemat( os.path.join(savePath,'normalEq.mat'),
                     {'G':G.ExtractView().T,'b':bNormal.ExtractView().T} )
    return A,b,G,bNormal


