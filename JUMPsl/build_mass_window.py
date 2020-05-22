import JUMPsl.library_resampling as lr, JUMPsl.spectral_data as sd, glob, os, os.path, numpy.random as nr, numpy as np, numpy.linalg as nl, sys, scipy.io as sio, os, torch as tc, scipy.optimize as so, JUMPsl.binning as bi, scipy.sparse as ss, sys, h5py, time
from matplotlib.pyplot import *

def median_prec_mass( rdr, peptide ):
    return np.median([rdr.prec_mass(p) for p in rdr.spectra_by_peptide(peptide)])

def applyMagic( A ):
    A = A.multiply(A > 1e-4)
#    A.data **= .5
    return A

def getBpred( libCollection, trainCollection, binsz, el ):
    libBinner = bi.Binner( binsz, max(libCollection.maxMZ(),trainCollection.maxMZ()), libCollection )
    trainBinner = bi.Binner( binsz, max(libCollection.maxMZ(),trainCollection.maxMZ()), trainCollection )
    Dc = ss.dia_matrix( (np.vstack((np.ones((1,libBinner.nbins)),np.ones((1,libBinner.nbins))*.5,
                                    np.ones((1,libBinner.nbins))*.5)),(0,1,-1)),
                         (libBinner.nbins,libBinner.nbins) )

    Aelb = applyMagic(buildBinnedMatFromH5Reader( libBinner, trainCollection,
                                                  [s for (s,t),v in el] )[0].tocsr())*Dc
    Aerb = applyMagic(buildBinnedMatFromH5Reader( trainBinner, trainCollection,
                                                  [t for (s,t),v in el] )[0].tocsc())*Dc
    Dl = ss.dia_matrix( (1./np.sqrt(Aelb.multiply(Aelb).sum(1)).T,0), (Aelb.shape[0],Aelb.shape[0]) )
    Dr = ss.dia_matrix( (1./np.sqrt(Aerb.multiply(Aerb).sum(1)).T,0), (Aerb.shape[0],Aerb.shape[0]) )
    bpred = np.array(Dr*Dl*Aelb.multiply(Aerb).sum(1))[:,0]
    return bpred

def evalBinned( libCollection, trainCollection, binsz, trainEl, evalEl ):
    
    b = np.array([v for ((s,t),v) in trainEl])
    be = np.array([v for ((s,t),v) in evalEl])
    taub = optimizeCutoff( getBpred( libCollection, trainCollection, binsz, trainEl ),
                           b )
    bpred = getBpred( libCollection, trainCollection, binsz, evalEl )

    bfp = ((bpred >= taub)*(be == 0)).sum()
    bfn = ((bpred < taub)*(be != 0)).sum()
        
    return bpred,taub,bfp,bfn
    
def optimizeTSVDLeastSqClassifier( A, b, G, bn ):
    s = time.time()
    U,S,Vt = tc.svd(G)
    e = time.time()
    print( 'SVD time is {}'.format(e - s) )
    s = time.time()
    res = so.minimize_scalar( lambda alpha: crossValidateMC( A, lambda B: pinv(U, S, Vt, 10**-alpha, tc.transpose(tc.matmul(tc.transpose(B,0,1),A),0,1)), b, maxiter=2 ),
                                  method='bounded', bounds=[1,5], options={'disp':3}, tol=1e-3 )
    e = time.time()
    print( 'MC time is {}'.format(e - s) )

    x = pinv(U,S,Vt,10**-res.x,bn)
    tau = optimizeCutoff( tc.matmul(A,x), b )
    return (x,tau)

def optimizeCutoff( pred, b ):
    err = lambda tau: ((pred >= tau)*(b == 0)).sum().item() + ((pred < tau)*(b != 0)).sum().item()
    res = so.minimize_scalar( err, method='bounded', bounds=[0,1], options={'disp':3}, tol=1e-3 )
    return res.x

def buildBinnedMatFromH5Reader( binner, collection, peptides ):
    r = []
    c = []
    v = []
    dl = []
    for i,idx in enumerate(peptides):
        cv = binner( idx )
        c.append( cv[0] )
        v.append( cv[1] )
        r.append( np.ones(cv[0].shape[0])*i )
        dl.append( collection.get_name( idx ) )
        
    return (ss.coo_matrix( (np.hstack(v),(np.hstack(r),np.hstack(c))), shape=(len(dl),binner.nbins) ),
            dl)

def pinv( U, S, Vt, tol, b ):
    S = S.reshape((-1,1))
    return tc.matmul(Vt,(S**-1*(S > (tol*S[0]).item()).to(tc.double))*tc.matmul(tc.transpose(U,0,1),b))

def tinv( U, S, Vt, alpha, b ):
    S = S.reshape((-1,1))
    return tc.matmul(Vt,((S + alpha)**-1)*tc.matmul(tc.transpose(U,0,1),b))

def crossValidateMC( A, genX, b, tol=1e-3, maxiter=8 ):
    errLast = 0
    err = [0]
    done = False
    itern = 0
    blksz = 2500
    bnz =  b.nonzero()[:,0]
    nnz = bnz.shape[0]
    c = tc.transpose(tc.LongTensor().new_ones((nnz//2,blksz),device=A.device)*tc.arange(blksz,dtype=tc.long,device=A.device),0,1).reshape((1,-1))
#    for i in range(blksz): c[i*blksz:(i+1)*blksz] *= i
    
    while not done and itern < maxiter:
        # r = tc.cat([bnz[nr.choice(nnz,nnz//2,replace=False)]
        #             for i in range(blksz)]).reshape((1,-1))
        r = bnz[tc.randint(0,nnz,(1,blksz*(nnz//2)),dtype=tc.long,device=A.device)]
            
        B = tc.sparse.DoubleTensor(tc.cat((r,c),0),
                                   b[bnz[0],0]**2*tc.DoubleTensor().new_ones((nnz//2*blksz,),device=A.device),
                                   (A.shape[0],blksz),device=A.device)
        AX = tc.matmul(A,genX(B))
        R = tc.div(tc.nn.functional.hardtanh(AX,min_val=0,max_val=1)+1,2.) - (b/b.max())
        R[r,c] = 0
        err.append( R.abs_().sum().item() )
        if np.abs(np.mean(err) - errLast) < tol:
            done = True
        else: 
            errLast = np.mean(err)
        itern += 1
            
    return err[-1]
    
def splitList( el, t=.5 ):
    trainEl = []
    evalEl = []
    for i in range(0,len(el)):
        if nr.uniform(0,1) > t:
            evalEl.append( el.pop() )
        else:
            trainEl.append( el.pop() )

    return (trainEl,evalEl)

if __name__ == '__main__':
    rcParams['interactive'] = sys.flags.interactive
    nsamples = 61
    seed = 2
    nedges = 1

    nr.seed(seed)
    libCollection = sd.CSRSpectralDataReader( os.path.join(os.environ['BINROOT'],
                                                           'best_replicate_csr.h5') )
    trainCollection = sd.CSRSpectralDataReader( os.path.join(os.environ['BINROOT'],
                                                              
                                                         'pooled_csr.h5') )

    libPeptides = set(libCollection.peptides()).intersection(trainCollection.peptides())
    medMass = float(sys.argv[2])
    twindow_idx,twindow_mass = trainCollection.window_by_precmass(medMass,medMass,5)
    tpepNameSet = set([trainCollection.idx2name(i) for i in twindow_idx])
    lwindow_idx,lwindow_mass = libCollection.window_by_precmass(medMass,medMass,5)
    lpepNameSet = set([libCollection.idx2name(i) for i in lwindow_idx])
    pepNameSet = sorted(list(tpepNameSet.intersection(lpepNameSet)))

    assert len(pepNameSet) > 128
#    validationPeptides = set([p for p in pepNameSet if np.abs(median_prec_mass(trainCollection,p) - medMass) < .125])
    # pepSet = set()
    # for nom in pepNameSet:
    #     pepSet.update(set([p for p in trainCollection.nindex[nom]]))
    #     pepSet.update(set([p for p in libCollection.nindex[nom]]))
    
    # trainCollection.shrink(pepSet)

    # libCollection.shrink(pepSet.intersection(myLibPeps))

    def cos( x, y ):
        x = x.todense()
        y = y.todense()
        return np.dot(x,y)[0,0]/(nl.norm(x)*nl.norm(y))

    xDict = {}
    errDict = {}
#    AlibBinned = ss.lil_matrix((AlibFourier.shape[0],libBinner.nbins))
    dimMap = []
    print( 'processing {} peptides'.format(len(pepNameSet)) )
    pepLib = libCollection
    # make sure that the library has only one peptide
    #    pepLib.shrink([pepId])
    
    print( 'building edge list for...', end='', flush=True )
    el = lr.buildAsymEdgeListN( nsamples, pepLib, trainCollection, pepNameSet, pepNameSet, nedges,
                                 seed=seed+len(xDict) )
    print( str.format( 'edge list has {} edges', len(el) ),
           flush=True )
    assert np.all([int(libCollection.idx2name(s) ==
                       trainCollection.idx2name(t)) == val for (s,t),val in el])

    trainEl,evalEl = splitList( el, .8 )
    nTrainPos = sum([e[1] for e in trainEl])
    nEvalPos = sum([e[1] for e in evalEl])
    # while nTrainPos / float(nTrainPos + nEvalPos) < .4 or nTrainPos / float(nTrainPos + nEvalPos) > .6:
    #     print( 'failed balanced split with {}, trying again...'.format(nTrainPos / float(nTrainPos + nEvalPos)) )
    #     trainEl,evalEl = splitList( trainEl + evalEl )
    #     nTrainPos = sum([e[1] for e in trainEl])
    #     nEvalPos = sum([e[1] for e in evalEl])

    print( 'training mat has {} positive instances, eval has {}.'.format(sum([e[1] for e
                                                                              in trainEl]),
                                                                         sum([e[1] for e
                                                                              in evalEl])) )
    
    print( 'building training matrix...', flush=True, end='' )
    Al,dl = lr.buildLibMatFromH5Reader( nsamples, libCollection, [s for ((s,t),v) in trainEl] )
    Ar,dr = lr.buildLibMatFromH5Reader( nsamples, trainCollection, [t for ((s,t),v) in trainEl] )
    Al = lr.pointMultiply(Ar,lr.conjugate(Al))
    del Ar
    mu = Al.mean(0)
    Al -= mu
    b = tc.DoubleTensor([v for ((s,t),v) in trainEl]).reshape((-1,1)).to(Al.device)


    G,bn = lr.normalEq( Al, b, bias=1 )#np.sqrt(len(trainEl)/
                                            #    float(sum([e[1] for e
                                             #              in trainEl]))) )
    print( 'done.' )

    
    x,tau = optimizeTSVDLeastSqClassifier( Al, b, G, bn )

    outf = h5py.File( sys.argv[1], 'w' )
    outf.create_group(str(medMass)).create_group('fourier').create_dataset( 'x', data=np.array(x.cpu()) )
    outf[str(medMass)]['fourier'].attrs['tau'] = tau
    outf[str(medMass)].create_group('peptides')

    binsz = 1.
#     pepNameIndices = [libCollection.get_index(q,0) for q in pepNameSet]
#     libBinner = bi.Binner( binsz, max(libCollection.maxMZ(),trainCollection.maxMZ()), libCollection )
#     trainBinner = bi.Binner( binsz, max(libCollection.maxMZ(),trainCollection.maxMZ()), trainCollection )
#     Dc = ss.dia_matrix( (np.vstack((np.ones((1,libBinner.nbins)),np.ones((1,libBinner.nbins))*.5,
#                                     np.ones((1,libBinner.nbins))*.5)),(0,1,-1)),
#                          (libBinner.nbins,libBinner.nbins) )

#     Aelb = applyMagic(buildBinnedMatFromH5Reader( libBinner, trainCollection,
#                                                   pepNameIndices )[0].tocsr())*Dc
#     Dl = ss.dia_matrix( (1./np.sqrt(Aelb.multiply(Aelb).sum(1)).T,0), (Aelb.shape[0],Aelb.shape[0]) )

#     Ael,dll = lr.buildLibMatFromH5Reader( nsamples, libCollection, [s for ((s,t),v) in evalEl], blocksize=20 )
#     Aer,dlr = lr.buildLibMatFromH5Reader( nsamples, trainCollection, [t for ((s,t),v) in evalEl], blocksize=20 )
#     #    Aerb = buildBinnedMatFromH5Reader( trainBinner, trainCollection, [t for ((s,t),v) in evalEl] )[0]
#     be = tc.DoubleTensor([v for ((s,t),v) in evalEl]).to(Al.device)
#     pred = tc.matmul(lr.pointMultiply(Aer,lr.conjugate(Ael))-mu,x)[:,0]
#     fp = ((pred >= tau)*(be == 0)).sum().item()
#     fn = ((pred < tau)*(be != 0)).sum().item()
#     bpred,taub,bfp,bfn = evalBinned( libCollection, trainCollection, binsz, trainEl, evalEl )
#     outf[str(medMass)].create_dataset( 'pred', data=np.array(pred.cpu()) ).attrs['tau'] = tau
#     outf[str(medMass)].create_dataset( 'bpred', data=np.array(bpred) ).attrs['tau'] = taub

#     for p in sorted(validationPeptides):
#         outf[str(medMass)]['peptides'].create_group(sd.remap(p))
#         for t in sorted(set([e[0][1] for e in evalEl if trainCollection.idx2name(e[0][1]) == p])):
#             print( 'building evaluation matrix for {}-{}...'.format(p,t), flush=True, end='' )
#             dataGroup = outf[str(medMass)]['peptides'][sd.remap(p)].create_group(str(t))

#             el = [((s,t),libCollection.idx2name(s) == trainCollection.idx2name(t))
#                   for s in pepNameIndices]
                  
#             Ael,dll = lr.buildLibMatFromH5Reader( nsamples, libCollection, [s for ((s,t),v) in el], blocksize=20 )
#             Aer,dlr = lr.buildLibMatFromH5Reader( nsamples, trainCollection, [t for ((s,t),v) in el], blocksize=20 )
#             #    Aerb = buildBinnedMatFromH5Reader( trainBinner, trainCollection, [t for ((s,t),v) in evalEl] )[0]
#             be = tc.DoubleTensor([v for ((s,t),v) in el]).to(Al.device)
#             dataGroup.create_dataset( 'Ael', data=np.array(Ael.cpu()) ).attrs['dll'] = str.join( '\n', dll )
#             dataGroup.create_dataset( 'Aer', data=np.array(Ael.cpu()) ).attrs['dlr'] = str.join( '\n', dlr )
#             dataGroup.create_dataset( 'be', data=np.array(be.cpu()) )
#             print( 'done.' )
    
        
#             #     Aelx = lr.pointMultiply(lr.conjugate(Ael),tc.transpose(x,0,1))
#             pred = tc.matmul(lr.pointMultiply(Aer,lr.conjugate(Ael))-mu,x)[:,0]
#             dataGroup.create_dataset( 'pred', data=np.array(pred.cpu()) )

            
#             Aerb = applyMagic(buildBinnedMatFromH5Reader( trainBinner, trainCollection,
#                                                           [t for (s,t),v in el] )[0].tocsc())*Dc
#             Dr = ss.dia_matrix( (1./np.sqrt(Aerb.multiply(Aerb).sum(1)).T,0), (Aerb.shape[0],Aerb.shape[0]) )
#             bpred = np.array(Dr*Dl*Aelb.multiply(Aerb).sum(1))[:,0]
#             dataGroup.create_dataset( 'bpred', data=bpred )
            
#             pred2 = tc.matmul(Al,x)

#             # figure(1)
#             # hist(np.array(pred2.cpu())[np.array(1-b.cpu()).nonzero()[0]],bins=20,histtype='step')
#             # hist(np.array(pred2.cpu())[np.array(b.cpu()).nonzero()[0]],bins=20,histtype='step')
#             # xlabel('score')
#             # ylabel('density')
#             # title('Fourier: training set')
#             # savefig('ftrain-{}-{}.png'.format(sd.remap(p),t))

#             # figure(2)
#             # hist(np.array(pred.cpu())[np.array(1-be.cpu()).nonzero()[0]],bins=20,histtype='step')
#             # hist(np.array(pred.cpu())[np.array(be.cpu()).nonzero()[0]],bins=20,histtype='step')
#             # xlabel('score')
#             # ylabel('density')
#             # title('Fourier: evaluation set')
#             # savefig('feval-{}-{}.png'.format(sd.remap(p),t))
            
#             # figure(3)
#             # hist(np.array(bpred)[np.array(be.cpu()).nonzero()[0]],bins=20,histtype='step')
#             # hist(np.array(bpred)[np.array(1-be.cpu()).nonzero()[0]],bins=20,histtype='step')
#             # xlabel('score')
#             # ylabel('density')
#             # title('Binning: evaluation set')
#             # savefig('beval-{}-{}.png'.format(sd.remap(p),t))

    
# #     AlibFourier[len(dimMap),:] = np.array(lr.conjugate(Ael))
# #     XlibFourier[len(dimMap),:] = np.array(x[:,0])
# #     AlibBinned[len(dimMap),:] = ss.lil_matrix(buildBinnedMatFromH5Reader( libBinner, libCollection,
# #                                                                               [s for (s,t),v in evalEl[:1]] )[0])

# # #    dimMap.append( pepName )
    
# #     print( 'fourier: {} false positives, {} false negatives out of {}'.format(fp,fn,max(be.shape)) )
# #     print( 'baseline: {} false positives, {} false negatives out of {}'.format(bfp,bfn,max(be.shape)) )
# #     xDict[pepName] = x
# #     errDict[pepName] = np.array((fp,fn))

#     # print( 'writing data to {}'.format(os.path.join(savePath,'dist-{}'.format(pepName))) )
#     # dstore = factory.newShelf( pepName )
#     # dstore.put( 'A', np.array(A) )
#     # dstore.put( 'x', np.array(x) )
#     # dstore.put( 'Aelx', np.array(Aelx) )
#     # dstore.put( 'Aer', np.array(Aer) )
#     # dstore.put( 'trainEl', trainEl )
#     # dstore.put( 'evalEl', evalEl )
#     # dstore.put( 'pred', np.array(pred) )
#     # dstore.put( 'bpred', np.array(bpred) )
#     # dstore.put( 'tau_taub', np.array((tau,taub)) )
    
#     # Aerb = Aerb.tocoo()
#     # dstore.put( 'Aerb', np.vstack( (Aerb.row.reshape((1,-1)),Aerb.col.reshape((1,-1)),
#     #                                 Aerb.data.reshape((1,-1))) ) )
#     # dstore.put( 'Aerb_shape', Aerb.shape )
                

#     # shared = factory.newShelf( 'peptideList' )
#     # shared.put( 'peptideList', pepNameSet )

#     # shared = factory.newShelf( 'errDict' )    
#     # for k,v in errDict.items():
#     #     shared.put( k, v )

#     # shared = factory.newShelf( 'AlibFourier' )
#     # shared.put( 'A', AlibFourier )
#     # shared.put( 'X', XlibFourier )
#     # shared.put( 'dimMap', dimMap )
    
#     # AlibBinned = AlibBinned.tocoo()
#     # shared = factory.newShelf( 'AlibBinned' )
#     # shared.put( 'AlibBinned', np.vstack( (AlibBinned.row.reshape((1,-1)),AlibBinned.col.reshape((1,-1)),AlibBinned.data.reshape((1,-1))) ) )
#     # shared.put( 'shape', np.array(AlibBinned.shape) )

#     # print( 'done with all {} peptides.'.format(len(pepNameSet)) )

#     # factory.putAttr( 'peptides', str.join('\n',pepNameSet) )
