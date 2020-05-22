import spectral_data as sd, heapq, numpy as np, h5py, tempfile, os.path

class OutOfCoreMerger:
    def __init__( self, bufsz=1048576 ):
        self.tdir = tempfile.TemporaryDirectory()
        self.h5file = h5py.File(os.path.join(self.tdir.name,'_cache.h5'),'w')
        self.data_buf = []
        self.bufsz = bufsz

    def push_record( self, idxval, array_data, kvdata ):
        if len(self.data_buf) >= self.bufsz:
            self.flush()
        self.data_buf.append( (idxval,array_data,kvdata) )
    
    def flush( self ):
        g = self.h5file.create_group(str(len(self.h5file)))
        self.data_buf.sort(key = lambda t: t[0])
        g.create_dataset( 'idx', data=[t[0] for t in self.data_buf] )
        
        keys = self.data_buf[0][1].keys()
        for k in keys:
            gs = g.create_group(k)
            gs.create_dataset( 'data', data=np.hstack([self.data_buf[i][1][k] for i in range(len(self.data_buf))]) )
            accum = 0
            offsets = []
            for i in range(len(self.data_buf)):
                offsets.append( accum )
                accum += self.data_buf[i][1][k].shape[0]

            offsets.append( accum )
            gs.create_dataset( 'offsets', data=offsets )

        for i in range(len(self.data_buf)):
            g.attrs[str(i)] = repr(self.data_buf[i][2])

        self.data_buf = []

    def idx_value( self, i, j ):
        return self.h5file[str(i)]['idx'][j]

    def nitems( self, i ):
        return self.h5file[str(i)]['idx'].shape[0]

    def read( self, i, j ):
        arr_data = {}
        for k in self.h5file[str(i)]:
            if k != 'idx':
                start = self.h5file[str(i)][k]['offsets'][j]
                end = self.h5file[str(i)][k]['offsets'][j+1]
                arr_data[k] = np.array(self.h5file[str(i)][k]['data'][start:end])

        kvdata = eval(self.h5file[str(i)].attrs[str(j)])
        return (arr_data,kvdata)

    def __iter__( self ):
        self.flush()

        for i in range(len(self.h5file)):
            heapq.heappush(self.data_buf,(self.idx_value(i,0),i,0))

        return self

    def __next__( self ):
        if len(self.data_buf) == 0:
            raise StopIteration

        t = heapq.heappop(self.data_buf)
        j = t[-1]
        if j + 1 < self.nitems(t[1]):
            heapq.heappush(self.data_buf,(self.idx_value(t[1],j+1),t[1],j+1))

        arr_data,kvdata = self.read(t[1],j)
        return (t[0],arr_data,kvdata)
