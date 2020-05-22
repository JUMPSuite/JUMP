import bisect

class NonUniformIntervalMap:
    def __init__( self, listOfLists ):
        self.thresholds = [len(listOfLists[0])]
        for l in listOfLists[1:]:
            self.thresholds.append(len(l) + self.thresholds[-1])
        self.lol = listOfLists

    def __len__( self ):
        return self.thresholds[-1]

    def __getitem__( self, x ):
        return self.__call__( x )

    def __call__( self, x ):
        index = bisect.bisect(self.thresholds,x)
        offset = x
        if index > 0:
            offset -= self.thresholds[index-1]
        return (index,offset,
                self.lol[index][offset])

class UniformIntervalMap:
    def __init__( self, listOfLists ):
        self.n = len(listOfLists)
        self.m = len(listOfLists[0])
        self.lol = listOfLists

    def __len__( self ):
        return self.n*self.m

    def __getitem__( self, x ):
        return self.__call__( x )

    def __call__( self, x ):
        index = x // self.m
        offset = x % self.m
        return (index,offset,self.lol[index][offset])

class UniformMatrixView:
    def __init__( self, listOfLists ):
        self.n = len(listOfLists)
        self.l = listOfLists

    def __len__( self ):
        return self.n*self.n

    def __getitem__( self, x ):
        return self.__call__( x )

    def __call__( self, x ):
        index = x // self.n
        offset = x % self.n
        return (self.l[index],self.l[offset])
    
class DualNonUniIntervalMap:
    def __init__( self, listOfLists1, listOfLists2 ):
        self.map1 = NonUniformIntervalMap( listOfLists1 ) 
        self.map2 = NonUniformIntervalMap( listOfLists2 )

    def __len__( self ):
        return len(self.map1)*len(self.map2)
        
    def __getitem__( self, x ):
        return self.__call__( x )

    def __call__( self, x ):
        idx1,off1,item1 = self.map1( x )
        idx2,off2,item2 = self.map2( off1 )
        return (item1,item2,idx1 == idx2)
