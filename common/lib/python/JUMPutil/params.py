import re

class Params:
    def __init__( self, file ):
        for l in open(file,'r').read().strip().split('\n'):
            l = re.sub('#.*','',l).replace('\r','').strip()
            k,v = re.split('\\s*=\\s*',l)
            if v.isnumeric():
                setattr( self, k, int(v) )
            else:
                try:
                    setattr( self, k, float(v) )
                except ValueError:
                    setattr( self, k, v )

