import subprocess

class Config:
    @staticmethod
    def get( key ):
        return subprocess.check_output(['jump','-config',key])

    @staticmethod
    def put( key, value ):
        subprocess.call('jump -config {}:={}'.format(key,value), shell=True, stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL )


