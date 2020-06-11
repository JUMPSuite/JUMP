import subprocess, os

class JobManager:
    """
    Interface into JUMP's job management system
    """
    def __enter__( self ):
        self.jobs = []
        return self

    def __exit__( self, exc_type, exc_value, exc_traceback ):
        pipe = subprocess.Popen('jump -jobManager',stdin=subprocess.PIPE, stdout=None, stderr=None, env=os.environ, shell=True, text=True)
        pipe.communicate(input=str.join('\n',[str.join(' ',t) for t in self.jobs]))

    def add_job( self, tool_type, cmds ):
        self.jobs.append( (tool_type,cmds) )
