#!/usr/bin/env python3
import glob
import subprocess
import fileinput
import re
import os
import sys
import time
import math
import argparse
#mzXML_path = sys.argv[1]  #give the fiel path here
#comet_params = sys.argv[2]  #give the path for parameter file

def msg(name=None):                                                            
    return '''\n\n\njump -comet parameterFile file.mzXML/s\n\n\n '''

parser = argparse.ArgumentParser(description="Search using Comet", prog="cometsearch.py",usage=msg())
parser.add_argument("parameterfile", help="comet parameter file")
parser.add_argument("mzXML",help="single or list of mzXML files",nargs='+')
args = parser.parse_args()

mzXMLs = args.mzXML

mzXML_path = os.getcwd()

#mzXMLs = glob.glob(mzXML_path+"/*.mzXML")
comet_params = args.parameterfile

def fromListToFolder(filelist):
    folder = []
    for mzFile in filelist:
        dirToMake = mzXML_path+"/"+os.path.basename(mzFile).split(".")[0]
        os.system("mkdir "+dirToMake)
        os.system("mv "+mzFile+" "+dirToMake)
        os.system("cp "+comet_params+" "+dirToMake+"/comet.params")

def waitUntilFinish(n):
    hold = 0
    while hold!=1: 
        if len(glob.glob(mzXML_path+"/*/log/log.out")) < n:
            time.sleep(60)
        else:
            hold = 1
    return

def create_job_file(mzxml_file):
    fileroot = mzxml_file.split("/")[-2]+"_"+mzxml_file.split("/")[-1].split(".mzXML")[0] 
    log_dir = mzXML_path+"/"+os.path.basename(mzxml_file).split(".")[0]+"/log"
    os.system("mkdir "+log_dir)
    job_header = "#!/bin/bash\n#BSUB -P TestComet\n#BSUB -J comet\n#BSUB -oo "+log_dir+"/log.out\n#BSUB -eo "+log_dir+"/error.err\n#BSUB -n 8\n#BSUB -N spoudel1@stjude.org\n"
    job_body1 = "comet -P"+comet_params+" "+mzxml_file
    
    jobfile = fileroot+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+job_body1)
    return jobfile


def submit_job(jobf):
  os.system('bsub -q gpu_rhel7 -R "rusage[mem=8000]" < '+jobf)


#comet_params = "comet.params"
fromListToFolder(mzXMLs)

new_input_mzxmls = glob.glob(mzXML_path+"/*/*.mzXML")

if len(mzXMLs) < 40:
    for mz_file in new_input_mzxmls:
        submit_job(create_job_file(mz_file))
        print ("\n\nJob is submitted for COMET search on "+mz_file+"\n\nPLEASE WAIT PATIENTLY\n\n")
else:
    rounds = math.ceil(len(mzXMLs)/40.0)
    a=0
    for total_jobs in range(1,rounds+1):
        total_work = total_jobs*40
        for mz_file in new_input_mzxmls[a:total_work]:
            submit_job(create_job_file(mz_file))
            print ("\n\nJob is submitted for COMET search on "+mz_file+"\n\nPLEASE WAIT PATIENTLY\n\n")
        #while total_work >= len(mzXMLs[1:]):
        waitUntilFinish(total_work)
        a = total_work

hold2 = "in"
while hold2!="out":
    cnt = 0
    log = glob.glob(mzXML_path+"/*/log/log.out")
    if len(log) < len(mzXMLs):
        time.sleep(60)
        if len(log) == 0:
            print ("Comet has not completed a single file, please be patient\n\n")
        else:
            files_searched = []
            for logs in log:    
                done_file = logs.split("/")[-3]+".mzXML"
                print ("Comet has finished searching "+done_file+"\nA total of "+str(len(log))+" mzXML files are searched, respective pepXML files are being written\n\nPlease be patient\n\n")
    else:
        print ("Comet has finished searching a total of "+str(len(log))+" mzXML files, respective pepXML files are being written\n\nPlease be patient\n\n")
        time.sleep(60)
        hold2 = "out"
        print ("CONGRATULATIONS !!! All searches have end\nPlease Wait for 5 more minutes for the completion of pepXML generation\n\n\n")
    #for pep_file in pep_list:
     #   fileobj=open(pep_file,"rb+")
      #  if fileobj.closed:
       #     cnt+=1    
        #    if cnt < len(mzXMLs):
         #       time.sleep(200)
          #  else:
           #     hold2 = "out"

#ios.system("mv "+comet_params+" comet.params")


#bj=open(filename,"wb+")

#if not fileobj.closed:
#    print("file is already opened")
