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
    return '''\n\n\ncometsearch.py parameterFile file.mzXML/s\n\n\n '''

parser = argparse.ArgumentParser(description="Search using Comet", prog="cometsearch.py",usage=msg())
parser.add_argument("parameterfile", help="comet parameter file")
parser.add_argument("mzXML",help="single or list of mzXML files",nargs='+')
parser.add_argument("--queue",dest="queue", default="queue=standard")
parser.add_argument("--mem",dest="mem", default="mem=8000")

args = parser.parse_args()
mzXMLs = args.mzXML
queue = args.queue.split("=")[-1]
mem = args.mem.split("=")[-1]

mzXML_path = os.getcwd()
#comet = "/home/spoudel1/PengProteomics/comet"
#mzXMLs = glob.glob(mzXML_path+"/*.mzXML")
comet_params = args.parameterfile

#line 38 and 42 added by suresh 07/09/2020
def fromListToFolder(filelist):
    folder = []
    for mzFile in filelist:
        dirToMake = mzXML_path+"/"+os.path.basename(mzFile).split(".")[0]
        folder.append(dirToMake)
        os.system("mkdir "+dirToMake)
        os.system("mv "+mzFile+" "+dirToMake)
        os.system("cp "+comet_params+" "+dirToMake+"/comet.params")
    return folder

#suresh added this section 07/09/2020

def waitUntilFinish(n, fol_list):
    hold=0
    while hold!=1:
        log = []
        for mzFol2 in fol_list:
            filenameMS2 = mzFol2.split("/")[-1]
            new_log = glob.glob(mzFol2+"/log/log.out")
            if len(new_log) == 0:
                continue
            else:
                log.append(new_log)
        if len(log) < n:
            time.sleep(60)
        else:
            hold = 1
    return

def create_job_file(mzxml_file):
    fileroot = mzxml_file.split("/")[-1].split(".mzXML")[0] 
    log_dir = mzXML_path+"/"+fileroot+"/log"
    old_results = mzXML_path+"/"+fileroot+"/oldsearch"
    os.system("mkdir "+log_dir)
    os.system("mkdir "+old_results)
    try:
        log_outFile = log_dir+"/log.out"
        error_outFile = log_dir+"/error.err"
        old_pepFile = mzXML_path+"/"+fileroot+"/"+fileroot+".1.pep.xml"

        cmd_rm = "rm "+log_outFile+" "+error_outFile
        cmd_mv = "mv "+old_pepFile+" "+old_results
        os.system(cmd_rm)
        os.system(cmd_mv)

    except:
        print ("Creating log directory\n")

    job_header = "#!/bin/bash\n#BSUB -P TestComet\n#BSUB -J comet\n#BSUB -oo "+log_dir+"/log.out\n#BSUB -eo "+log_dir+"/error.err\n#BSUB -n 8\n"
    #activate this if you need 2018 comet
    job_body1 = "comet -P"+comet_params+" "+mzxml_file
    #job_body1 = comet+" -P"+comet_params+" "+mzxml_file
    
    jobfile = fileroot+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+job_body1)
    return jobfile


def submit_job(jobf,queue,mem):
  cmd = 'bsub -q '+queue+' -R "rusage[mem='+mem+']" < '+jobf
  os.system(cmd)
  #os.system('bsub -q gpu -R "rusage[mem=8000]" < '+jobf)

def updateCometWork(new_folders):
    done_job = []
    running_job = []
    for mzFol2 in new_folders:
        filenameMS2 = mzFol2.split("/")[-1]

        new_log = glob.glob(mzFol2+"/log/log.out")        
        pepXmlFormed = glob.glob(mzFol2+"/*.pep.xml")

        if len(pepXmlFormed) != 0:
            running_job.append(pepXmlFormed[0])
            print ("Comet is searching searching a total of ",len(running_job),"respective pepXML files are being written\n")
         
        if len(new_log) != 0:
            done_job.append(new_log[0])
            print ("Comet has finished searching a total of ",len(done_job),"respective pepXML files are written\n")

            print ("Comet is searching searching a total of ",len(running_job),"respective pepXML files are being written\n")

    return len(done_job)



#comet_params = "comet.params"
#fromListToFolder(mzXMLs)

new_folders = fromListToFolder(mzXMLs)
new_input_mzxmls = []
for mzFol in new_folders:
    filenameMS = mzFol.split("/")[-1]
    new_ms2 = glob.glob(mzFol+"/"+filenameMS+".mzXML")[0]
    new_input_mzxmls.append(new_ms2)

#new_input_mzxmls = glob.glob(mzXML_path+"/*/*.mzXML")

if len(mzXMLs) < 40:
    for mz_file in new_input_mzxmls:
        submit_job(create_job_file(mz_file),queue,mem)
        print ("\n\nJob is submitted for COMET search on "+mz_file+"\n\nPLEASE WAIT PATIENTLY\n\n")
else:
    rounds = math.ceil(len(mzXMLs)/40.0)
    a=0
    for total_jobs in range(1,rounds+1):
        total_work = total_jobs*40
        if total_work < len(mzXMLs):
            for mz_file in new_input_mzxmls[a:total_work]:
                submit_job(create_job_file(mz_file),queue,mem)
                print ("\n\nJob is submitted for COMET search on "+mz_file+"\n\nPLEASE WAIT PATIENTLY\n\n")
        #while total_work >= len(mzXMLs[1:]):
            waitUntilFinish(total_work, new_folders)
            a = total_work
        else:
            for mz_file in new_input_mzxmls[a:len(mzXMLs)]:
                submit_job(create_job_file(mz_file),queue,mem)
                #print ("good")


hold2 = "in"
while hold2!="out":
    doneFiles =  updateCometWork(new_folders)
    if doneFiles < len(mzXMLs):
        time.sleep(60)
    else:
        print ("Comet has finished searching a total of "+str(doneFiles)+" mzXML files, respective pepXML files are being written\n\nPlease be patient\n\n")
        time.sleep(60)
        hold2 = "out"
        print ("CONGRATULATIONS !!! All searches have end\n\n\n\n")
