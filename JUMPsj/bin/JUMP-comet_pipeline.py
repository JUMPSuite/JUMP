#!/usr/bin/env python

import collections
import sys
import re
import os
import pandas as pd
import numpy as np
from datetime import datetime
import glob
import argparse
import subprocess
import os.path
import time
import math
from os.path import dirname
import fileinput
#mzXML_path = sys.argv[1]  #give the fiel path here
#comet_params = sys.argv[2]  #give the path for parameter file

def msg(name=None):
    return '''\n\n\npython JUMP-comet.py.py jump parameterFile file.mzXML/s\n\n\n '''

parser = argparse.ArgumentParser(description="Search using JUMP and Comet", prog="JUMP-comet",usage=msg())
parser.add_argument("jump_parameterfile", help="jump parameter file")
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
jump_params = args.jump_parameterfile

out = subprocess.check_output(["comet","-p"])

comet_params = "comet.params.new"
#required for ppi matching
proton = 1.00727646677


def storeJUMPParams(paramFile):
    dict1 = {}
    with open(paramFile,"r") as f:
        for line in f:
            if (line.startswith("#") == False) and ("=" in line.strip()):
    
                te_line = line.strip().split("=")
                key = te_line[0].strip()

                if "#" in te_line[1]:
                    value = te_line[1].split("#")[0].strip()
                else:
                    value = te_line[1].strip()
                dict1[key]= value
    return dict1


# In[186]:


jumpParamsDict = storeJUMPParams(jump_params)


# In[202]:


def storeCometParams(paramFile):
    cometComments = {}
    dict1 = {}
    fileLines = []
    with open(paramFile,"r") as f:
        for line in f:
            if (line.startswith("#") == False) and ("=" in line.strip()):
                te_line = line.strip().split("=")
                key = te_line[0].strip()
                if "#" in te_line[1]:
                    te_line2 = te_line[1].split("#")
                    value1 = te_line2[0].strip()
                    valueComments = line.strip().split("#")[1]
                    cometComments[key] = "#"+valueComments
                else:
                    cometComments[key] = ""
                    value1 = te_line[1].strip()
                dict1[key]= value1
                fileLines.append(key)
            else:
                fileLines.append(line.strip())
    
    return fileLines,cometComments,dict1


# In[203]:


cometParamLines, cometComments, cometParamsDict = storeCometParams(comet_params)

if (jumpParamsDict['add_Nterm_peptide'] == '229.162932') and ("dynamic_S" in jumpParamsDict.keys()) and (jumpParamsDict["dynamic_S"] == "79.96633"):
    data = "TMThhpho"
elif (jumpParamsDict['add_Nterm_peptide'] == '229.162932'):
    data = "TMThh"
elif (jumpParamsDict['add_Nterm_peptide'] == '304.2071453') and ("dynamic_S" in jumpParamsDict.keys()) and (jumpParamsDict["dynamic_S"] == "79.96633"):
    data = "TMTpro_pho"
elif (jumpParamsDict['add_Nterm_peptide'] == '304.2071453'):
    data = "TMTpro"
elif (jumpParamsDict['mass_correction'] == '2') and (jumpParamsDict['MS2_deisotope'] == '1'):
    data = "HH"
elif (jumpParamsDict['mass_correction'] == '1') and (jumpParamsDict['MS2_deisotope'] == '0'):
    data = "HL"
else:
    print ("Please Check your JUMP parameters again to make sure you have correct parameters")


# In[226]:


#check enzyme in jump
if "Tryptic" in jumpParamsDict["enzyme_info"]:
    cometParamsDict["search_enzyme_number"] = "1"
elif "LysC" in jumpParamsDict["enzyme_info"]:
    cometParamsDict["search_enzyme_number"] = "3"
elif "ArgC" in jumpParamsDict["enzyme_info"]:
    cometParamsDict["search_enzyme_number"] = "5"
elif "GluC" in jumpParamsDict["enzyme_info"]:
    cometParamsDict["search_enzyme_number"] = "8"
else:
    print ("Invalid enzyme selection in JUMP parameters")


# In[228]:


if jumpParamsDict["digestion"] == "full":
    cometParamsDict["num_enzyme_termini"]= "2"
elif jumpParamsDict["digestion"] == "partial":
    cometParamsDict["num_enzyme_termini"] = "1"
else:
    print ("Incorrect enzyme digestion parameters for JUMP")


# In[229]:


#ion series check
#ion_series = 1 1 0 0 0 0 0 1 0									# a, b, c, d, v, w, x, y and z ions, respectively

allIons = jumpParamsDict["ion_series"].split()

cometParamsDict["use_A_ions"] = allIons[0]
cometParamsDict["use_B_ions"] = allIons[1]
cometParamsDict["use_C_ions"] = allIons[2]
cometParamsDict["use_X_ions"] = allIons[-3]
cometParamsDict["use_Y_ions"] = allIons[-2]
cometParamsDict["use_Z_ions"] = allIons[-1]


# In[230]:


if (data == "TMThhpho") or (data == "TMTpro_pho"):
    cometParamsDict["variable_mod02"] = "79.966331 STY 0 3 -1 0 0"
else:
    cometParamsDict["variable_mod02"] = "0.0 X 0 3 -1 0"


# In[234]:


if (data == "TMThhpho") or (data == "TMTpro_pho") or (data == "HH") or (data=="TMThh")  or (data == "TMTpro"):
    cometParamsDict["fragment_bin_tol"] = "0.02"              # binning to use on fragment ions
    cometParamsDict["fragment_bin_offset"] = "0.0"              # offset position to start the binning (0.0 to 1.0)
    cometParamsDict["theoretical_fragment_ions"] = "0"          # 0=use flanking peaks, 1=M peak only
    cometParamsDict["spectrum_batch_size"] = "10000"


# In[245]:


cometParamsDict["database_name"] = jumpParamsDict["database_name"][0:-4]
cometParamsDict["allowed_missed_cleavage"] = jumpParamsDict["max_mis_cleavage"]
cometParamsDict["digest_mass_range"] = jumpParamsDict["min_peptide_mass"]+" "+jumpParamsDict["max_peptide_mass"]
cometParamsDict["peptide_mass_tolerance"] = jumpParamsDict["peptide_tolerance"]
cometParamsDict["peptide_mass_units"] = jumpParamsDict["peptide_tolerance_units"]
cometParamsDict["add_Nterm_peptide"] = jumpParamsDict["add_Nterm_peptide"]
cometParamsDict["add_K_lysine"] = jumpParamsDict["add_K_Lysine"]
cometParamsDict["add_C_cysteine"] = jumpParamsDict["add_C_Cysteine"]
cometParamsDict["variable_mod01"] = jumpParamsDict["dynamic_M"]+" M 0 3 -1 0 "
cometParamsDict["max_variable_mods_in_peptide"] = jumpParamsDict["max_modif_num"]


cometParamsDict["num_threads"] = "4"


# In[248]:

rmparam = "rm "+comet_params
try:
    os.system(rmparam)
except:
    print ("No comet.params.new file generated")

cometFlyParams = "comet.params"
with open(cometFlyParams, "w") as paramC:
    for line in cometParamLines:
        if line.startswith("#"):
            paramC.write(line+"\n")
#             print (line)      
        elif line == "":
            paramC.write(line+"\n")
        
        else:
            try:
                paramC.write(line+" = "+cometParamsDict[line]+"\t"+cometComments[line]+"\n")
#                 print (line," = ",cometParamsDict[line],"\t",cometComments[line])
            except:
                paramC.write(line+"\n")
#jump -deisotope jump_ss_HH_tmt10_mouse.params mix_ratio.mzXML
#/home/spoudel1/CometDeisotopePipeline/mix_ratio/mix_ratio.1/jump.params

def waitUntilFinish(mzXML):
    hold=0
    while hold!=1:
        base = os.path.dirname(mzXML)
        sample = mzXML.split("/")[-1].split(".mzXML")[0]
        jumpParams = base+"/"+sample+"/"+sample+".1/jump.params"
        if len(glob.glob(jumpParams))==0:
            time.sleep(60)
        else:
            hold = 1

def precMZCalc(MH, z): #MH = MH+ from dta file, z = charge and proton = H+
    proton = 1.00727646677
    precmz = (float(MH)+((int(z)-1)*proton))/int(z)
    return precmz

def pepXMLOutFileConversionOld(pepxmlFile, basefile, outFileDict):
    #suffix = pepxmlFile.split(".1.")[0]
    suffix = basefile
    update_variable_before = ''
    update_variable_after = ''
    x = ''
    y = ''
    cnt = 0
 
    with fileinput.FileInput(pepxmlFile, inplace=True, backup='.cometOriginal') as f:
        for line in f:
            update_variable_before = x
            update_variable_after = y

            if "<spectrum_query spectrum=" in line:
                pattern = '<spectrum_query spectrum="('+suffix+'\.\d+\.\d+\.\d+)"'
                x=re.match(pattern, line.strip()).group(1)
                outfileSplit = x.split(".")

                if update_variable_before == x:
                    cnt+=1
                    #add_value = str(int(update_variable_after.split(".")[2])+1)
                    add_value = outFileDict[x][cnt]
                    y = outfileSplit[0]+"."+outfileSplit[1]+"."+add_value+"."+outfileSplit[-1]
                else:
                    cnt=0
                    add_value = outFileDict[x][cnt]
                    #y = outfileSplit[0]+"."+outfileSplit[1]+".1."+outfileSplit[-1]
                    y = outfileSplit[0]+"."+outfileSplit[1]+"."+add_value+"."+outfileSplit[-1]
                print (line.rstrip().replace(x,y))
            else:
                print (line.rstrip())

def mvLogsParams(mzFol2, basefile):
    pepxml_new = glob.glob(mzFol2+"/"+basefile+".pep.xml")
    # print (pepxml_new)
    # dtas = glob.glob(mzFol2+"/*.dtas")
    # sample = mzXML.split("/")[-1].split(".mzXML")[0]

    dtas = glob.glob(mzFol2+"/"+basefile+".*dtas")
    pattern =  mzFol2+"/"+basefile
    suffixes = []
    for dta in dtas:
        suffix =  re.search(pattern+".(\d+).dtas",dta).group(1)
        suffixes.append(int(suffix))
    foldersuffix = max(suffixes)

    # foldersuffix =  re.search(pattern+".(\d+).dtas",dtas).group(1)
    # foldersuffix = len(pepxml_all)
    # print (foldersuffix)
    pepxlm_moved = pepxml_new[0].split(".pep.xml")[0]+"."+str(foldersuffix)+".pep.xml"
    # print (pepxlm_moved)
    newfolder = mzFol2.split("/")[-1]
    createFolder = mzFol2+"/"+newfolder+"."+str(foldersuffix)
    #cmd = "mkdir "+createFolder
    #os.system(cmd)
    log_folder = mzFol2+"/log"
    paramsFile = mzFol2+"/comet.params"
    #cmd2 = "mv "+log_folder+" "+createFolder

    #os.system(cmd2)
    cmd3 = "cp "+paramsFile+" "+createFolder
    # os.system(cmd3)

    cmd4 = "mv "+pepxml_new[0]+" "+pepxlm_moved
    os.system(cmd4)

    #This part is added to change the nomenclature of the outfile so that we have all JUMP-f consider all precursor ions
    reqd_dta = glob.glob(mzFol2+"/"+basefile+"."+str(foldersuffix)+".dtas")[0]
    
    outFileDict, outFileDictPrec = dtaIonDictMap(reqd_dta)
    pepXMLOutFileConversion(pepxlm_moved, basefile, outFileDict, outFileDictPrec)
   
    cometOrigPep = glob.glob(mzFol2+"/*.cometOriginal")[0]
    
    #push the backup comet original pepxml file to log folder for backup
    cmd5  = "mv "+cometOrigPep+" "+log_folder
    os.system(cmd5)





def fromListToFolder(filelist):
    folder = []
    for mzFile in filelist:
        dirToMake = mzXML_path+"/"+os.path.basename(mzFile).split(".")[0]
        folder.append(dirToMake)
        #os.system("mkdir "+dirToMake)
        #os.system("mv "+mzFile+" "+dirToMake)
        os.system("cp comet.params "+dirToMake+"/comet.params")
    return folder



#dtas = glob.glob(work_path+"/*/*.1.dtas")[0]
#new_ms2 = dtas.split(".1.dtas")[0]+".ms2"

#def isclose(a, b, rel_tol=1e-09, abs_tol=0.0): #This is checking the closeness of the precurso ions m/z
 #   return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

def isclose(a, b, tol=1):
    a = float(a)
    b = float(b)
  #   massError = abs(a-b)*1e6/a
    massError = (b-a)*1e6/a  #calculates ppm error of the a(theoretical) from b(observed)
    if abs(massError) < tol:
        return True
    else:
        return False

def dta_to_ms2(dtas, new_ms2):
    f = open(dtas,"r")
    line=f.readlines()
    dtas_dict = {}
    count_outfile = 0
    for x in range(0, len(line),3):
        dta = line[x]
        mass_ms2 = line[x+1]
        ms2_int = line[x+2]
        dta_info = dta.split(".")
        file = dta_info[0]
        scan = dta_info[1]
        ppi = dta_info[2]
        charge = dta_info[3]
       

        neutral_mass = dta.split()[-2] #[M+H]+1

       
        check = dta_info[2]
        prec_mz = str(precMZCalc(neutral_mass, charge))
        
        count_outfile+=1
        if int(scan) not in dtas_dict:
            dtas_dict[int(scan)] = [[dta,mass_ms2,ms2_int,dta_info,file,scan, charge,neutral_mass, prec_mz]]
        else:
            dtas_dict[int(scan)].append([dta,mass_ms2,ms2_int,dta_info,file,scan, charge,neutral_mass, prec_mz])
    now = datetime.now()
    print("now =", now)
    # dd/mm/YY H:M:S
    dt_string = now.strftime("%m/%d/%Y %H:%M %p")
    year = now.strftime("%Y")
    date_time = dt_string.split() #date is 0 index, time is 1 index, AM/PM is 2 index

    header_ms2 = "H\tCreationDate\t"+dt_string+"\nH\tExtractor\tMakeMS2\nH\tExtractorVersion\t1.0\nH\tComments\tMakeMS2 written by Suresh Poudel, "+year+"\nH\tExtractorOptions\tMS2/MS1\n"

    #new_ms2 = mzxml_file.split(".mzXML")[0]+".ms2"
    with open(new_ms2,"w") as new_ms2:
        new_ms2.write(header_ms2)
        od = collections.OrderedDict(sorted(dtas_dict.items()))
        count_key = 0
        for key, all_values in od.items():
            count_key+=1
#             if len(all_values) > 1:
#                 print ("\nThis scan has multiple charges, the scan = ", key)
            for mul_values in all_values:
                new_ms2.write("S\t"+str(key)+"\t"+str(key)+"\t"+mul_values[-1]+"\nZ\t"+mul_values[-3]+"\t"+mul_values[-2]+"\n")
                all_mass_list = mul_values[1].split()
                all_int_list = mul_values[2].split()
                for index, val in enumerate(all_mass_list):
                    new_ms2.write(val+"\t"+all_int_list[index]+"\n")
              
    print ("\nTotal dta keys = ",count_key)
#jump -deisotope jump_ss_HH_tmt10_mouse.params mix_ratio.mzXML

def dtaIonDictMap(dtas):
    f = open(dtas,"r")
    line=f.readlines()
    outFileDict = {} #This is the outfile dictionary that have the ppi information 
    outFileDictPrec = {} #This stores the precursor mass to get the accuate ppi information
    dtas_dict = {}
    count_outfile = 0
    for x in range(0, len(line),3):
        dta = line[x]
        mass_ms2 = line[x+1]
        ms2_int = line[x+2]
        dta_info = dta.split(".")
        file = dta_info[0]
        scan = dta_info[1]
        ppi = dta_info[2]
        charge = dta_info[3]
        scan5digit = str('%05d' % int(scan)) #convert to 5 digit for example 1000 scan will be 01000
        scanKey = file+"."+scan5digit+"."+scan5digit+"."+charge #resembles the outfile format of comet pep.xml file 
        
        #if int(charge) != 1: #comet does not search +1 charge for .ms2 so excluding that
        if scanKey not in outFileDict.keys(): #this collects ppi in same order as dtas so these can be used  
            outFileDict[scanKey] = [ppi]
        else:
            outFileDict[scanKey].append(ppi)

        neutral_mass = dta.split()[-2] #[M+H]+1

        if scanKey not in outFileDictPrec.keys(): #this collects ppi in same order as dtas so these can be used
            outFileDictPrec[scanKey] = [neutral_mass]
        else:
            outFileDictPrec[scanKey].append(neutral_mass)

     
    return outFileDict, outFileDictPrec

def create_job_file( reqd_dta,mzxml_file, sample, comet_params):
    #fileroot = mzxml_file.split("/")[-1].split(".ms2")[0]
    #dta_to_ms2(reqd_dta, mzxml_file)
    log_dir = mzXML_path+"/"+sample+"/log"
    cmd1 = "rm -r "+log_dir
    try:
        os.system(cmd1)
    except:
        print ("No log directory for ", sample)
    os.system("mkdir "+log_dir)
    job_header = "#!/bin/bash\n#BSUB -P TestComet\n#BSUB -J comet\n#BSUB -oo "+log_dir+"/log.out\n#BSUB -eo "+log_dir+"/error.err\n#BSUB -n 8\n"
    #activate this if you need 2018 comet
    #cnt_pepxml = len(glob.glob(mzXML_path+"/"+fileroot+"/"+fileroot+"*.pep.xml"))
    #-Ntest.3
    #suffix = "."+str(cnt_pepxml+1)
    dta_to_ms2(reqd_dta, mzxml_file)
    job_body1 = "comet -P"+comet_params+" "+mzxml_file

    #job_body1 = comet+" -P"+comet_params+" "+mzxml_file

    jobfile = sample+".sh"
    with open(jobfile,"w") as new_file:
        new_file.write(job_header+job_body1)
    return jobfile

def submit_job(jobf,queue,mem):
  cmd = 'bsub -q '+queue+' -R "rusage[mem='+mem+']" < '+jobf
  os.system(cmd)

def pepXMLOutFileConversion(pepxmlFile,basefile, outFileDict, outFileDictPrec):
    suffix = basefile
    with fileinput.FileInput(pepxmlFile, inplace=True, backup='.cometOriginal') as f:

        for line in f:
            if "<spectrum_query spectrum=" in line:
                pattern = '<spectrum_query spectrum="('+suffix+'\.\d+\.\d+\.\d+)".+precursor_neutral_mass="(\d+\.\d+)"'

                allPat=re.match(pattern, line.strip())
                x = allPat.group(1)
                neutralMass = allPat.group(2) 
                MH_mass = float(neutralMass) + proton 
                outfileSplit = x.split(".")

                for val in range(0,len(outFileDictPrec[x])):
                    if isclose(MH_mass, float(outFileDictPrec[x][val])):
                        break
    #             if x == "HL_human.07462.07462.3":
    #                 print (val,"\t",outFileDict[x][val])
    #                 print (len(outFileDictPrec[x]))
    #                 print (outFileDictPrec[x])
                add_value = outFileDict[x][val]
                y = outfileSplit[0]+"."+outfileSplit[1]+"."+add_value+"."+outfileSplit[-1]
                print (line.rstrip().replace(x,y))
            else:
                print (line.rstrip())

def mzXMLtoMS2(mzXML):
    sample = mzXML.split("/")[-1].split(".mzXML")[0]
    dtas = glob.glob(mzXML_path+"/"+sample+"/"+sample+".*dtas")
    pattern =  mzXML_path+"/"+sample+"/"+sample
    suffixes = []
    for dta in dtas:
        suffix =  re.search(pattern+".(\d+).dtas",dta).group(1)
        suffixes.append(int(suffix))
    value = max(suffixes)

    reqd_dta = mzXML_path+"/"+sample+"/"+sample+"."+str(value)+".dtas"

    new_ms2 = reqd_dta.split("."+str(value)+".dtas")[0]+".ms2"
    jobfile = create_job_file(reqd_dta, new_ms2,sample, cometFlyParams)
    submit_job(jobfile,queue,mem)
    

cmd = "jump -deisotope "+jump_params+" "+" ".join(mzXMLs)
os.system(cmd)

# print ("################################################################\n")
# print ("#       ****  Starting COMET search                 ****       #\n")
# print ("################################################################\n\n")
print ("        ****  Starting COMET search                 ****        \n\n")

#cometParams = "comet_HH_tmt10_mouse.params"
if (mzXMLs == ["*.mzXML"]) or (glob.glob(mzXML_path+"/*.mzXML")==[]):
    mzXMLs = glob.glob(mzXML_path+"/*/*.mzXML")

if len(mzXMLs) < 40:
    for mzXML in mzXMLs:
        
        mzXMLtoMS2(mzXML)


        print ("\n\nJob is submitted for COMET search on "+mzXML+"\n\nPLEASE WAIT PATIENTLY\n\n")
else:
    rounds = math.ceil(len(mzXMLs)/40.0)
    a=0
    for total_jobs in range(1,rounds+1):
        total_work = total_jobs*40
        if total_work < len(mzXMLs):
            for mzXML in mzXMLs[a:total_work]:
                mzXMLtoMS2(mzXML)
 

                print ("\n\nJob is submitted for COMET search on "+mzXML+"\n\nPLEASE WAIT PATIENTLY\n\n")
        #while total_work >= len(mzXMLs[1:]):
            new_log = glob.glob(mzXML_path+"/*/log/log.out")
            #print ("Total log files = ",new_log)
            hold=0
            while hold!=1:
                
                if len(new_log) < total_work:
                    time.sleep(60)
                    new_log = glob.glob(mzXML_path+"/*/log/log.out")
                    for logF in new_log:
                        print (dirname(dirname(logF))+".mzXML search is complete\n")        
                    #print ("Total log files = ",len(new_log))
                else:
                    a = total_work
                    hold = 1
        else:
            for mzXML in mzXMLs[a:len(mzXMLs)]:
                mzXMLtoMS2(mzXML)
 
new_folders = fromListToFolder(mzXMLs)
# print (new_folders)
#print ("Total log files = ",new_log)

finish_log = glob.glob(mzXML_path+"/*/log/log.out")
while len(finish_log) != len(mzXMLs):
    finish_log = glob.glob(mzXML_path+"/*/log/log.out")
    time.sleep(30)


for mzFol2 in new_folders:
  filenameMS2 = mzFol2.split("/")[-1]
  mvLogsParams(mzFol2, filenameMS2)



print ("CONGRATULATIONS !!! All searches have been completed.\n\n")
#cmdCp = "cp "+comet_params+" "+mzXML_path+"/comet.params"
#os.system(cmdCp)







