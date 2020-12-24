import glob
import subprocess
import fileinput
import re
import os
import sys
import configparser

mzxml_path = sys.argv[1]

#comet = "/home/spoudel1/JUMP_CometParams/JUMPparams/comet"
#path = os.getcwd()
subprocess.call(["comet","-p"])

def find_replace(line_replace, pat1_list, pat2_list):
    for index, value in enumerate(pat1_list):
        if value in line_replace:
            line_replace =  re.sub(pat1_list[index],pat2_list[index],line_replace.rstrip())
        else:
            line_replace = line_replace.strip()
    return line_replace
    

def modify_create_new_params(file1, file2, pat1, pat2):
    with open(file1,"r") as f, open(file2,"w") as g:
        for line in f:
            if "commonly adjusted parameters" in line:
              next_line = f.readline()
              if "search_engine = comet" in next_line:
                g.write(line+"\n"+next_line+"\n")
              else:
                line = line.strip()+"\n"+"search_engine = comet"
                g.write(line+"\n"+next_line+"\n")
            else:
              line = find_replace(line, pat1, pat2)
              g.write(line+"\n")
    return file2

#for TMThh 
#copy jump-s parameter files as jump_sc_HH_tmt10_human.params
dataTypes = ["TMTpro","TMTpro_pho", "HH", "HL", "TMThh","TMThhpho"]

for data in dataTypes:
    paramsFile = glob.glob("ParameterFiles/"+data+"/jump_sj_*_human.params")[0]
    paramsSplit = paramsFile.split("_sj_")
    jumpcometParams = paramsSplit[0]+"_sc_"+paramsSplit[1]
    cmdcopy = "cp "+paramsFile+" "+jumpcometParams

    os.system(cmdcopy) 

    with fileinput.FileInput(jumpcometParams, inplace=True, backup='.bak') as file:
        for line in file:
            print(line.replace("search_engine = JUMP", "search_engine = COMET"), end='')




pattern1TMThh = ["database_name = /some/path/db.fasta","num_threads = 0","output_suffix =","add_Nterm_peptide = 0.0","add_K_lysine = 0.0000"]
pattern2TMThh = ["database_name = /hpcf/authorized_apps/proteomics_apps/database/20200422/HUMAN/human_ft_mc2_c0_TMT_K229.fasta","num_threads = 4","output_suffix =","add_Nterm_peptide = 229.162931","add_K_lysine = 229.162932"]

comet_params = "comet.params.new"
comet_params_newTMThh = "ParameterFiles/TMThh/comet_HH_tmt10_human.params"
#comet_source_newTMThh = "ParameterFiles/TMThh/comet.params"

modify_create_new_params(comet_params,comet_params_newTMThh,pattern1TMThh,pattern2TMThh)
#os.system("cp "+comet_params_newTMThh+" "+comet_source_newTMThh)
#for HH and HL

pattern1hh = ["database_name = /some/path/db.fasta","num_threads = 0","output_suffix ="]
pattern2hh = ["database_name = /hpcf/authorized_apps/proteomics_apps/database/20200422/HUMAN/human_ft_mc2_c0_TMT_K229.fasta","num_threads = 4","output_suffix = .1"]


comet_params_newHH = "ParameterFiles/HH/comet_HH_human.params"
comet_params_newHL = "ParameterFiles/HL/comet_HL_human.params"
#comet_source_newHH = "ParameterFiles/HH/comet.params"
#comet_source_newHL = "ParameterFiles/HL/comet.params"

modify_create_new_params(comet_params,comet_params_newHH,pattern1hh,pattern2hh)
modify_create_new_params(comet_params,comet_params_newHL,pattern1hh,pattern2hh)

#os.system("cp "+comet_params_newHH+" "+comet_source_newHH)
#os.system("cp "+comet_params_newHL+" "+comet_source_newHL)

#for TMThhpho

pattern1TMThhpho = ["variable_mod02 = 0.0 X 0 3 -1 0 0","use_NL_ions = 0"]
pattern2TMThhpho = ["variable_mod02 = 79.966331 STY 0 3 -1 0 0","use_NL_ions = 1"]


comet_params_newTMTpho = "ParameterFiles/TMThhpho/comet_HH_pho_tmt10_human.params"
#comet_source_newTMTpho = "ParameterFiles/TMThhpho/comet.params"

modify_create_new_params(comet_params_newTMThh,comet_params_newTMTpho,pattern1TMThhpho,pattern2TMThhpho)
#os.system("cp "+comet_params_newTMTpho+" "+comet_source_newTMTpho)

subprocess.call(["rm",comet_params])

jumpf_params_source = "ParameterFiles/TMThh/jump_fj_HH_tmt10_human.params"
jumpf_paramsTMThhComet = "ParameterFiles/TMThh/jump_fc_HH_tmt10_human.params"


filetypes = ('*.mzXML','*.raw','*.RAW')
files_grabbed = []
for files in filetypes:
  files_grabbed.extend(glob.glob(files))

#mzXMLs_all = []
#for all_files in files_grabbed:
#  if ".mzXML" not in all_files:
#    mzXMLs_all.append(all_files.split(".")[0]+".mzXML") 
#  else:
#    mzXMLs_all.append(all_files)   

mzXMLs = list(set(files_grabbed))  #mzXML files are expected in the same place as the script = "test1_Comet:/home/spoudel1/PengProteomics/TestComet/FTLD_f42.1" #you can have any name for the results directory here, I use an example as test1

if len(mzXMLs) >= 1:
  input_mzXMLs_line = "test1_Comet:"
  for files in mzXMLs:
    path2add = files.split(".")[0]+"/"+files.split(".")[0].split("/")[-1]+".1"
    input_mzXMLs_line+=path2add+"\n"
else:
  input_mzXMLs_line = "test1_Comet:/home/spoudel1/PengProteomics/TestComet/FTLD_f42.1"

print (input_mzXMLs_line)
pattern1_jumpf = ["HH_tmt10_human_jump:/hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.0/SampleData/TMThh/HH_tmt10_human_jump/HH_tmt10_human_jump.1","pit_file = 0","min_XCorr = 10","mass_accuracy = 1","output_pepXML = 1"]  #currently present in jump_fj_params
pattern2_jumpfTMThhComet = [input_mzXMLs_line,"pit_file = /hpcf/authorized_apps/proteomics_apps/database/20200422/HUMAN/human_ft_mc2_c0_TMT_K229.pit","min_XCorr = 1","mass_accuracy = 0","output_pepXML = 0"]  #values to be modified


modify_create_new_params(jumpf_params_source,jumpf_paramsTMThhComet,pattern1_jumpf,pattern2_jumpfTMThhComet)



jumpf_paramsTMThhphoComet = "ParameterFiles/TMThhpho/jump_fc_HH_pho_tmt10_human.params"


pattern1_TMThhComet_source = ["unique_protein_or_peptide = protein","initial_outfile_fdr = 5","one_hit_wonders_removal = 2","mods = 0",
            "min_outfile_num_for_XCorr_filter = 500","one_hit_wonders_min_dCn = 0"]
pattern1_TMThhphoComet_new = ["unique_protein_or_peptide = peptide","initial_outfile_fdr = 10","one_hit_wonders_removal = 0","mods = STY",
            "min_outfile_num_for_XCorr_filter = 200","one_hit_wonders_min_dCn = 0.1"]

modify_create_new_params(jumpf_paramsTMThhComet,jumpf_paramsTMThhphoComet,pattern1_TMThhComet_source,pattern1_TMThhphoComet_new)



jumpf_paramshhComet = "ParameterFiles/HH/jump_fc_HH_human.params"
jumpf_paramshlComet = "ParameterFiles/HL/jump_fc_HL_human.params"

pattern1_TMThhComet_source2 = ["one_hit_wonders_removal = 2","one_hit_wonders_min_dCn = 0"]
pattern1_HHComet_new = ["one_hit_wonders_removal = 0","one_hit_wonders_min_dCn = 0.1"]

modify_create_new_params(jumpf_paramsTMThhComet,jumpf_paramshhComet,pattern1_TMThhComet_source2,pattern1_HHComet_new)

modify_create_new_params(jumpf_paramsTMThhComet,jumpf_paramshlComet,pattern1_TMThhComet_source2,pattern1_HHComet_new)

idtxt_path = mzxml_path+"/sum_test1_Comet"
idtxt = idtxt_path+"/ID.txt"

pattern1_jumpq = ["idtxt = /hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.0/SampleData/TMThh/sum_HH_tmt10_human_jump/ID.txt"
                  ,"save_dir = HH_tmt10_human_sequest"
                  ,"impurity_matrix = /hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.0/JUMPq/TMT10.ini"
                  ,"tmt_reporters_used = sig126; sig127N; sig127C; sig128N; sig128C; sig129N; sig129C; sig130N; sig130C; sig131","sig131 = S10"]
pattern2_jumpq = ["idtxt = "+idtxt,"save_dir = test1_Comet","impurity_matrix = /hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.0/JUMPq/TMT11.ini"
                  ,"tmt_reporters_used = sig126; sig127N; sig127C; sig128N; sig128C; sig129N; sig129C; sig130N; sig130C; sig131N; sig131C",
                  "sig131N = S10\nsig131C = S11"]


# In[35]:
jumpq_params_source = "ParameterFiles/TMThh/jump_qs_HH_tmt10_human.params"


#This creates a new parameter file for jump quantification with required modifications

jumpq_params_tmthh = "ParameterFiles/TMThh/jump_qc_HH_tmt10_human.params"

modify_create_new_params(jumpq_params_source,jumpq_params_tmthh,pattern1_jumpq,pattern2_jumpq)

jumpq_params_tmthhpho = "ParameterFiles/TMThhpho/jump_qc_HH_pho_tmt10_human.params"

IDmod = idtxt_path+"/IDmod.txt"


pattern1_jumpqTMT = ["idtxt = /hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.0/SampleData/TMThhpho/sum_HH_pho_tmt10_human_jump_mod/IDmod.txt","save_dir = HH_tmt10_human_sequest","impurity_matrix = /hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.0/JUMPq/TMT10.ini",
                  "min_intensity_method_1_2_psm = 1, 4","min_intensity_value_1_2_psm = 2000, 10000","percentage_trimmed = 25"]
pattern2_jumpqTMTpho = ["idtxt = "+IDmod,"save_dir = test1_Comet_mod","impurity_matrix = /hpcf/authorized_apps/proteomics_apps/pipeline/release/version1.13.0/JUMPq/TMT11.ini","min_intensity_method_1_2_psm = 0","min_intensity_value_1_2_psm = 0","percentage_trimmed = 50"]


modify_create_new_params(jumpq_params_source,jumpq_params_tmthhpho,pattern1_jumpqTMT,pattern2_jumpqTMTpho)




#make params file for stagewise PTM analysis 

config = configparser.ConfigParser()

comment1 = "# Stage of round of PTM. For example, if the JUMP-s has been done, the stage is round 1. If first round of PTM is completed, stage is round 2 and so on. Mostly in our case, stage 1 is phospho PTM and stage 2 Ubiquitin. This changes for phosphoenriched data. If the JUMP-s and QC is done on phosphoenriched data search for STY, then round 1 will be another PTM such as Ubiquitin\n"
comment2 = "#this is very important step as the mzfile_path changes with stage1 and stage 2. For stage 1, it is always mzXML file as .ms2 is first created at this stage\n"
comment3 = "#please add the desire directory name where you want to write your .ms2 files\n"
comment4 = "#if this is stage 1, there is no idtxt_ref at this point, this is after later stage where complete path for ID.txt file that you want to work with example /research/projects/penggrp/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples/PTM_paper_2020/Final_Results/ms2_round1/sum_result_FDR1_5_on_mod_phopho_search_mod/IDmod.txt\n"
comment5 = "#if this is stage 1, you need to use the QC folder that comes after the using JUMP-qc /research/projects/penggrp/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples/PTM_paper_2020/Final_Results/Pho_2020/comprehensive_search/qc_qcFirstRound_Pho_Final .....  if you are in other stages, you can simply keep the value as 0\n"
comment6 = "#if this is stage 1, you can keep this as 0 but it this is stage 2,3,4 ... this step is crucial. With comet search it assigns the symbol *,#,@ to the amino acid requires modification. It depends on which variable mod number you assign, if we keep M +15.99 as first, we get * but we may get # if kept second and so on, so we need to give what mod was done on previous round for example it could be K'#'and this key is used to filter for next ms2 generationi. NOTE since * can be problem replace write '0' if you have * in mod\n" 
comment7 = "#this is the id_uni_pep.txt file containing mod of the round for the development of rules for the appropriate cutoffs to get good matches: for example /research/rgs01/project_space/penggrp/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples/ms2_round1/sum_Peptide_GroupSize1000_FDR1_5_on_mod_phopho_search_mod/publications/id_uni_pep.txt\n"
comment8 = "#this is the file from whole proteome search id_uni_pep.txt so that we can compute the psm ratio of peptides containing MOD with respect to psm number of proteins that this peptide maps to from unmodified search\n" 
comment9 = "#this is the file from enriched phosphorylation study. This file is used as the reference to see the overlap between the whole proteome phosphorylation study. The presence absence column in the output file are based on whether the peptide is present in the enriched data or not\n"


config['makeMS2'] = {comment1+"stage": "1\n",
                     comment2+"mzfile_path":"/research/rgs01/project_space/penggrp/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples\n",
                    comment3+"out_dir":"psm_rec99_ms2\n",
                    comment4+"idtxt_ref":"0\n",
                    comment5+"qc_path":"/path/containing/qc/results\n",
                    comment6+"prev_mod":"#\n",
                    comment7+"id_mod_file":"/research/rgs01/project_space/penggrp/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples/ms2_round1/sum_Peptide_FDR1_5_on_mod_phopho_search_mod/publications/id_uni_pep.txt\n",
                    comment8+"whl_unmod_id":"/research/rgs01/project_space/penggrp/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples/sum_YuxinTest_pooled_whole3/publications/id_uni_pep.txt",
                    comment9+"id_ref_pho_file":"/research/rgs01/project_space/penggrp/Alzheimer_BannerSunInstitute/penggrp/proteomics/batch0_pooledSamples/PTM_paper_2020/Final_Results/Pho_2020/comprehensive_search/sum_human_pho1_5FDR_mod/publications/id_uni_pep.txt"}

with open('ParameterFiles/jump_PT_ms2.params', 'w') as configfile:
    config.write(configfile)


