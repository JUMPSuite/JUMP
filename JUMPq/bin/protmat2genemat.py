# In[1]:

import pandas as pd
import glob
import numpy as np
import argparse
# In[2]:


pd.options.mode.chained_assignment = None

def msg(name=None):
    return '''\n\n\npython proteogenomics_software_version3.py --input_file\n\n\n '''

parser = argparse.ArgumentParser(description="Gene Map to Protein", prog="proteogenomics_software_version3.py",usage=msg())
parser.add_argument("inputfile", help="tab delimited output file from jump q")
#parser.add_argument("mzXML",help="single or list of mzXML files",nargs='+')
args = parser.parse_args()

in1 = args.inputfile

if "_prot.txt" in in1:
    ou1 = in1.split("_prot.txt")[0]+"_gene.txt"
if "_prot_quan.txt" in in1:
    ou1 = in1.split("_prot_quan.txt")[0]+"_gene_quan.txt"

if "_pep.txt" in in1:
    ou1 = in1.split("_pep.txt")[0]+"_pep_gene.txt"
if "_pep_quan.txt" in in1:
    ou1 = in1.split("_pep_quan.txt")[0]+"_pep_gene_quan.txt"


def clean_input_file(infile,outfile):
    a=0
    with open(infile,"r") as f, open(outfile,"w") as g:
        for line in f:
            while a==0:
                if ("#" in line.strip()) & ("\t" in line.strip()):
                    a=1
                else:
                    line = f.readline()
            g.write(line.strip()+"\n")

in1_clean = in1.split(".txt")[0]+"_clean.txt"

clean_input_file(in1,in1_clean)


# in1 = "/mnt/c/Users/spoudel1/Desktop/jump-f-batch/Gene_Quantification_12102019/combined_norm_uni_prot.txt"
# ou1 = "/mnt/c/Users/spoudel1/Desktop/jump-f-batch/Gene_Quantification_12102019/combined_norm_uni_gene.txt"


# In[7]:


def average_PSM(row, df):
    header = list(df.columns)
    select_PSMs=[]
    for names in header:
        if "PSM" in names:
            select_PSMs.append(names)
    series = row[select_PSMs].astype('float')
    avg_psm = series.sum(skipna=True)/len(select_PSMs)
    return avg_psm

# In[8]:
def average_intensity(row, df):
    header = list(df.columns)
    select_channel=[]
    for names in header:
        if "sig" in names:
            select_channel.append(names)
#     avg_intensity = row[select_channel].sum(skipna=True)/len(select_channel)
    series = row[select_channel].astype('float')
    avg_intensity = series.sum(skipna=True)/len(select_channel)
    return avg_intensity

# In[9]:

'''This function looks for the priority sp| signature in the uniprot dataset and generate a unique gene to protein map based on
following criteria
1. if the mapping is to only one protein -- that is automatically taken
2. if the mapping is to more than one protein 
    a. it first consiers the canonical forms
    b. if there are multiple canonical forms matches it looks for psm
    c. if there are same psm matches, it looks for average intensity and chooses the protein that has maximum average intensity
    d. similarly it looks at isoforms and chooses the isoform with lowest number
    e. if multiple isoforms arise with same low number arises it does step 2b, 2c for isoforms too'''

def sp_only(select_val_sp, psm_list,intensity_list,sp_index):
    final_value = ""
    if len(select_val_sp) == 1:
        final_value = select_val_sp[0]
    else:
        canonical = []
        canonical_index = []
        isoforms = []
        iso_index = []

        for i,new_val in enumerate(select_val_sp):
            if "-" not in new_val:
                canonical.append(new_val)
#                 canonical_index.append(select_val_sp.index(new_val))
                canonical_index.append(sp_index[i])
            else:
                isoforms.append(new_val)
#                 iso_index.append(select_val_sp.index(new_val))
                iso_index.append(sp_index[i])
        if len(canonical) == 1:
            final_value = canonical[0]
        if len(canonical) > 1:
            new_psm_list_canonical = [psm_list[index] for index in canonical_index]
            max_psm = max(new_psm_list_canonical)
            cnt_psm = new_psm_list_canonical.count(max_psm)
            
            if cnt_psm == 1:
                max_psm_index = new_psm_list_canonical.index(max_psm)
                final_value = select_val_sp[max_psm_index]
            else:
                new_intensity_list_canonical = [intensity_list[index] for index in canonical_index]
                max_intensity_index = new_intensity_list_canonical.index(np.max(new_intensity_list_canonical))
                final_value = select_val_sp[max_intensity_index]
                
        if (len(canonical) == 0) & (len(isoforms)==1):
            final_value = isoforms[0]
            
        if (len(canonical) == 0) & (len(isoforms)>1):
            
            iso_number = []
            
            for isoValue in isoforms:
                
                iso_num=isoValue.split("-")[1][0]
                iso_number.append(iso_num)
                
            lowest_iso = iso_number.index(min(iso_number))
            cnt_low_iso = iso_number.count(min(iso_number))
                                        
            
            if cnt_low_iso == 1:
                final_value = isoforms[lowest_iso]
            else:
                new_psm_list_isoforms = [psm_list[index] for index in iso_index]
                max_psm_iso = np.max(new_psm_list_isoforms)
                cnt_psm_iso = new_psm_list_isoforms.count(max_psm_iso)
                if cnt_psm_iso == 1:
                    max_psm_iso_index = new_psm_list_isoforms.index(max_psm_iso)
                    final_value = select_val_sp[max_psm_iso_index]
                else:
                    new_intensity_list_isoforms = [intensity_list[index] for index in iso_index]
                    max_intensity_iso_index = new_intensity_list_isoforms.index(np.max(new_intensity_list_isoforms))
                    final_value = select_val_sp[max_intensity_iso_index] 
    #print (final_value)    
    return final_value

# In[10]:


'''This function looks for the protein signature besides swissprot sp| signature and generate a unique gene to protein map based on
following criteria
1. if the mapping is to only one protein -- that is automatically taken
2. if the mapping is to more than one protein 
    a. if there are multiple  matches it looks for psm
    b. if there are same psm matches, it looks for average intensity and chooses the protein that has maximum average intensity
    '''

def other_prot_only(select_val_tr, psm_list1, intensity_list1):
    final_value1 = ""
    if len(select_val_tr) == 1:
        final_value1 = select_val_tr[0]
    else:
        
        canonical1 = []
        canonical_index1 = []

        for i1,new_val1 in enumerate(select_val_tr):
            canonical1.append(new_val1)
            canonical_index1.append(select_val_sr.index(new_val1))
                
        if len(canonical1) > 1:
            new_psm_list_canonical1 = [psm_list1[index] for index in canonical_index1]
            max_psm1 = np.max(new_psm_list_canonical1)
            cnt_psm1 = new_psm_list_canonical1.count(max_psm1)
            
            if cnt_psm1 == 1:
                max_psm_index1 = new_psm_list_canonical1.index(max_psm1)
                final_value1 = select_val_sp[max_psm_index1]
            else:
                new_intensity_list_canonical1 = [intensity_list1[index] for index in canonical_index1]
                max_intensity_index1 = new_intensity_list_canonical1.index(np.max(new_intensity_list_canonical1))
                final_value1 = select_val_sp[max_intensity_index1]
        
    return final_value1
                


# In[16]:


def select_final_valtest(row):
    val_sp = []
    val_sp_index = []
    
    val_tr = []
    val_tr_index = []
    
    val_co = []
    val_co_index = []
    
    val_cu = []
    val_cu_index = []
    
    final_value_main = ""
    
    intensity_list_avg = list(row["average_Intensity"])
    psm_list_avg = list(row["average_PSMs"])
    access_list =list(row["Protein Accession #"])
    
    for ind,val in enumerate(access_list):
        if "sp|" in val:
            val_sp.append(val)
            val_sp_index.append(ind)
            
        if "tr|" in val:
            val_tr.append(val)
            val_tr_index.append(ind)
            
        if "co|" in val:
            val_co.append(val)
            val_co_index.append(ind)
            
        if "cu|" in val:
            val_cu.append(val)
            val_cu_index.append(ind)
        
        if len(val_sp) != 0:
            final_value_main = sp_only(val_sp, psm_list_avg, intensity_list_avg, val_sp_index)
        
        elif len(val_tr)!=0:
            final_value_main = sp_only(val_tr, psm_list_avg, intensity_list_avg, val_tr_index)
        
        elif len(val_co)!=0:
            final_value_main = sp_only(val_co, psm_list_avg, intensity_list_avg, val_co_index)
        
        else:
            final_value_main = sp_only(val_cu, psm_list_avg, intensity_list_avg, val_cu_index)
    print (row["GN"], "\t------\t", final_value_main)
    
    return final_value_main


# In[ ]:


def merge_df(df1,df2):
    df3 = pd.merge(df1,df2,how="inner",left_on="accession",right_on="Protein Accession #")
    return df3


# In[74]:
print ("\n\nThis program will map one gene with one protein.The best matches are printed in the screen\n\n\nGene Name\t------\tSelected Protein Accession ------------- (It takes some time to print for very large files)")

df_quant = pd.read_csv(in1_clean, delimiter = "\t", low_memory=False)
df_quant_clean=df_quant.dropna(subset=["GN"])

df_quant_clean["average_PSMs"]=df_quant_clean.apply(average_PSM, df=df_quant_clean, axis=1)
df_quant_clean["average_Intensity"]=df_quant_clean.apply(average_intensity, df=df_quant_clean, axis=1)


# In[75]:


df_quant_groupGN=df_quant_clean.groupby(["GN"])[["Protein Accession #","average_PSMs","average_Intensity"]].agg(lambda x: x.tolist()).reset_index()


# In[76]:


df_quant_groupGN["accession"]=df_quant_groupGN.apply(select_final_valtest, axis=1)


# In[77]:


merge_gene_map = merge_df(df_quant_groupGN,df_quant_clean)


# In[78]:

all_columns = list(merge_gene_map.columns)

extract_col = [all_columns[0],all_columns[4],all_columns[5]]+all_columns[7:-2]

headers = []
for x in extract_col:
    if ("_" in x) or (".1" in x):
        headers.append(x[0:-2])
    else:
        headers.append(x)

for pos, head in enumerate(headers):
    if head == "accession":
        headers[pos] = "Protein Accession #"



final_result=merge_gene_map[extract_col]


# In[81]:


final_result.to_csv(ou1,sep="\t",index=None,header=headers)

