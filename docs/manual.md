# 1. Introduction
JUMP SUITE is a software integrating database creation, database search, identification filtering and protein quantification. 
The system starts with a raw file and ends with a list protein of interest. 
It is capable of identifying multiple candidate peptides from mixture spectra and producing de novo sequence tags and generates a Jscore that merges the local tag score and the global pattern matching score. 
JUMP is a sensitive and specific database search method for peptide identification and outperforms than other search engines such as SEQUEST, Mascot, InsPecT, and PEAKS DB.
Currently, JUMP SYSTEM has four main components: database creation, database search, identification filtering, and protein quantification.

## 1.1. JUMP database creation
JUMP database creation is a tool designed to build a database for JUMP search and a corresponding protein inference table (PIT) which is necessary for the protein grouping when running jump filtering. 
It requires a parameter file, `jump_d.params`, in which some parameters for the database generation and a jump search parameter file (e.g. `jump_s.params`) should be specified.

## 1.2 JUMP search
JUMP search is a new hybrid database search algorithm that combines database search and de novo sequencing for protein identification by tandem mass spectrometry. 
Unlike other existing search algorithms, JUMP generates amino acid tags and ranks peptide spectrum matches by the tags as well as pattern matching, and is able to use all potential sequence tags, as short as only one amino acid, in database search. 
JUMP-derived tags allow partial de novo sequencing and facilitate the unambiguous assignment of modified residues. 
This search engine also provides additional features, such as identification of co-eluted peptides from mixture MS/MS spectra, and assignment of modification sites by tags.   

## 1.3 JUMP filtering
JUMP filtering is a program that filters peptide identifications after database search. 
It utilizes either a heuristic grouping method or a machine learning Linear Discriminant Analysis (LDA) method to discriminate correct and incorrect peptide-spectrum matches (PSMs). 
The grouping method performs FDR filtering by grouping PSMs by peptide length, trypticity, miscleavage, protein modifications, ion charge, and matching scores from JUMP search (Jscore). 
For the LDA method, it performs FDR filtering by calculating an LDA score based on search score, peptide length, Jscore and ppm.

## 1.4 JUMP quantification
JUMP quantification performs quantification for Tandem Mass Tag (TMT) tagging experiments. 
The program utilizes relative intensity of reporter ions to determine the relative quantification for each sample. 
The program also performs summarization and differential analysis between groups of samples. 

## 1.5 License
The current version of JUMP is optimized for the analysis of high resolution MS/MS spectra and is written in perl for high performance parallel computing systems but also can be run on a single machine. 
The source code of JUMP can be downloaded from Dr. Peng’s lab website (http://www.stjuderesearch.org/site/lab/peng) at St. Jude Children’s Research Hospital and used according to the terms of the GNU General Public License (GPL). 

# 2.	How to run JUMP software suite
## 2.1 Installation 
The installation procedure is well described in README.md document. JUMP software suite works under Unix-like operating system such as Linux. Windows 10 users can install and use JUMP software suite using Windows Subsystem for Linux (WSL). For the details about WSL, please refer [https://docs.microsoft.com/en-us/windows/wsl/install-win10](https://docs.microsoft.com/en-us/windows/wsl/install-win10). 

## 2.2 Preparation of mzXML file(s) 
JUMP software suite takes as input mzXML file(s). Users can convert raw data files from MS instrument to mzXML files using software such as ReAdW, Proteowizard and so on.

## 2.3 Run JUMP software suite
### 2.3.1. Run JUMP database creation
To run the JUMP program, first set up the parameter in `jump_d.params` file first (users can change the name of parameter file). For parameters setup, please go to the section of parameter setup.  To create the database, type the following command in a terminal, 
```
jump -d jump_d.params
```
 
### 2.3.2 Run JUMP search 
To perform the search function, first set up the parameter in your parameter file (for parameters setup, please go to the section of parameter setup), then type the search command `jump -s` followed by a space, then the name of the parameter file (e.g. `jump_sj_HH_tmt10_human.params`), and a space, and finally, the name of the .mzXML file(s) (e.g. `HH_tmt10_human.mzXML`). For example,
```
jump -s jump_sj_HH_tmt10_human.params HH_tmt10_human.mzXML
``` 
For the storage of results, you can either create a directory where output files will be located or just press the enter key. If you choose to create a directory, the results of jump search will be stored in the directory. If you press the enter key without creating a specifically named directory, the results will be stored in the directory which has the same name as your input file (e.g. `./HH_tmt10_human/`) and a sub-directory which has the same name as the input file with a numerical extension added to the end of the directory name (starts with 1 and then increase incrementally depending on how many times you have run the search task for the input file, like `./HH_tmt10_human/HH_tmt10_human.1/`).  
Jump search may take a while depending on the size and number of input .mzXML file(s). 
 
### 2.3.3 Run JUMP filtering 
To perform the filtering function, first set up the parameter in your parameter file (for parameters setup, please go to the section of parameter setup), then type the filtering command `jump –f` followed by the filtering parameter file name, e.g. `jump_f_HH_tmt10_human.params`, as follows. 
```
jump -f jump_f_HH_tmt10_human.params
```
All the results of the filtering task will be stored in a directory which has the prefix of "sum" (representing "summarization") and the name of output directory specified in the parameter file (e.g. `./sum_HH_tmt10_human`) under the current working directory.

### 2.3.4 Run JUMP quantification 
To run the JUMP quantification, first set up the parameter in your parameter file (for setting parameters, please see the section of parameter setup), then type the quantification command `jump –q`, followed by the quantification parameter file name `jump_q_HH_tmt10_human.params`.
```
jump -q jump_q_HH_tmt10_human.paramms
```  
JUMP quantification results will be stored in the directory with a prefix "quan" (representing "quantification"),  for example, `./quan_HH_tmt10_human/` under the working directory.  

### 2.3.5 Run JUMP localization
JUMP localization can be run by setting up the parameters in your parameter file, then type the localization command `jump –l`, followed by the localization parameter file. For example,
```
jump -l jump_l.params
```
The result will be located in the directory which has the prefix of "loc" (representing "localization") and the output directory name specified in the parameter file.

## 2.4 Parameter Setup 
A parameter file is required for each component of JUMP software suite. Sample parameter files can be obtained by simply running the command `jump -params` in your working directory. You can find a directory of `./ParameterFiles` in which sample parameter files are located. For the simple routine tasks, it is not necessary to change default parameters except path(s) of input file(s). Advanced users can change those parameters according to their needs. Whenever the parameters are re-defined/changed, please save the changes to the parameter file. Otherwise, Jump software suite will use default ones.

### 2.4.1 Setup parameters for Jump database creation
The parameters for database creation are used to create database files for JUMP search and/or a corresponding protein inference table (PIT) for the protein grouping to be performed when running JUMP filtering. All the modification information for the database depends on the `jump.params` file (please see below for more details).

1. `input_database1 = /home/user/MOUSE.fasta`  
The full (absolute) path of the input .fasta file used for creating a database. At least one .fasta file should be specified.
2. `input_database2 = /home/user/HUMAN.fasta`  
Multiple .fasta files can be used to create a database. When using multiple .fasta files, please do NOT forget numbering for the input databases (i.e. `input_database1 = ...`, `input_database2 = ...`, etc.).
3. `output_prefix = mouse`  
The prefix for newly generated database and PIT files. The suffix will be automatically generated according to the conditions for the database. For example, if a user defines `output_prefix = mouse` with the condition of fully tryptic digestion up to 2 miscleavage with no static cystein modification, then the new database will have the name such as `mouse_ft_mc2_c0.fasta.xxx`.
4. `include_contaminants = 1`  
If a user wants to include contaminant sequences to the database, this parameter needs to be set to 1. Otherwise, set to 0. If a user uses `input_database` which already include contaminant sequences, `include_contaminants` should be set to 0 to avoid the duplication of contaminant sequences in the database.
5. `input_contaminants = /data1/database/contaminats.fasta`  
The absolute path of contaminant sequences (should be .fasta format).
6. `decoy_generation = 1`  
If a user wants to include decoy sequences of proteins (and contaminants, if exist), this parameter needs to be set to 1. Otherwise, set to 0.
7. `decoy_generation_method = 1`  
If `decoy_generation = 1`, a user should specify the method of generating decoy sequences (1 = simply reversing target protein sequences, 2 = reversing protein sequences and then swapping every K and R with their preceding amino acid).
8. `jump.params = /home/user/jump_search.params`  
A user _should_ prepare a JUMP search parameter file (e.g. `jump_search.params`) in which static and dynamic modifications to be used to generate a new database. Please specify the full (absolute) path of the JUMP search parameter file.
9. `bypass_db_generation = 0`  
A user can create PIT file, but bypass the generation of database (time-consuming). If that is the case, please set the parameter to 1 (1 = bypass the database generation, 0 = no).
10. `list_protein_abundance = /home/user/human_abudance.txt`  
Optional parameter. This protein abundance information is used to generate PIT for grouping and sorting proteins. If this parameter is used, the full (absolute) path needs to be specified. The file for protein abundance information should be a tab-delimited text format with the following rules: Column1 = Uniprot accession of a protein (e.g. P12345), Column2 = abundance of the protein (numeric value), Column3, 4, 5, ... = any information/description (will be ignored).  
Multiple protein abundance information can be used as following.  
    ```
    list_protein_abundance1 = /home/user/mouse_abunance.txt
    list_protein_abundance2 = /home/user/rat_abundance.txt
    ```
11. There are some more optional parameters related with protein annotation. They will be used for PIT generation.
    ```
    list_oncogenes = /path/oncogenes_from_literatures.txt
    list_TFs = /path/tfs_from_TRANSFAC.txt
    list_kinases = /path/kinases_from_pkinfam.txt
    list_GPCRs = /path/gpcrs.txt
    list_epigenetic_factors = /path/epigenetic_regulators.txt
    list_spliceosomal_proteins = /path/spliceosomal_proteins.txt
    ```

### 2.4.2 Setup the parameters for JUMP search
1. For the database search, full (absolute) paths of database file (.mdx) and corresponding PIT file should be specified. 
    ```
    database_name = /home/user/database/HH_tmt10_human.fasta.mdx  
    pit_file = /home/user/database/HH_tmt10_human.pit
    ```
2. `peptide_tolerance = 6`  
Precursor mass tolerance (MH+). Default is 6 ppm
3. `peptide_tolerence_units = 2`  
Peptide tolerance unit: 1 = Da, and 2 = ppm.
4. `first_scan_extraction = 5000`  
It defines the first scan number for search. For the search of full scans, it should be 1 (i.e. the 1st scan).
5. `last_scan_extraction = 10000`  
It defines the last scan number for search. For the search of full scans, it should be a large number (e.g 10E6).
6. `isolation_window = 1.4`  
It is defined as +/- (isolation window)/2 based on MS2 isolation window (e.g. 1.4 m/z).
7. `isolation_window_offset = 0.25`  
It is defined as +/- isolation window offset based on MS2 isolation window (e.g. 0.25 m/z).
8. `isolation_window_variation = 0.2 `  
It is defined as +/- isolation window variation based on MS2 isolation window (e.g. 0.2 m/z).
9. `interscanppm = 15`  
It defines the mass tolerance for interscan precursor identification.
10. `intrascanppm = 10`  
It defines the mass tolerance for intrascan isotopic decharging.
11. `max_num_ppi = 0`  
It defines max precursor ions selected for mixed MS2 search (please note that "ppi" means precursor peak intensity and "num" means the number of peak series, not the peak number). If set to 0, it is disabled.
12. `percentage_ppi = 50`  
It defines the minimal percentage of precursor peak intensity (ppi) when max_num_ppi = 0 and the intensity is the sum of the intensities of all isotopic peaks.
13. `ppi_charge_0 = 1`  
If set to 0, it will discard uncharged MS1 (charge = 0). If set to 1, it will perform manual assignment of the charge (+2 and +3).
14. `ppi_charge_1 = 1`  
If set to 0, it will discard MS1 with charge +1. If set to 1, it will allow MS1 to have original charge +1.
15. `mass_correction = 2`  
0 = no mass correction, 1 = MS1-based mass correction (using the peak of 445), 2 = MS2-based mass correction (using y1 ion of K and/or R, and TMT reporter ion) 3 = manual mass correction.
16. `prec_window = 3`  
0 = disabled. If set to 1-10 (Da), it defines a m/z windows for removing precursor ions.
17. `MS2_deisotope = 1`  
0 = disable, 1 = enable to do deisotope. 
18. `ppm = 10`  
It defines the mass tolerance for MS2 decharging and deisotoping.
19. `charge12_ppm = 15`  
It defines the mass tolerance for merging different charged ions with the same mass.
20. `ms2_consolidation = 10`  
It defines a maximum number of peaks retained within each 100-Da window for processing MS2 spectra. 
21. `TMT_data = 0`  
It specifies the definition of a dataset. 0 = non-TMT data, 1 = TMT data.
22. `tag_generation = 1`  
Whether or not generate tags for peptide identification. 0 = disable, 1 = enable to generate tags.
23. `tag_tolerance = 10`  
It defeines the mass tolerance for measuring peak distance for generating tags.
24. `tag_tolerance_unit = 2`  
Tag tolerance unit: 1 = Da, 2 = ppm.
25. `tag_select_method = comb_p`  
It defines the method of ranking tags: comb_p, hyper_p or rank_p 
26. `ion_series = 1 1 0 0 0 0 0 1 0`  
It defines the type of ions. For example, 1 1 0 0 0 0 0 1 0 means a, b, c, d, v, w, x, y and z ions, respectively.
27. `frag_mass_tolerance = 15`  
It defines the mass tolerance for matching MS2 fragment ions.
28. `frag_mass_tolerance_unit = 2`  
The unit of `frag_mass_tolerance`: 1 = Da, 2 = ppm.
29. `ion_losses_MS2 = 0 0 0 0`  
It defines the type of neutral losses: 0 = disable, 1 = enable for neutral loss of H2O, HPO3, H3PO4 and NH3, respectively.
30. `ion_losses_MS1 = 0`  
0 = disabled. If set to 1, it indicates the use of precursor ion phosphate neutral loss to estimate the number of S/T phosphorylation.
31. `ion_scoring = 1`  
If set to 1, it indicates the simultaneous scoring of product ions. If set to 2, it indicates the separate scoring of ion series and the determination of charge states.
32. `matching_method = hyper_p`  
It defines the method of PSM scoring (comb_p, hyper_p, rank_p).
33. `tag_search_method = 2`  
1 = exit when found; 2 = exhaustive search using tags defined by `max_number_tag_for_search`. 
34. `max_number_tags_for_search = 50`  
It defines the maximum number of tags used for search unless the total number of tags is smaller than this defined value.
35. `number_of_selected_result = 5`  
It defines the maximum tentative number of PSMs in .spout file ranked by Jscore.
36. `number_of_detailed_result = 5`  
It defines the maximum tentative number of PSMs in .spout file ranked by pattern matching score.
37. `second_search = 1`  
0 = disable, 1 = enable. For PSMs with FDR > 0, perform the another round of search by relaxing monoisotopic mass by including M-2, M-1, M, M+1, M+2.
38. `dynamic_M = 15.99492`  
It defined dynamic modification to an amino acid. For other modifications, add each dynamic modification by one line starting with `dynamic_` (M: 15.99492; C: 57.02146 carbamidomethylation or 71.0371 acrylamide; SILAC: K:4.02511, 6.02013, 8.01420; SILAC: R:6.02013, 10.00827). Note that JUMP requires a new database for database search with dynamic modifications.   
For phosphoproteome analysis, the following dynamic modifications should be added to the parameter file.
    ```
    dynamic_S = 79.96633
    dynamic_T = 79.96633
    dynamic_Y = 79.96633
    ```
39. `enzyme_info = Tryptic KR P`  
It is used to provide the information of the enzymes (Tryptic KR P, LysC K ; ArgC R ; GluC DE).
40. `digestion = full`  
It is used to specify if the digestion is full or partial.
41. `max_mis_cleavage = 2`  
It is used to define the maximum number of miscleavage sites allowed for each peptide.
42. `min_peptide_mass = 400.0000`  
It is used to define the minimum mass of peptide database.
43. `max_peptide_mass = 6000.0000`  
It is used to define the maximum mass of peptide database.
44. `max_modif_num = 3`  
It is used to define the maximum number of modifications allowed for each peptide.
45. Parameters specifying static modifications of amino acids.
    ```
    add_Nterm_peptide = 229.162932  (for TMT-based experiment. 0.0000 for non-TMT)
    add_Cterm_peptide = 0.0000
    add_A_Alanine = 0.0000
    add_B_avg_NandD = 0.0000
    add_C_Cysteine = 0.0000
    add_D_Aspartic_Acid = 0.0000
    add_E_Glutamic_Acid = 0.0000
    add_F_Phenylalanine = 0.0000
    add_H_Histidine = 0.0000
    add_I_Isoleucine = 0.0000
    add_J_user_amino_acid = 0.0000
    add_K_Lysine = 229.162932       (for TMT-based experiment. 0.0000 for non-TMT)
    add_L_Leucine = 0.0000
    add_M_Methionine = 0.0000
    add_N_Asparagine = 0.0000
    add_O_Ornithine = 0.0000
    add_P_Proline = 0.0000
    add_Q_Glutamine = 0.0000
    add_R_Arginine = 0.0000
    add_S_Serine = 0.0000
    add_T_Threonine = 0.0000
    add_U_user_amino_acid = 0.0000  (SEQUEST search only)
    add_V_Valine = 0.0000
    add_W_Tryptophan = 0.0000
    add_X_LorI = 0.0000             (SEQUEST search only)
    add_Y_Tyrosine = 0.0000
    add_Z_avg_QandE = 0.0000        (SEQUEST search only)
    ```
46. `simulation = 0`  
It is used for testing the target-decoy strategy (0 = disable, 1 = enable).
47. `sim_MS1 = 1000`  
It defines the ppm addition for MS1 decoys.
48. `sim_MS2 = 5`  
It defines Da window for randomized MS2 peaks.
49. `temp_file_removal = 1`  
It specifies whether keeping the temporary files or removing them (0 = disable (i.e. keep temporary files), 1 = enable (remove temporary files)).

### 2.4.3 Setup the parameters for JUMP filtering
1. `[name of output directory]:[absolute path of JUMP search result]`  
To filter the JUMP search result (and perform protein grouping and so on), a user should specify the full paths of search result directory which may contain several output files such as .dta and/or .out/.spout files and .pit file which corresponds to the database used in JUMP search. For example, `HH_human_jump:/home/user/test/HH_tmt10_human/HH_tmt10_human_jump.1`
2. `unique_protein_or_peptide = protein`  
It determines the level of filtering, i.e. either protein or peptide.
3. `initial_outfile_fdr = 5`  
It defines the percentage (%) of the initial FDR for filtering scores (default = 5%).
4. `multistep_FDR_filtering = 1`  
Multistep FDR filtering (0 = disable, 1 = enable).
5. `FDR = 2`  
It defines the % of FDR for filtering peptides or one-hit-wonder proteins (fixed < 1% FDR for proteins matched by two or more precursors). Final FDR will be around 1%, when summarizing all groups of proteins/peptides.
6. `one_hit_wonders_removal = 0`  
It determines whether to keep or remove one hit wonders (-1 = removal all, 0 = no filter, 1 = partial+fully, 2 = fully).
7. `mods = 0`  
It specifies modified peptides (0 = no modification, K = Lys, STY = Phosphorylation).
8. `modpairs = 0`  
It specifies whether to keep pairs of modified and unmodified peptides or not (0 = only modified peptides, 1 = pairs).
9. `pit_file = 0`  
It specifies whether to use a custom .pit file or not (0 = use .pit file used in JUMP search, otherwise the absolute path of a custome file).
10. `min_peptide_length = 7`  
It defines the minimum peptide length (6 can be used for a small database.
11. `max_peptide_mis = 2`  
It defines the maximum number of miscleavages allowed for one peptide (default = 2).
12. `max_peptide_mod = 3`  
It defines the maximum number of modifications allowed for one peptide (M = 2, SILAC (KR) = 4, Ub = 3, Pho (STY) = 5).
13. `peptide_mod_removal = 0`  
It defines whether to remove peptides containing specific modification(s) (0 = do not remove any peptides, C = remove all C-modified peptides, STY = remove all STY-modifed peptides).
14. `peptide_aa_removal = 0`  
It defines whether to remove peptides containing specific amino acid(s) (0 = do not remove any peptides, M = remove all M-containing peptides).
15. `min_XCorr = 10`  
The minimum threshold of search scores (similar to XCorr in SEQUEST,  default = 10).
16. `min_dCn = 0`  
The minimum threshold of delta scores (similar to dCn in SEQUEST).
17. `mix_label = 0`  
It defines whether to remove mixed labeled peptides or not (0 = do not remove any peptides, KR = remove peptide labeled with SILAC, C = remove peptide labeled with ICAT, etc.).
18. `filter_contaminants = 0`  
It defines whether to remove listed contaminants or not (0 = do not remove contaminants, 1 = remove listed contaminants named with "CON_").
19. `12combinations = 1 1 1 1 1 1 1 1 0 0 0 0`  
Trypticity and charge => FT1 FT2 FT3 FT4 PT1 PT2 PT3 PT4 NT1 NT2 NT3 NT4 (1 = yes, 0 = no).
20. `bypass_filtering = 0`  
It defines whether to bypass all mass accuracy and dCn/XCorr filtering or not (0 = do all mass accuracy and dCn/XCorr filtering, 1 = bypass all filtering).
21. `mass_accuracy = 1`  
It defines whether to perform mass accuracy-based filtering or not (0 = no mass accuracy-based filtering, 1 = do filtering).
22. `mass_consideration = 1`  
It specifies the mass consideration for accuracy filtering (1 = (MH), 2 = (MH, MH+1), 3 = (MH, MH+1, MH+2), 4 = (MH, MH+1, MH+2, MH+3), 5 = (MH, MH+1, MH+2, MH+3, MH+4), 6 = (MH-1, MH, MH+1, MH+2), 7 = (MH-2, MH-1, MH, MH+1, MH+2)).
23. The parameter related with mass shift, "sd_or_static", specifies wheter mass accuracy cutoff is based on experimental standard deviation (sd_or_static = sd) or static ppm value (sd_or_static = static).
    ```
    sd_or_static = sd
    sd = 5
    ```
    ```
    sd_or_static = static
    static_cutoff = 6
    ```
24. `static_cutoff_without_mass_calib = 10`  
If there are not enough good scans, use this threshold value for ppm cut without mass calibration.
25. `FDR_filtering_method = group`  
It defines the filtering methods. Select one of the two filtering methods (LDA or group).
26. `min_outfile_num_for_XCorr_filter = 200`  
It defines the number of outfiles in each group for filtering scores; any number between 500 and 1000 is recommended.
27. `one_hit_wonders_min_XCorr_z1 = 100`  
It defines the minimum score threshold (in fact Jscore from JUMP search) for peptides with charge state of 1.
28. `one_hit_wonders_min_XCorr_z2 = 25`  
It defines the minimum score threshold for peptides with charge state of 2.
29. `one_hit_wonders_min_XCorr_z3 = 35`  
It defines the minimum score threshold for peptides with charge state of 3 or above.
30. `one_hit_wonders_min_dCn = 0.1`  
It defines the minimum delta score threshold for one-hit-wonders.
31. `one_hit_wonders_mis = 1`  
It specifies the number of miscleavages allowed for one-hit-wonders.
32. `one_hit_wonders_mods = 1`  
It specifies the number of modifications allowed for one-hit-wonders (M = 1, Ub = 2, SILAC (KR) = 3, Pho (STY) = 4).
33. `output_pepXML = 0`  
It specifies whether to enable HTML-based access or not (0 = disable HTML generation, 1 = enable HTML generation).

### 2.4.4 Setup the parameters for JUMP quantification
JUMP quantification performs the quantification of proteins/peptides from TMT-labeled experiments.

1. `idtxt = /home/user/sum_HH_tmt10_human_jump/ID.txt`  
It specifies the path of JUMP filtering result to be quantified by JUMP quantification.
2. `save_dir = HH_tmt10_human_jump`  
The name of the directory where JUMP quantification results are stored (prefix "quan-" will be added).
3. `ppi_filter = 70`  
It specifies the threshold of precursor peak intensity percentage of PSMs (don't be confused with search ppi).
4. `min_intensity_method = 1, 4`  
It defines the method of minimum intensity-based filtering of PSMs (0 = no use of the filter, 1 = minimum (intensity), 2 = maximum, 3 = mean, 4 = median. Multiple methods can be used).
5. `min_intensity_value = 1000, 5000`  
It defines the threshold value(s) of the above filter(s).
6. `min_intensity_method_1_2_psm = 1, 4`  
It defines the method of minimum intensity-based filtering of proteins identified from one or two PSMs (0 = no use of the filter, 1 = minimum (intensity), 2 = maximum, 3 = mean, 4 = median. Multiple methods can be used).
7. `min_intensity_value_1_2_psm = 2000, 10000`  
It defines the threshold value(s) of the above filter(s).
8. `impurity_correction = 1`  
It defines whether to correct the impurity of TMT reporters (0 = no, 1 = yes).
9. `impurity_matrix = /home/user/TMT10.ini`  
When `impurity_correction` is set to 1, a file containing TMT reporter impurity information should be specified by this parameter.
10. `loading_bias_correction = 1`  
It defines whether to correct the loading-biases over samples during quantification (0 = no, 1 = yes).
11. `loading_bias_correction_method = 1`  
It specifies the method of the loading-bias correction (1 = mean intensity-based, 2 = median intensity-based).
12. `SNratio_for_correction = 10`  
It defines the minimum signal-to-noise (SN) ratio for the loading-bias correction.
13. `percentage_trimmed = 25`  
It defined the percentage of most variable intensities to be trimmed for the loading-bias correction.
14. `interference_removal = 0`  
It defines whether to remove the interference in TMT quantification (i.e. ratio suppression) or not. If set to 1, it will take a significantly long time.
15. `tmt_reporters_used = sig126; sig127N; ...`  
It specifies TMT reports to be used in the quantification. For example, when all reporters of TMT-10plex are used, `tmt_reporters_used = sig126; sig127N; sig127C; sig128N; sig128C; sig129N; sig129C; sig130N; sig130C; sig131`. The reporter names should be separated by semicolon (;).
16. `tmt_peak_extraction_second_sd = 8`  
It is related with the tolerance for finding TMT reporter ion peaks.
17. `tmt_peak_extraction_method = 1`  
It defines the method of TMT reporter ion extraction method (1 = strongest peak, 2 = closest peak).
18. Sample labels need to be specified to the corresponding TMT reporters used in the quantification. Do NOT use any special character in sample labels other than underscore. For example,
    ```
    sig126 = WT_rep1
    sig127N = WT_rep2
    sig127C = WT_rep3
    sig128N = WT_rep4
    sig128C = KO1_rep1
    sig129N = KO1_rep2
    sig129C = KO1_rep3
    sig130N = KO2_rep1
    sig130C = KO2_rep2
    sig131 = KO2_rep3  
    ```
19. `comparison_analysis = 1`  
It defines whether any comparison analysis should be performed during the quantification (0 = no, 1 = yes).
20. `comparison_group_<comparison name> = ...`  
It describes how the comparison analysis should be performed. For defining comparison analyses, the prefix `comparison_groups_` and sample labels defined above should be used. Multiple comparisons are allowed. In each comparison, groups are separated by colon (:) and sample labels are separated by comma (,).  
For example, `comparison_groups_WtvsKO = WT_rep1, WT_rep2 : KO1_rep1, KO1_rep2`

# 3. References
1. Peng, J. et al. Evaluation of multidimensional chromatography coupled with tandem mass spectrometry (LC/LC−MS/MS) for large-scale protein analysis: the yeast proteome. J Proteome Res. 2003; 2: 43-50.
2. Xu, P. et al. Systematical optimization of reverse-phase chromatography for shotgun proteomics. J Proteome Res. 2009; 8: 3944-50.
3. Wang, X. et al. JUMP: a tag-based database search tool for peptide identification with high sensitivity and accuracy. Mol & Cell Proteomics. 2014;13: 3663-73.
4. Li, Y. et al. JUMPg: an integrative proteogenomics pipeline identifying unannotated proteins in human brain and cancer cells. J Proteome Res. 2016; 15: 2309-20.
5. Niu, M. et al. Extensive peptide fractionation and y1 ion-based interference detection method for enabling accurate quantification by isobaric labeling and mass spectrometry. Anal Chem. 2017; 89(5): 2956-63.
6. Taus, T. et al. Universal and confident phosphorylation site localization using phosphoRS. J Proteome Res. 2011; 10: 5354-62.
