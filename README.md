This repository contains codes and files used for the project.
This readme explain how each part of the repository is used.
NOTE REMEMBER TO CHANGE PATH TO DIRECTORY TO ENABLE CODES RUN PROPERLY IN:
mfe_output_processing.py in the nuapck output folders

####STRUCTURAL ANALYSIS####
Codes to perform secondary structure analysis are in:
MRes_SSB_thesis/mrna_stability_analysis/SS/

monocistronic denotes codes that use data from RegulonDB
kim2008 denotes codes that use data from Kim et al. 2008

#JUST OUTPUT#
If just want to see output, then the "_nupack_input" folders can be ignored
So can "mfe_output_processing.py" scripts in each "nupack_output" folder
Run the "SSvsHL.py" scripts in the "_output" folders. All necessary datasets are present in these folders.

#GENERATING DATSETS FOR NUPACK
If want to see how NUPACK inputs were generated navigate to the "_nupack_input" folders
Run the single python script in there to generate a csv of UTR sequences
1 csv will be generated in the 'kim_2008_nupack_input" move this to the "mfe" folder
Then run "csv_to_txt.py" to generate .IN files which contain sequence infromaiton, which can be input into NUPACK in a linux server

2 csvs will be generated in the 'monocistronic_nupack_input"
A 5'utr csv and 3' UTR csv move these to the "5" and "3" folders respecitvely
The run "csv_to_txt.py" to generate .IN files which contain sequence infromaiton, which can be input into NUPACK in a linux server
Then run "csv_to_txt.py" to generate .IN files which contain sequence infromaiton, which can be input into NUPACK in a linux server

Once the output has been generated from NUPACK move the reuslting files into:
Kim sequences into "kim2008_output/mfe/" and run "mfe_output_processing.py" in "kim2008_output/" to generate a csv of secodnary sturctrues and mfes
3' UTR monocistronic sequence into "monocistronic_output/3/mfe/" and run "mfe_output_processing.py" in "monocistronic_output/" to generate a csv of secodnary sturctrues and mfes
5' UTR monocistronic sequence into "monocistronic_output/5/mfe/" and run "mfe_output_processing.py" in "monocistronic_output/" to generate a csv of secodnary sturctrues and mfes

####TIR####
#NATURAL GENE OUTPUT#
To get output of the natural transcript mRNA TIR vs half life
Navigate to mrna_stability_analysis/TIR/
Run "TIRvsHL.py"
If want to see how the datasets for this were generated downlaod the CSV files  from the links in "LINK TO GOOGLE DRIVE FOR RBS DATASETS" and place in this folder
And run RBS_calc_output_processing.py
NOTE: this script processes the output from the RBS calculator, NOT generate input csvs for calculator

#SYNTHETIC GENE OUTPUT#
To get output of the syntheitc mRNA TIR vs half life
Navigate to mrna_stability_analysis/TIR/Cambray/
Downlaod the CSV files  from the links in "UPLOAD LINK FOR CAMBRAY DATASETS" and place in this folder
And TIRvsHL.py

####ANALYSING FRAMEWORK####
#USING COPASI#
To perform COPASI analysis download COPASI
Navigate to system_model_analysis/COPASI_analysis/copasi_model/
To do ON state analysis transfer the circularisation.xml file to copasi
To do OFF state analysis transfer the circularisation_linearisation.xml file to copasi

#ANALYSING COPASI OUTPUT#
If want to see analysis of COPASI output 
Navigate to system_model_analysis/COPASI_analysis/copasi_data/
PSA: copasi_PSA.py
Parameter scan: copasi_scan.py
Sampling: copasi_sample.py

#TRANSLATION MODEL#
To run infintie ORF and mRNA translation models
Navugate to system_model_analysis/mRNA_translation_model/python/
And run the python scripts inside

#FRAMEWORK MODEL#
To run stochastic and deterministic simualtions of the framework
Navigate to system_model_analysis/Python_models_of_system
And run the python scripts inside


