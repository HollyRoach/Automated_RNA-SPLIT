####################################################################################################################
#SCRIPT USED TO PERFORM RBNA-SPLIT DATA ANALYSIS PIPELINE
#please contact Holly Roach at hmroach@hotmail.co.uk if you have any questions
####################################################################################################################
#IMPORTANT - before running this script make sure you have done the following:
# 1 - In the Watershed_Algorithm_Results folder make sure you have a folder corresponding to the new of your new cell line
#     this file should contain the raw data files for each cropped cloud image that were created by the algorithm on the micron server
#     within this folder there should be subfolders for the Data_Set, Phase and Timepoints used in the analysis
# 2 - makes sure here, tidyverse, stringr, RColorBrewer and tcltk packages have been installed
#     if these are not installed you can run - install.packages("package_name") - in the console
# 3 - make sure you have ran the ImageJ macro (ImageJ_Cloud_Volume.ijm) prior to running this script
# 4 - within the Pulse_Chase_Analysis folder make sure you have created a specific for your new cell line in the Nearest_Neighbour subfolder
#     this folder should contain a nascent_Xist_dynamics and a Xist_turnover_on_chromatin folder
#     within each dataset folder, there should be an EdU_Control and Random_Control folder which contain files corresponding to these controls
#     this only needs to be done for the NNA as subfolders for all other parameters are generated for you by the script
# 5 - make sure that the user inputs in lines 22-96 have been made compatible for your analysis

#MAKE SURE YOU HAVE FOLLOWED THE ABOVE INSTRUCTIONS BEFORE ATTEMPTING TO RUN THIS SCRIPT

#To run this script, press Ctrl+Shift+S (or Command+Shift+S on Mac)

################################################################################
#USER INPUT REQUIRED

#set name of new cell line - should be spelled the same as folder name within the directory
#make sure not to use any spaces or symbols other than an underscore (_), dash (-), or full stop (.)
#eg "Mettl3_dTAG" or "SPEN_RRM_del" or "Ciz1_KO"
New_Line_Name <- "Test"           

#define type of cells used in experiment - either "mESCs" or "NPCs"
Cell_Type <- "mESCs"

#define time points used in both datasets
#make sure that the lists include all the possible time-points for your data
#it is ok if the lists contain extra time-points, as these lists are a superset of time points
#it doesn't matter if initiation/maintenance have different time points
Dynamic_Time <- c(10, 20, 30, 40, 50, 60)
Turnover_Time <- c(0, 60, 80, 100, 120, 140, 160, 180, 200, 220)

#set name of reference/control line 
#make sure spelling of this name is the same as what appears in the Pulse_Chase_Analysis\mESCs\Parameter\All_Cell_Lines folder
Reference_Line_Name <- "WT"

Cell_Line_3 <- "Ciz1_KO"

Cell_Line_4 <- "SPEN_RRM_del"  

#define the threshold used to determine a significant result during wilcoxon test
p <- 0.05

#define colours for each cell line 
#set colour for the reference line (pink)
Colour_1 <- "#F768A1"        #sets colour for reference line
Colour_1_EdU <- "#FDE0DD"    #sets colour for EdU control of Reference_Line_Name
Colour_1_exp <- "#FA9FB5"    #sets colour for expansion phase Reference_Line_Name
Colour_1_ss <- "#DD3497"     #sets colour for steady state phase Reference_Line_Name
Colour_1_Random <- "#FCC5C0" #sets colour for Random control of Reference_Line_Name

#set colour for the new cell line (blue)
Colour_2 <- "#67A9CF"        #sets colour for new cell line
Colour_2_EdU <- "#ECE2F0"    #sets colour for EdU control of new cell line
Colour_2_exp <- "#A6BDDB"    #sets colour for expansion phase new cell line
Colour_2_ss <- "#3690C0"     #sets colour for steady state phase new cell line
Colour_2_Random <- "#D0D1E6" #sets colour for Random control of new cell line

#set colour for the cell line 3 (orange)
Colour_3 <- "#FD8D3C"     #sets colour for cell_line_3
Colour_3_exp <- "#FEB24C" #sets colour for expansion phase new cell_line_3
Colour_3_ss <- "#FC4E2A"  #sets colour for steady state phase new cell_line_3

#set colour for the cell line 4 (green)
Colour_4 <- "#78C679"     #sets colour for cell_line_4
Colour_4_exp <- "#ADDD8E" #sets colour for expansion phase new cell_line_3
Colour_4_ss <- "#41AB5D"  #sets colour for steady state phase new cell_line_3

#set name of how cell lines should be presented in the plot - ensures name consistency with other papers
#if name contains delta/triangle symbol use - "SPEN^"~Delta*"RRM"
Name_1 <- bquote("METTL3_FKBP12"^"F36V") #sets name for new cell line
Name_2 <- bquote("WT")                   #sets name for Reference_Line_Name
Name_3 <- bquote("Ciz1_KO")              #sets name for cell_line_3
Name_4 <- bquote("SPEN"^~Delta*"RRM")    #sets name for cell_line_4

#select which time points to present for transcription dynamics
Timepoints <- c("20", "40", "60")

#define which dataset you want to present for the nearest neighbour analysis
New_Line_Data <- "Turnover" #should be "Turnover" or "Dynamic"
Ref_Line_Data <- "Dynamic" #should be "Turnover" or "Dynamic"

################################################################################
#load libraries
library(here)
library(tidyverse)
library(stringr)
library(tcltk)
library(RColorBrewer)

#load functions used in scripts
source(here("utils.R"))

################################################################################
#OTHER INPUTS
#the below inputs should not need changing

#creates list of datasets being used for compile
Data_Set_List <- c("nascent_Xist_dynamics",
                   "Xist_turnover_on_chromatin")

#creates list of phases being used for compile
Phase_List <- c("Initiation", "Maintenance")

#define file path to where Watershed Algorithm "results" are found - needs to be located in documents due to long file path names
Input_File_Path <- choose_dir(caption = "Select (open on Mac) Watershed_Algorithm_Results folder, where raw data is found")

#define file path to where "Pulse_Chase_Analysis" is located - results from this scripts will be stored here
Output_File_Path <- choose_dir(caption = "Select (open on Mac) Pulse_Chase_Analysis folder, where compiled data will be stored")


################################################################################
#Run scripts to compile and manipulate all the data
source(here("Compilation_Manipulation", "C1_C2-Centroids_Compilation_Manipulation.R"))
source(here("Compilation_Manipulation", "Cloud_Volume_Compilation_Manipulation.R"))
source(here("Compilation_Manipulation", "Density_Compilation_Manipulation.R"))
source(here("Compilation_Manipulation", "NNA_Compilation_Manipulation.R"))

################################################################################
#Run scripts to inspect all the data
source(here("Inspection", "Centroid_Inspection.R"))
source(here("Inspection", "Cloud_Volume_Inspection.R"))
source(here("Inspection", "Density_Inspection.R"))
source(here("Inspection", "NNA_Inspection.R"))
source(here("Inspection", "Total_Molecule_Count_Inspection.R"))
source(here("Inspection", "Transcription_Dynamics_Inspection.R"))

################################################################################
#Run scripts to plot and prefrom unpaired, two-tailed wilcoxon stats test on all the data
print("Overall outputs from stats tests")

source(here("Analysis", "Cloud_Volumes_Box_Plot_Wilcoxon.R"))
source(here("Analysis", "Density_Box_Plot_Wilcoxon.R"))
source(here("Analysis", "Molecule_Count_Box_Plot_Wilcoxon.R"))
source(here("Analysis", "Transcription_Box_Plot_Wilcoxon.R"))
source(here("Analysis", "NNA_Violin_Plot_Wilcoxon.R"))
