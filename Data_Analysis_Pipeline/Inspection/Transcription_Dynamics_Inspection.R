################################################################################
#SCRIPT INSCPECTS TRANSCRIPTION DYNAMICS DATA - gives quick insight into the data 
#please dontact Holly Roach at hmroach@hotmail.co.uk if you have any questions
################################################################################
#returns basic summary stats for transcription dynamics
#how do these values compare to the WT/reference cell line
#inspect stats summary tables made for new cell line and compare to WT/hypothesis
#     if an unusual/ unexpected result occurs use the New_Cell_Line to locate the cloud name/image that relates to the unusual data point
#     visualise cloud corresponding to unusual data point in Fiji/ImageJ


# #load libraries
# library(tidyverse)
# library(stats)
# library(tcltk)
# 
# #create function to replace choose.dir so that it is compatible with mac OS
# choose_dir <- function(caption = 'Select data directory') {
#   if (exists('utils::choose.dir')) {
#     choose.dir(caption = caption)
#   } else {
#     tk_choose.dir(caption = caption)
#   }
# }


################################################################################
#USER INPUT REQUIRED

# #set name of new cell line - should be spelled the same as folder name within the directory
# #make sure not to use any spaces or symbols other than an underscore (_), dash (-), or full stop (.)
# #eg "Mettl3_dTAG" or "SPEN_RRM_del" or "Ciz1_KO"
# New_Line_Name <- "Test"  
# 
# #set name of reference/control line 
# #make sure spelling of this name is the same as what appears in the Pulse_Chase_Analysis\mESCs\Molecule_Count\All_Cell_Lines folder
# Reference_Line_Name <- "WT"
# 
# #define type of cells used in experiment - either "mESCs" or "NPCs"
# Cell_Type <- "mESCs"


################################################################################
#OTHER INPUTS
#the below inputs should not need changing

#define file path to where "Pulse_Chase_Analysis" is located
File_Path <- Output_File_Path


################################################################################
#VALIDATE INPUTS

#checks name of cell type
if (!(Cell_Type == "mESCs" | Cell_Type == "NPCs")) {
  stop("Invalid name entered in Cell_Type - must be mESCs or NPCs")
}

#checks correct folder has been selected for file path
if (!(str_sub(File_Path, -20, -1) == "Pulse_Chase_Analysis")) {
  stop("Pulse_Chase_Analysis folder not selected for File_Path")
}

#checks that Pulse_Chase_Analysis folder contains folder for new cell line
if (!dir.exists(paste(File_Path, Cell_Type, "Transcription_Dynamics", New_Line_Name, sep="/"))) {
  stop(paste("Folder for the new cell line (", New_Line_Name, " - New_Line_Name) does not exist in Pulse_Chase_Analysis folder", sep=""))
}

################################################################################
#STEP 1:upload Molecule_Count data for new cell line and the reference/control sample

Reference_Line <- read_csv(paste(File_Path, Cell_Type, "Transcription_Dynamics", "All_Cell_Lines", paste(Reference_Line_Name, "Transcription_Dynamics_Compile.csv", sep="_"), sep="/"))

New_Cell_Line <- read_csv(paste(File_Path, Cell_Type, "Transcription_Dynamics", "All_Cell_Lines", paste(New_Line_Name, "Transcription_Dynamics_Compile.csv", sep="_"), sep="/"))


################################################################################
#STEP 2: create stat summary for the number of Molecule_Count present in each pulse per phase and time point

#stat summary for reference dynamic initiation dataset
transcription_stats <- Reference_Line %>% 
  bind_rows(New_Cell_Line) %>% 
  group_by(Cell_Line, Phase, Time) %>%                     #group data
  summarise(n=n(),                                         #find number of data points
            Mean=mean(No_Centroids),                       #finds the mean 
            Median = median(No_Centroids),                 #find the median
            Upper_Quartile = quantile(No_Centroids, 0.75), #find the upper quartile
            Lower_Quartile = quantile(No_Centroids, 0.25), #find the lower quartile
            SD=sd(No_Centroids)) %>%                       #finds standard deviation
  mutate(SE = SD/sqrt(n),                                  #finds standard error of mean
         CI_95 = SE * qt((1-0.05)/2 + .5, n-1)) %>%        #finds 95% confidence limit
  ungroup() %>% 
  arrange(Phase, Time)                                     #orders data
                                      

################################################################################
#STEP 3: look at generated stats table --> if any result seems unusual use the New_Cell_Line table to identify the cloud + view in Fiji/ImageJ

################################################################################
#STEP 4: save summary stats table for new cell line

Save_Path <- paste(File_Path, Cell_Type, "Transcription_Dynamics", New_Line_Name, sep="/")
File_Name <-paste(New_Line_Name, "Transcription_Dynamics_Summary_Stats.csv", sep="_")

write_csv(transcription_stats, paste(Save_Path, File_Name, sep="/"))
