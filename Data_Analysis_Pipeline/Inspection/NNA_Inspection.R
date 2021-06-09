################################################################################
#SCRIPT INSCPECTS ALL NNA DATA - gives quick insight into the data 
################################################################################
#want to know the median, min, max, n, UQ, LQ of new cell line 
#how do these values compare to the WT/reference cell line
#inspect stats summary tables made for new cell line and compare to WT/hypothesis
#     if an unusual/ unexpected result occurs use the New_Line to locate the cloud that relates to the unusual data point
#     visualise cloud corresponding to unusual data point in Fiji/ImageJ


#load libraries
library(tidyverse)
library(stats)


################################################################################
#USER INPUT REQUIRED

#set name of new cell line
New_Line_Name <- "Mettl3_dTAG"           #name should be spelled the same as file name within the directory

#set name of reference/control line
Reference_Line_Name <- "WT"

#define type of cells used in experiment
Cell_Type <- "mESCs"                     #either "mESCs" or "NPCs"

#define file path to where "Pulse_Chase_Analysis" is located - this is where the compiled data is stored
File_Path <- choose.dir(default = "", caption = "Select Pulse_Chase_Analysis folder, where compiled data is stored")

################################################################################
#STEP 1:upload all NNA data for a reference/control sample

All_NNA_Data <- read_csv(paste(File_Path, Cell_Type, "Nearest_Neighbour", "All_Cell_Lines_Merged", "New_All_Cell_Lines_NNA_Compile.csv", sep="/"))

#choose which cell lines to compare too
Reference_Line <- All_NNA_Data %>%
  filter(Cell_Line ==Reference_Line_Name) 


################################################################################
#STEP 2: load all NNA data relating to new cell line - chooses files which contain original data and cloud names

#contains Cell_line, Data Set, Phase, Pulse, Time for each data-point
New_Cell_Line <- read_csv(paste(File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, paste(New_Line_Name, "NNA_Compile.csv", sep="_"), sep="/"))

################################################################################
#STEP 3: calculate statistics for new and WT cell line based on
#if unusual results occur, use Cloud_Names_New_Cell_Line tibble to identify and inspect the clouds responsible

#calculate summary statistics table for new cell line NNA data grouped by phase
Phase_Stats_New_Cell_Line <- New_Cell_Line %>%      #use New_Cell_Line as contains data for all 4 phases
  group_by(Data_Set, Phase) %>%                     #groups the data based on phase and dataset
  mutate(Median = median(Distance),                 #finds the median NNA distance for each phase
         Max = max(Distance),                       #finds the maximum NNA distance for each phase
         Min = min(Distance),                       #finds the minimum NNA distance for each phase
         Upper_Quartile = quantile(Distance, 0.75), #finds the upper quartile NNA distance for each phase
         Lower_Quartile = quantile(Distance, 0.25), #finds the lower quartile NNA distance for each phase
         No_Data_Points = n()) %>%                  #finds the number of data points in each phase
  select(-Distance, -Cloud_Name, -Time) %>%         #removes distances column 
  distinct()                                        #ensures just one set of statistics for each phase

#calculate summary statistics table for new cell line NNA data grouped by time
Time_Stats_New_Cell_Line <- New_Cell_Line %>%       #use New_Cell_Line as contains data for all 4 Times
  group_by(Data_Set, Time) %>%                      #groups the data based on Time and dataset
  mutate(Median = median(Distance),                 #finds the median NNA distance for each Time
         Max = max(Distance),                       #finds the maximum NNA distance for each Time
         Min = min(Distance),                       #finds the minimum NNA distance for each Time
         Upper_Quartile = quantile(Distance, 0.75), #finds the upper quartile NNA distance for each Time
         Lower_Quartile = quantile(Distance, 0.25), #finds the lower quartile NNA distance for each Time
         No_Data_Points = n()) %>%                  #finds the number of data points in each Time
  select(-Distance, -Cloud_Name) %>%                #removes distances column 
  distinct()  

#calculate summary statistics table for WT/comparison cell line
Stats_Reference_Line <- Reference_Line %>%
  group_by(Phase) %>%                     #groups the data based on phase and dataset
  mutate(Median = median(Distance),                 #finds the median NNA distance for each phase
         Max = max(Distance),                       #finds the maximum NNA distance for each phase
         Min = min(Distance),                       #finds the minimum NNA distance for each phase
         Upper_Quartile = quantile(Distance, 0.75), #finds the upper quartile NNA distance for each phase
         Lower_Quartile = quantile(Distance, 0.25), #finds the lower quartile NNA distance for each phase
         No_Data_Points = n()) %>%                  #finds the number of data points in each phase
  select(-Distance) %>%                             #removes distances column 
  distinct()                                        #ensures just one set of statistics for each phase

################################################################################
#STEP 5: save summary stats tables

#save in csv format
Save_Path <- paste(File_Path, Cell_Type, "Nearest_Neighbour", "Stats", sep="/")

write_csv(Time_Stats_New_Cell_Line, paste(Save_Path, paste(New_Line_Name, "_NNA_Summary_Stats_Compare_Time_Points.csv", sep=""), sep="/"))
write_csv(Phase_Stats_New_Cell_Line, paste(Save_Path, paste(New_Line_Name,"_NNA_Summary_Stats_Compare_Pulses.csv", sep="" ), sep="/"))
