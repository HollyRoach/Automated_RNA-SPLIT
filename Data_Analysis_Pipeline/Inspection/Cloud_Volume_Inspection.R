################################################################################
#SCRIPT INSCPECTS ALL Cloud Volume DATA - gives quick insight into the data 
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
#STEP 1:upload all Cloud_Volume data for a reference/control sample

All_Cloud_Volume_Data <- read_csv(paste(File_Path, Cell_Type, "Cloud_Volume", "All_Cell_Lines_Merged", "New_All_Cell_Lines_Cloud_Volume_Compile.csv", sep="/"))

#choose which cell lines to compare too
Reference_Line <- All_Cloud_Volume_Data %>%
  filter(Cell_Line ==Reference_Line_Name) 


################################################################################
#STEP 2: load all Cloud_Volume data relating to new cell line - chooses files which contain original data and cloud names

#contains Cell_line, Data Set, Phase, Pulse, Time for each data-point
New_Cell_Line <- read_csv(paste(File_Path, Cell_Type, "Cloud_Volume", New_Line_Name, paste(New_Line_Name, "Cloud_Volume_Compile.csv", sep="_"), sep="/"))


################################################################################
#STEP 3: calculate statistics for new and WT cell line based on
#look at tables generated 
#if unusual results occur, use New_Cell_Line tibble to identify and inspect the clouds responsible

#calculate summary statistics table for new cell line Cloud volume data grouped by phase and data set
Phase_Stats_New_Cell_Line <- New_Cell_Line %>%          
  group_by(Phase, Data_Set) %>%                         #groups the data based on phase and data set
  mutate(Median = median(Cloud_Volume),                 #finds the median  Cloud_Volume for each phase
         Max = max(Cloud_Volume),                       #finds the maximum  Cloud_Volume for each phase
         Min = min(Cloud_Volume),                       #finds the minimum  Cloud_Volume for each phase
         Upper_Quartile = quantile(Cloud_Volume, 0.75), #finds the upper quartile  Cloud_Volume for each phase
         Lower_Quartile = quantile(Cloud_Volume, 0.25), #finds the lower quartile  Cloud_Volume for each phase
         No_Data_Points = n()) %>%                      #finds the number of data points in each phase
  select(-Cloud_Name, -Time, -Cloud_Volume) %>%         #removes unnessasry columns
  distinct()                                            #ensures just one set of statistics for each phase

#calculate summary statistics table for new cell line Cloud volume data grouped by phase and data set
Time_Stats_New_Cell_Line <- New_Cell_Line %>% 
  group_by(Time, Phase, Data_Set) %>%                     #groups the data based on Time and data set
  mutate(Median = median(Cloud_Volume),                   #finds the median  Cloud_Volume for each Time
         Max = max(Cloud_Volume),                         #finds the maximum  Cloud_Volume for each Time
         Min = min(Cloud_Volume),                         #finds the minimum  Cloud_Volume for each Time
         Upper_Quartile = quantile(Cloud_Volume, 0.75),   #finds the upper quartile  Cloud_Volume for each Time
         Lower_Quartile = quantile(Cloud_Volume, 0.25),   #finds the lower quartile  Cloud_Volume for each Time
         No_Data_Points = n()) %>%                        #finds the number of data points in time
  select(-Cloud_Name, -Cloud_Volume) %>%                  #removes unnessasry columns
  distinct() %>%                                          #ensures just one set of statistics for each Time
  arrange(Time) %>%                                       #puts data in time order
  arrange(Data_Set)                                       #puts data in order of Data Set

#calculate summary statistics table for WT/comparison cell line
Stats_Reference_Line <- Reference_Line %>%
  group_by(Phase) %>%                                       #groups the data based on phase
  mutate(Median = median(Cloud_Volume),                     #finds the median  Cloud_Volume for each phase
         Max = max(Cloud_Volume),                           #finds the maximum  Cloud_Volume for each phase
         Min = min(Cloud_Volume),                           #finds the minimum  Cloud_Volume for each phase
         Upper_Quartile = quantile(Cloud_Volume, 0.75),     #finds the upper quartile  Cloud_Volume for each phase
         Lower_Quartile = quantile(Cloud_Volume, 0.25),     #finds the lower quartile  Cloud_Volume for each phase
         No_Data_Points = n()) %>%                          #finds the number of data points in each phase) %>%                 
  select(-Cloud_Volume) %>%                                 #removes Cloud_Volumes column 
  distinct()                                                #ensures just one set of statistics for each phase

################################################################################
#STEP 5: save summary stats tables for new cell line

Save_Path <- paste(File_Path, Cell_Type, "Cloud_Volume", "Stats", sep="/")

write_csv(Time_Stats_New_Cell_Line, paste(Save_Path, paste(New_Line_Name, "_Cloud_Volume_Summary_Stats_Compare_Time_Points.csv", sep=""), sep="/"))
write_csv(Pulse_Stats_New_Cell_Line, paste(Save_Path, paste(New_Line_Name,"_Cloud_Volume_Summary_Stats_Compare_Pulses.csv", sep="" ), sep="/"))
