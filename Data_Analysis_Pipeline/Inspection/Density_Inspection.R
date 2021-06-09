################################################################################
#SCRIPT INSCPECTS ALL DENSITY DATA - gives quick insight into the data 
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
New_Line_Name <- "Test"           #name should be spelled the same as file name within the directory

#set name of reference/control line
Reference_Line_Name <- "WT"

#define type of cells used in experiment
Cell_Type <- "mESCs"                     #either "mESCs" or "NPCs"

#define file path to where "Pulse_Chase_Analysis" is located - this is where the compiled data is stored
File_Path <- choose.dir(default = "", caption = "Select Pulse_Chase_Analysis folder, where compiled data is stored")

################################################################################
#STEP 1:upload all Density data for a reference/control sample

All_Density_Data <- read_csv(paste(File_Path, Cell_Type, "Density", "All_Cell_Lines_Merged", "New_All_Cell_Lines_Density_Compile.csv", sep="/"))

#choose which cell lines to compare too
Reference_Line <- All_Density_Data %>%
  filter(Cell_Line ==Reference_Line_Name) 


################################################################################
#STEP 2: load all Density data relating to new cell line - chooses files which contain original data and cloud names

#contains Cell_line, Data Set, Phase, Pulse, Time for each data-point
New_Line <- read_csv(paste(File_Path, Cell_Type, "Density", New_Line_Name, paste(New_Line_Name, "Density_Compile.csv", sep="_"), sep="/"))


################################################################################
#STEP 3: calculate statistics for new and reference/control cell line based on

#calculate summary statistics table for new cell line Density data grouped by time and data set
Time_Stats_New_Cell_Line <- New_Line %>% 
  group_by(Data_Set, Phase, Time) %>%               #groups the data based on Phase, Time and data set
  mutate(Median = median(Distance),                 #finds the median Distance for each Time per data set
         Max = max(Distance),                       #finds the maximum Distance for each Time per data set
         Min = min(Distance),                       #finds the minimum Distance for each Time per data set
         Upper_Quartile = quantile(Distance, 0.75), #finds the upper quartile Distance for each Time per data set
         Lower_Quartile = quantile(Distance, 0.25), #finds the lower quartile Distance for each Time per data set
         No_Data_Points = n()) %>%                  #finds the number of data points in each time per data set     
  select(-Cloud_Name,-Distance, -Pulse) %>%         #removes Distances/Cloud Names column 
  distinct() %>%                                    #ensures just one set of statistics for each Time
  arrange(Time) %>%                                 #puts table in time order
  arrange(Data_Set)                                 #separates data sets within table

#calculate summary stats for each pulse per data set

Pulse_Stats_New_Cell_Line <- New_Line %>% 
  group_by(Data_Set, Phase, Pulse) %>%              #groups the data based on Phase, Pulse and data set
  mutate(Median = median(Distance),                 #finds the median Distance for each Pulse per data set
         Max = max(Distance),                       #finds the maximum Distance for each Pulse per data set
         Min = min(Distance),                       #finds the minimum Distance for each Pulse per data set
         Upper_Quartile = quantile(Distance, 0.75), #finds the upper quartile Distance for each Pulse per data set
         Lower_Quartile = quantile(Distance, 0.25), #finds the lower quartile Distance for each Pulse per data set
         No_Data_Points = n()) %>%                  #finds the number of data points in each Pulse per data set 
  select(-Cloud_Name,-Distance, -Time) %>%          #removes Distances/Cloud Names column 
  distinct() %>%                                    #ensures just one set of statistics for each Pulse
  arrange(Phase) %>%                                #puts table in Pulse order
  arrange(Data_Set)                                 #separates data sets within table


#calculate summary statistics table for reference/control/comparison cell line
Stats_Reference_Line <- Reference_Line %>%
  mutate(Median = median(Distance),                 #finds the median  Distance over all time points
         Max = max(Distance),                       #finds the maximum  Distance over all time points
         Min = min(Distance),                       #finds the minimum  Distance over all time points
         Upper_Quartile = quantile(Distance, 0.75), #finds the upper quartile  Distance over all time points
         Lower_Quartile = quantile(Distance, 0.25), #finds the lower quartile  Distance over all time points
         No_Data_Points = n()) %>%                  #finds the number of data points in each phase  
  select(-Distance) %>%                             #removes Distances column 
  distinct()                                        #ensures just one set of statistics over all time points

################################################################################
#STEP 4: look at generated stats table --> if any result seems unusual use the New_Cell_Line table to identify the cloud + view in Fiji/ImageJ

################################################################################
#STEP 5: save summary stats tables

Save_Path <- paste(File_Path, Cell_Type, "Density", "Stats", sep="/")

write_csv(Time_Stats_New_Cell_Line, paste(Save_Path, paste(New_Line_Name, "_Density_Summary_Stats_Compare_Time_Points.csv", sep=""), sep="/"))
write_csv(Pulse_Stats_New_Cell_Line, paste(Save_Path, paste(New_Line_Name,"_Density_Summary_Stats_Compare_Pulses.csv", sep="" ), sep="/"))

