################################################################################
#SCRIPT INSCPECTS ALL  Centriod DATA - gives quick insight into the data 
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
#STEP 1:upload all centroid data for a reference/control sample

All_Centroid_Data <- read_csv(paste(File_Path, Cell_Type, "Turnover", "All_Cell_Lines_Merged", "New_All_Cell_Lines_Centroids_Compile.csv", sep="/"))

#choose which cell lines to compare too
Reference_Line <- All_Centroid_Data %>%
  filter(Cell_Line ==Reference_Line_Name) 


################################################################################
#STEP 2: load all centoits data relating to new cell line - chooses files which contain original data and cloud names

New_Cell_Line <- read_csv(paste(File_Path, Cell_Type, "Turnover", New_Line_Name, paste(New_Line_Name, "Centroids_Compile.csv", sep="_"), sep="/"))


################################################################################
#STEP 4: calculate number of clouds with few centroids in pulse 1 or pulse 2

Ref_Few_Pulse_1_Centroids <- filter(Reference_Line, Pulse == "Pulse_1") %>%
  filter(No_Centroids < 5) 

Ref_Few_Pulse_2_Centroids <- filter(Reference_Line, Pulse == "Pulse_2") %>%
  filter(No_Centroids < 5) 

New_Line_Few_Pulse_1_Centroids <- filter(New_Cell_Line, Pulse == "Pulse_1") %>%
  filter(No_Centroids < 5) 

New_Line_Few_Pulse_2_Centroids <- filter(New_Cell_Line, Pulse == "Pulse_2") %>%
  filter(No_Centroids < 5) 

Few_Centroids <- New_Line_Few_Pulse_1_Centroids %>%    #combines all data points which record Few centroids
  bind_rows(New_Line_Few_Pulse_2_Centroids) %>%
  select(-Cloud_Name,) %>%
  bind_rows(Ref_Few_Pulse_1_Centroids,
            Ref_Few_Pulse_2_Centroids) 

#see which data sets have few centroids - is this different to the Ref?
Compare_Few_Centroids_Pulse_Stats <- Few_Centroids %>%
  group_by(Cell_Line, Data_Set, Phase, Pulse) %>%             #groups the data based on Cell_Line, Dataset, phase and pulse
  mutate(No_Data_Points = n()) %>%                            #finds the number of data points per phase in each Time per data set
  select(-No_Centroids, -Time) %>%                            #removes No_Centroids column 
  distinct() %>%                                              #ensures just one set of statistics for each Time per pulse/data set
  arrange(Data_Set,
          Phase,
          Pulse) 


################################################################################
#STEP 4: calculate statistics for new and Ref cell line based on

#calculate summary statistics table for new cell line Centroid data grouped by time, phase and data set
Time_Stats_New_Cell_Line <- New_Cell_Line %>% 
  group_by(Time, Data_Set, Phase) %>%                         #groups the data based on Time and data set
  mutate(Median = median(No_Centroids),                       #finds the median  No_Centroids for each Time per data set
         Max = max(No_Centroids),                             #finds the maximum  No_Centroids for each Time per data set
         Min = min(No_Centroids),                             #finds the minimum  No_Centroids for each Time per data set
         Upper_Quartile = quantile(No_Centroids, 0.75),       #finds the upper quartile  No_Centroids for each Time per data set
         Lower_Quartile = quantile(No_Centroids, 0.25),       #finds the lower quartile  No_Centroids for each Time per data set
         No_Data_Points = n()) %>%                            #finds the number of data points per phase in each Time per data set
  select(-No_Centroids, -Cell_Line, -Cloud_Name, -Pulse) %>%   #removes No_Centroids column 
  distinct() %>%                                              #ensures just one set of statistics for each Time per pulse/data set
  arrange(Time) %>%
  arrange(Data_Set)

#calculate summary stats for new cell line data based on phase, pulse, data set
Phase_Pulse_Stats_New_Cell_Line <- New_Cell_Line %>%  
  group_by(Phase, Pulse, Data_Set) %>%                        #groups the data based on Phase and Pulse
  mutate(Median = median(No_Centroids),                       #finds the median  No_Centroids for each Time per data set
         Max = max(No_Centroids),                             #finds the maximum  No_Centroids for each Time per data set
         Min = min(No_Centroids),                             #finds the minimum  No_Centroids for each Time per data set
         Upper_Quartile = quantile(No_Centroids, 0.75),       #finds the upper quartile  No_Centroids for each Time per data set
         Lower_Quartile = quantile(No_Centroids, 0.25),       #finds the lower quartile  No_Centroids for each Time per data set
         No_Data_Points = n()) %>%                            #finds the number of data points per phase in each Time per data set
  select(-No_Centroids, -Cell_Line, -Cloud_Name, -Time) %>%   #removes No_Centroids column 
  distinct()                                                  #ensures just one set of statistics for each Phase


#calculate summary statistics table for Ref/comparison cell line
Time_Stats_Reference_Line <- Reference_Line %>% 
  group_by(Time, Data_Set, Phase) %>%                         #groups the data based on Time and data set
  mutate(Median = median(No_Centroids),                       #finds the median  No_Centroids for each Time per data set
         Max = max(No_Centroids),                             #finds the maximum  No_Centroids for each Time per data set
         Min = min(No_Centroids),                             #finds the minimum  No_Centroids for each Time per data set
         Upper_Quartile = quantile(No_Centroids, 0.75),       #finds the upper quartile  No_Centroids for each Time per data set
         Lower_Quartile = quantile(No_Centroids, 0.25),       #finds the lower quartile  No_Centroids for each Time per data set
         No_Data_Points = n()) %>%                            #finds the number of data points per phase in each Time per data set
  select(-No_Centroids, -Cell_Line, -Pulse) %>%               #removes No_Centroids column 
  distinct() %>%                                              #ensures just one set of statistics for each Time per pulse/data set
  arrange(Time) %>%
  arrange(Data_Set)


Phase_Pulse_Stats_Reference_Line <-Reference_Line %>% 
  group_by(Phase, Pulse, Data_Set) %>%                        #groups the data based on Phase and Pulse
  mutate(Median = median(No_Centroids),                       #finds the median  No_Centroids for each Time per data set
         Max = max(No_Centroids),                             #finds the maximum  No_Centroids for each Time per data set
         Min = min(No_Centroids),                             #finds the minimum  No_Centroids for each Time per data set
         Upper_Quartile = quantile(No_Centroids, 0.75),       #finds the upper quartile  No_Centroids for each Time per data set
         Lower_Quartile = quantile(No_Centroids, 0.25),       #finds the lower quartile  No_Centroids for each Time per data set
         No_Data_Points = n()) %>%                            #finds the number of data points per phase in each Time per data set
  select(-No_Centroids, -Cell_Line, -Time) %>%                #removes No_Centroids column 
  distinct()                                                  #ensures just one set of statistics for each Phase

################################################################################
#STEP 5: look at generated stats table --> if any result seems unusual use the New_Cell_Line table to identify the cloud + view in Fiji/ImageJ

################################################################################
#STEP 5: save summary stats tables

Save_Path <- paste(File_Path, Cell_Type, "Turnover", "Stats", sep="/")

write_csv(Time_Stats_New_Cell_Line, paste(Save_Path, paste(New_Line_Name,"Time_Centroids_Summary_Stats.csv", sep="_"), sep="/"))
write_csv(Phase_Pulse_Stats_New_Cell_Line, paste(Save_Path, paste(New_Line_Name,"Phase_Pulse_Centroids_Summary_Stats.csv", sep="_"), sep="/"))
write_csv(Compare_Few_Centroids_Pulse_Stats, paste(Save_Path, paste(New_Line_Name, "Number_Data_Points_Few_Centroids_Stats.csv", sep="_"), sep="/"))


