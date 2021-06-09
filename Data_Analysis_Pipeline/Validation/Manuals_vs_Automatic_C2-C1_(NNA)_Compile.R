################################################################################
#MANUAL VS AUTOMATIC C2-C1 COMPILE VALIDATION
################################################################################
#for WT data only use nascent dynamics data set to calculate NNA!
#
#load libraries
library(extrafont)
font_import() #need to do once at beginning before the rest of the script can run
loadfonts(device = "win") #loads windows fonts to use
library(tidyverse)
library(stats)
library(RColorBrewer)
library(stringr)
library(fuzzyjoin)
library(readxl)

#create theme for plots - defines font, text size etc
theme <- theme(plot.title = element_text(family = "Calibri (Body)", face = "bold", size = (13)),
               legend.title = element_text(family = "Calibri (Body)", size = (11)), 
               legend.text = element_text(family = "Calibri (Body)", size = (11)), 
               axis.title = element_text(family = "Calibri (Body)", face = "bold", size = (11)),
               axis.text = element_text(family = "Calibri (Body)", face = "bold", size = (7)))

################################################################################
#USER INPUT REQUIRED

#set name of new cell line
New_Line_Name <- "WT_auto"           #name should be spelled the same as file name within the directory


#define file path to where "Pulse_Chase_Analysis" is located - this is where the compiled data is stored
File_Path <- "C:\\Users\\hmroa\\Documents\\RNA-SPLIT_Results\\Manual_vs_Auto\\Pulse_Chase_Analysis"

#define type of cells used in experiment
Cell_Type <- "mESCs"                      #either "mESCs" or "NPCs"


################################################################################
#STEP 1: upload all C2-C1 data for WT sample

setwd(paste(File_Path, Cell_Type, "Nearest_Neighbour", "All_Cell_Lines_Merged", "nascent_Xist_dynamics", sep="\\"))

# #loads files Lisa originally gave me
# WT_Original <- read_csv("Cleaned_WT_NNA_2to1_Compile.csv") %>%
#   filter(Cell_Line == "WT") %>%
#   mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_original")) %>% 
#   filter(Phase == "Initiation" | Phase == "Maintenance")



#loads new data for WT manual
WT_Manual_Line <- read_csv("WT_Original_Data_NNA_2to1_Compile.csv") %>%
  mutate(Cell_Line = "WT_manual") %>% 
  filter(Phase == "Initiation" | Phase == "Maintenance")            


################################################################################
#STEP 2: upload WT automated compile to inspect - tibble includes cloud names

setwd(paste(File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, "nascent_Xist_dynamics", sep="\\"))

WT_Auto_Initation <- read_csv("All_C2-C1_NNA_Data_Initiation_ALL_DATA.csv") %>%
  transmute(Cell_Line = "WT_auto",
            Phase = "Initiation",
            Time = Time,
            Cloud_Name = Cloud_Name,
            Distance = NN_dist_2to1)

WT_Auto_Maintenance <- read_csv("All_C2-C1_NNA_Data_Maintenance_ALL_DATA.csv") %>%
  transmute(Cell_Line = "WT_auto",
            Phase = "Maintenance",
            Time = Time,
            Cloud_Name = Cloud_Name,
            Distance = NN_dist_2to1)

WT_Auto_Line <- WT_Auto_Initation %>%
  bind_rows(WT_Auto_Maintenance)

#outcomes/observations
#the wT_Auto_Line contains 6,065 observations compared to 5,787 observation in the WT_Manual_Line (manual compiles)
#  --> wT_Auto_Line has 278 extra data points
#  --> extra data in wT_Auto_Line probably corresponds to data from clouds which were excluded during the manual compile


################################################################################
#STEP 3: What is the same and different about the manual (WT_Manual_Line) and automated (WT_Auto_Line) compiles?

#finds which rows are the same in the manual and auto compile
#  --> does create duplicates for when the same Distance ar recorded in the manual for the same, Time, Phase, Pulse and Dataset
same_to_4dp <- inner_join(mutate(WT_Auto_Line, Distance.auto = Distance,
                          Distance.auto.4dp = round(Distance, digits = 3)),
                   mutate(WT_Manual_Line, Distance.manual = Distance,
                          Distance.manual.4dp = round(Distance, digits = 3)),
                   by= c("Phase",
                         "Distance.auto.4dp" = "Distance.manual.4dp")) %>% 
  distinct(Phase,
           Time,
           Cloud_Name,
           Distance.auto, .keep_all=TRUE)


#finds which observations are only present in the automated compile
only_in_auto <- anti_join(mutate(WT_Auto_Line, Distance.auto = Distance,
                                 Distance.auto.4dp = round(Distance, digits = 4)),
                          mutate(WT_Manual_Line, Distance.manual = Distance,
                                 Distance.manual.4dp = round(Distance, digits = 4)),
                          by= c("Distance.auto.4dp" = "Distance.manual.4dp"))

#finds which observations are only present in the automated compile
only_in_manual <- anti_join(mutate(WT_Manual_Line, Distance.manual = Distance,
                                   Distance.manual.4dp = round(Distance, digits = 4)),
                            mutate(WT_Auto_Line, Distance.auto = Distance,
                                   Distance.auto.4dp = round(Distance, digits = 4)),
                            by= c("Distance.manual.4dp" = "Distance.auto.4dp"))

#find if there are any very similar datapoints between manual and auto compile
similar_to_1dp <- inner_join(mutate(only_in_auto, Distance.auto.1dp = round(Distance, digits = 1)),
                      mutate(only_in_manual, Distance.manual.1dp = round(Distance, digits = 1)),
                            by= c("Phase",
                                  "Distance.auto.1dp" = "Distance.manual.1dp")) %>% 
                      distinct(Phase,
                               Time,
                               Cloud_Name,
                               Distance.manual, .keep_all=TRUE) %>%
                      select(-Distance.x, -Distance.auto.4dp, -Distance.y, -Distance.manual.4dp) %>%
                      mutate(Difference = Distance.auto - Distance.manual) %>%    #calculates the differnce in original distance values (not rounded)
                      mutate(Difference_squ = Difference*Difference)              #need to remove negatives

#outcomes/observations
#out of the 6,065 observations in the WT_Auto_Line, 5,750 are the same as the WT_Manual_Line
#there are 341 observations that are only present in WT_Auto_Line
#  --> likely represents data from clouds which were excluded during the manual compile
#
#there are 63 observations that are only present in the WT_Manual_Line
#  --> however, 30 of these observations are found in the WT_Auto_Line when the data is at 1dp rather than 4dp
#  --> this could be due to differences in the rounded values stored in the manual vs automated compile
#  --> by comparing the differences in the raw distance values only 4 datpoints differ by 0.01 while 26 differ by at least 1e-6
#  --> this suggests that these 26 data points are the same in the WT_Auto_Line and WT_Manual_Line
#  --> the other 4 observations have no obvious explanation why they are present in the WT_Manual_Line
#
#therefore, 5,776 observations are the same in the WT_Auto_Line and WT_Manual_Line (5,750 + 26)


################################################################################
#STEP 4:compare data present in manual vs auto compile

#shows no data points for each phase/pulse for manual line
WT_Manual_Data_Points_Per_Phase <- WT_Manual_Line %>%
  group_by(Phase) %>%
  mutate(No_Data_Points.manual = n()) %>%
  select(-Distance) %>% 
  distinct() %>%
  ungroup()

#shows no data points for each phase/pulse for auto line
WT_Auto_Data_Points_Per_Phase <- WT_Auto_Line %>%
  select(-Time, -Cloud_Name) %>% 
  group_by(Phase) %>%
  mutate(No_Data_Points.auto = n()) %>%
  ungroup() %>% 
  select(-Distance) %>% 
  distinct()

#compare number of data points present for each phase per pulse
Compare_Phase <- WT_Manual_Data_Points_Per_Phase %>%
  mutate(No_Data_Points.auto = WT_Auto_Data_Points_Per_Phase$No_Data_Points.auto,
         Difference_M_A = No_Data_Points.manual - No_Data_Points.auto)
  
#outcomes/observations
#shows that:
#  - Initiation has 267 more observations in WT_Auto_Line
#  - Maintenance has 11 more observations in WT_Auto_Line
#  --> together suggests that the automated compile does contain data which was excluded in the manual compile
#try to identify which data was removed

# #shows no. data points for each cloud for auto line
# WT_Auto_Data_Points_Per_Cloud <- WT_Auto_Line %>%
#   group_by(Data_Set, Phase, Pulse, Time, Cloud_Name) %>%
#   mutate(No_Data_Points = n()) %>%
#   ungroup() %>% 
#   select(-Distance, -Pulse) %>% 
#   distinct()
# 
# #show number of clouds per time point --> doesn't work give n = 2,000 etc
# WT_Auto_Clouds_Per_Time <- WT_Auto_Line %>%
#   select(-Distance, -Pulse) %>% 
#   group_by(Data_Set, Phase, Time) %>%
#   ungroup() %>% 
#   # select(-Cloud_Name) %>% 
#   # group_by(Data_Set, Phase, Time) %>%
#   mutate(No_Clouds = n()) %>% 
#   distinct()


################################################################################
#STEP 5: identify which clouds were excluded in the WT_Manual_Line compile
#there are 1,173 data points only found in the WT_Auto_Line
#  --> assume this difference is due to specific clouds being excluded from the WT_Manual_Line compile
#  --> identify which cloud data is only present in WT_Auto_Line and remove these clouds

#find cloud names and number of data points that are used in the WT_Auto_Line compile
name_Auto_Line <- WT_Auto_Line %>%
  group_by(Phase, Cloud_Name) %>%
  mutate(No_Data_Points = n()) %>%
  select(-Distance) %>%
  distinct() 


#find cloud names which are included in both WT_Auto_Line and WT_Manual_Line
name_same <- same_to_4dp %>%
  group_by(Phase,Time, Cloud_Name) %>%
  mutate(No_Data_Points = n()) %>%
  select(Cell_Line.x, Phase, Time, Cloud_Name, No_Data_Points) %>%
  distinct() %>%
  rename(Cell_Line = Cell_Line.x)

#identify which clouds are only found in name_only_in_auto - will show which clouds were not included in manual compile
excluded_clouds <- anti_join(name_Auto_Line,
                             name_same,
                             by= c("Cell_Line",
                                   "Phase",
                                   "Time",
                                   "Cloud_Name")) %>%
  ungroup() %>% 
  distinct()

#outcomes/observations
#only 234 clouds are used for the manual compile analysis --> 240 clouds have been included in the auto compile
#there are 6 clouds which are not present in the manual compile
#  --> Dynamic Initiation 50min 12_Cloud-1
#  --> Dynamic Initiation 50min 15_Cloud-1
#  --> Dynamic Initiation 50min 18_Cloud-1
#  --> Dynamic Initiation 60min 11_Cloud-2
#  --> Dynamic Initiation 60min 13_Cloud-1
#  --> Dynamic Maintenance 40min 08_Cloud-2
#
#investigate what happens when these cloulds are removed from the analysis

################################################################################
#STEP 6: does removing the clouds improve the similarity between manual and automated compile?

#load this data
setwd(r"(C:\Users\hmroa\Documents\RNA-SPLIT_Results\Manual_vs_Auto\Pulse_Chase_Analysis\mESCs\Nearest_Neighbour\WT_auto\nascent_Xist_dynamics)")

WT_Auto_CR_Initation <- read_csv("All_C2-C1_NNA_Data_Initiation_Clouds_Removed.csv") %>%
  transmute(Cell_Line = "WT_auto",
            Phase = "Initiation",
            Time = Time,
            Cloud_Name = Cloud_Name,
            Distance = NN_dist_2to1)

WT_Auto_CR_Maintenance <- read_csv("All_C2-C1_NNA_Data_Maintenance_Clouds_Removed.csv") %>%
  transmute(Cell_Line = "WT_auto",
            Phase = "Maintenance",
            Time = Time,
            Cloud_Name = Cloud_Name,
            Distance = NN_dist_2to1)

WT_Auto_Clouds_Removed <- WT_Auto_CR_Initation %>%
  bind_rows(WT_Auto_CR_Maintenance)

#outcomes/observations
#removing the previously identified clouds reduces the number of observations from 6,065 to 5,975 - WT_Manual_Line has 5,787 obersvations
#  --> there is no stat difference between steady state phase of the WT_Auto_Clouds_Removed and WT_Manual_Line
#  --> there are statistically higher distances in the expansion phase of the WT_Auto_Clouds_Removed
#      this difference is likely due to some data still being included in the auto compile which should be excluded


################################################################################
#STEP 8:#REPEAT STEPS 3-5 TO IDENTIFY WHAT IS DIFFERENT BETWEEN THE WT_Auto_Clouds_Removed and WT_Manual_Line

##################
#is the difference due to more observations in the expansion phase?

#shows no data points for each phase/pulse for manual line
WT_Manual_Data_Points_Per_Phase <- WT_Manual_Line %>%
  group_by(Phase) %>%
  mutate(No_Data_Points.manual = n()) %>%
  select(-Distance) %>% 
  distinct() %>%
  ungroup()

#shows no data points for each phase/pulse for auto line
WT_Auto_Data_Points_Per_Phase <- WT_Auto_Clouds_Removed %>%
  select(-Time, -Cloud_Name) %>% 
  group_by(Phase) %>%
  mutate(No_Data_Points.auto = n()) %>%
  ungroup() %>% 
  select(-Distance) %>% 
  distinct()

#compare number of data points present for each phase per pulse
Compare_Phase <- WT_Manual_Data_Points_Per_Phase %>%
  mutate(No_Data_Points.auto = WT_Auto_Data_Points_Per_Phase$No_Data_Points.auto,
         Difference_M_A = No_Data_Points.manual - No_Data_Points.auto)

#outcomes/observations
#shows WT_Auto_Clouds_Removed expansion phase has 181 more observations than the WT_Manual_Line
#shows WT_Auto_Clouds_Removed steady state phase has 7 more observations than the WT_Manual_Line
#  --> most likely explains why a stat dif is only observed between expansion phase (p=0.0075)


##################
#find which data is only present in the WT_Auto_Clouds_Removed

#finds which rows are the same in the manual and auto compile
same_to_4dp <- inner_join(mutate(WT_Auto_Clouds_Removed, Distance.auto = Distance,
                                 Distance.auto.4dp = round(Distance, digits = 3)),
                          mutate(WT_Manual_Line, Distance.manual = Distance,
                                 Distance.manual.4dp = round(Distance, digits = 3)),
                          by= c("Phase",
                                "Distance.auto.4dp" = "Distance.manual.4dp")) %>% 
  distinct(Phase,
           Time,
           Cloud_Name,
           Distance.auto, .keep_all=TRUE)

#finds which observations are only present in the automated compile
only_in_auto <- anti_join(mutate(WT_Auto_Clouds_Removed, Distance.auto = Distance,
                                 Distance.auto.4dp = round(Distance, digits = 4)),
                          mutate(WT_Manual_Line, Distance.manual = Distance,
                                 Distance.manual.4dp = round(Distance, digits = 4)),
                          by= c("Distance.auto.4dp" = "Distance.manual.4dp"))

#finds which observations are only present in the automated compile
only_in_manual <- anti_join(mutate(WT_Manual_Line, Distance.manual = Distance,
                                   Distance.manual.4dp = round(Distance, digits = 4)),
                            mutate(WT_Auto_Clouds_Removed, Distance.auto = Distance,
                                   Distance.auto.4dp = round(Distance, digits = 4)),
                            by= c("Distance.manual.4dp" = "Distance.auto.4dp"))

#find if there are any very similar datapoints between manual and auto compile
similar_to_1dp <- inner_join(mutate(only_in_auto, Distance.auto.1dp = round(Distance, digits = 1)),
                             mutate(only_in_manual, Distance.manual.1dp = round(Distance, digits = 1)),
                             by= c("Phase",
                                   "Distance.auto.1dp" = "Distance.manual.1dp")) %>% 
  distinct(Phase,
           Time,
           Cloud_Name,
           Distance.manual, .keep_all=TRUE) %>%
  select(-Distance.x, -Distance.auto.4dp, -Distance.y, -Distance.manual.4dp) %>%
  mutate(Difference = Distance.auto - Distance.manual) %>%    #calculates the differnce in original distance values (not rounded)
  mutate(Difference_squ = Difference*Difference)              #need to remove negatives

#outcomes/observations
#there are 5750 observations that are the same between the WT_Auto_Clouds_Removed and WT_Manual_Line (to 4dp)
#there are 63 observations only present in the WT_Manual_Line
#  --> out of these 21 observations are likely to be included in the WT_Auto_Clouds_Removed but were not identified as the same due to difference in rounding
#      this is because when matched to 1dp 21 additional observations were the same and they only differ by at least +/- 1e-6 
#  --> there is no obvious explanation for the other 42 differences


##################
#which clouds have been excluded from the WT_Manual_Line

#find cloud names and number of data points that are used in the WT_Auto_Line compile
name_Auto_Line <- WT_Auto_Clouds_Removed %>%
  group_by(Phase, Cloud_Name) %>%
  mutate(No_Data_Points = n()) %>%
  select(-Distance) %>%
  distinct() 


#find cloud names which are included in both WT_Auto_Line and WT_Manual_Line
name_same <- same_to_4dp %>%
  group_by(Phase,Time, Cloud_Name) %>%
  mutate(No_Data_Points = n()) %>%
  select(Cell_Line.x, Phase, Time, Cloud_Name, No_Data_Points) %>%
  distinct() %>%
  rename(Cell_Line = Cell_Line.x)

#identify which clouds are only found in name_only_in_auto - will show which clouds were not included in manual compile
excluded_clouds <- anti_join(name_Auto_Line,
                             name_same,
                             by= c("Cell_Line",
                                   "Phase",
                                   "Time",
                                   "Cloud_Name")) %>%
  ungroup() %>% 
  distinct()

#outcomes/observations
#234 clouds were used in the manual and auto compile
#the auto_compile includes 1 extra cloud as it used 235 clouds
#this cloud is 11_Cloud-2 Initiation at 60mins 
#  --> the inclusion of this cloud likely explains the differnce in expansion phase between WT_Auto_Clouds_Removed and WT_Manual_Line

################################################################################
#STEP 9: Does excluding 11_Cloud-2 Initiation at 60mins  improve the similarity between manual and automated compile?

setwd(r"(C:\Users\hmroa\Documents\RNA-SPLIT_Results\Manual_vs_Auto\Pulse_Chase_Analysis\mESCs\Nearest_Neighbour\WT_auto\nascent_Xist_dynamics)")

WT_Auto_CR_Initation_2 <- read_csv("All_C2-C1_NNA_Data_Initiation_Clouds_Removed_2.csv") %>%
  transmute(Cell_Line = "WT_auto",
            Phase = "Initiation",
            Time = Time,
            Cloud_Name = Cloud_Name,
            Distance = NN_dist_2to1)

WT_Auto_CR_Maintenance_2 <- read_csv("All_C2-C1_NNA_Data_Maintenance_Clouds_Removed_2.csv") %>%
  transmute(Cell_Line = "WT_auto",
            Phase = "Maintenance",
            Time = Time,
            Cloud_Name = Cloud_Name,
            Distance = NN_dist_2to1)

WT_Auto_Clouds_Removed_2 <- WT_Auto_CR_Initation_2 %>%
  bind_rows(WT_Auto_CR_Maintenance_2)

#outcomes/observations
#removing the previously identified 11_Cloud-2 Initiation at 60mins  reduces the number of observations from 5,975 ro 5,952 - WT_Manual_Line has 5,787 obersvations
#  --> there is no stat difference between steady state phase of the WT_Auto_Clouds_Removed_2 and WT_Manual_Line
#  --> there are statistically higher distances in the expansion phase of the WT_Auto_Clouds_Removed_2
#      this difference is likely due to some data still being included in the auto compile which should be excluded

###############
#shows no data points for each phase/pulse for manual line
WT_Manual_Data_Points_Per_Phase <- WT_Manual_Line %>%
  group_by(Phase) %>%
  mutate(No_Data_Points.manual = n()) %>%
  select(-Distance) %>% 
  distinct() %>%
  ungroup()

#shows no data points for each phase/pulse for auto line
WT_Auto_Data_Points_Per_Phase <- WT_Auto_Clouds_Removed_2 %>%
  select(-Time, -Cloud_Name) %>% 
  group_by(Phase) %>%
  mutate(No_Data_Points.auto = n()) %>%
  ungroup() %>% 
  select(-Distance) %>% 
  distinct()

#compare number of data points present for each phase per pulse
Compare_Phase <- WT_Manual_Data_Points_Per_Phase %>%
  mutate(No_Data_Points.auto = WT_Auto_Data_Points_Per_Phase$No_Data_Points.auto,
         Difference_M_A = No_Data_Points.manual - No_Data_Points.auto)

#outcomes/observations
#shows WT_Auto_Clouds_Removed_2 expansion phase has 158 more observations than the WT_Manual_Line
#shows WT_Auto_Clouds_Removed_2 steady state phase has 7 more observations than the WT_Manual_Line
#  --> most likely explains why a stat dif is only observed between expansion phase (p=0.0128 --> this is smaller than last time (p=0.0075))


##################
#find which data is only present in the WT_Auto_Clouds_Removed_2

#finds which rows are the same in the manual and auto compile
same_to_4dp <- inner_join(mutate(WT_Auto_Clouds_Removed_2, Distance.auto = Distance,
                                 Distance.auto.4dp = round(Distance, digits = 3)),
                          mutate(WT_Manual_Line, Distance.manual = Distance,
                                 Distance.manual.4dp = round(Distance, digits = 3)),
                          by= c("Phase",
                                "Distance.auto.4dp" = "Distance.manual.4dp")) %>% 
  distinct(Phase,
           Time,
           Cloud_Name,
           Distance.auto, .keep_all=TRUE)

#finds which observations are only present in the automated compile
only_in_auto <- anti_join(mutate(WT_Auto_Clouds_Removed_2, Distance.auto = Distance,
                                 Distance.auto.4dp = round(Distance, digits = 4)),
                          mutate(WT_Manual_Line, Distance.manual = Distance,
                                 Distance.manual.4dp = round(Distance, digits = 4)),
                          by= c("Distance.auto.4dp" = "Distance.manual.4dp"))

#finds which observations are only present in the automated compile
only_in_manual <- anti_join(mutate(WT_Manual_Line, Distance.manual = Distance,
                                   Distance.manual.4dp = round(Distance, digits = 4)),
                            mutate(WT_Auto_Clouds_Removed_2, Distance.auto = Distance,
                                   Distance.auto.4dp = round(Distance, digits = 4)),
                            by= c("Distance.manual.4dp" = "Distance.auto.4dp"))

#find if there are any very similar datapoints between manual and auto compile
similar_to_1dp <- inner_join(mutate(only_in_auto, Distance.auto.1dp = round(Distance, digits = 1)),
                             mutate(only_in_manual, Distance.manual.1dp = round(Distance, digits = 1)),
                             by= c("Phase",
                                   "Distance.auto.1dp" = "Distance.manual.1dp")) %>% 
  distinct(Phase,
           Time,
           Cloud_Name,
           Distance.manual, .keep_all=TRUE) %>%
  select(-Distance.x, -Distance.auto.4dp, -Distance.y, -Distance.manual.4dp) %>%
  mutate(Difference = Distance.auto - Distance.manual) %>%    #calculates the differnce in original distance values (not rounded)
  mutate(Difference_squ = Difference*Difference)

#find cloud names and number of data points that are used in the WT_Auto_Line compile
name_Auto_Line <- WT_Auto_Clouds_Removed_2 %>%
  group_by(Phase, Cloud_Name) %>%
  mutate(No_Data_Points = n()) %>%
  select(-Distance) %>%
  distinct() 


#find cloud names which are included in both WT_Auto_Line and WT_Manual_Line
name_same <- same_to_4dp %>%
  group_by(Phase,Time, Cloud_Name) %>%
  mutate(No_Data_Points = n()) %>%
  select(Cell_Line.x, Phase, Time, Cloud_Name, No_Data_Points) %>%
  distinct() %>%
  rename(Cell_Line = Cell_Line.x)

#identify which clouds are only found in name_only_in_auto - will show which clouds were not included in manual compile
excluded_clouds <- anti_join(name_Auto_Line,
                             name_same,
                             by= c("Cell_Line",
                                   "Phase",
                                   "Time",
                                   "Cloud_Name")) %>%
  ungroup() %>% 
  distinct()

#outcomes/observations
#the number of clouds used in the auto compile is the same as used in the manual compile
#  --> difference in WT_Auto_Clouds_Removed_2 is likely due to some rows of data not being included for certain clouds during manual compile
#  --> this difference is likely due to human error during manual compile
#
#the difference between the WT_Auto_Clouds_Removed_2 and WT_Manual_Line in the expansion phase is 158 observations
#there are 228 observations that are unique to the WT_Auto_Clouds_Removed_2 and 63 for the manual
#  --> 228 - 63 = 165, suggesting that the majority of the observation only in the manual line are actually included in the auto compile but haven't been due to difference in rounds
#will be hard to identify exactly which rows of data were missed during manual compile but all the data suggests that the automated compile is valid


################################################################################
#STEP 7: Plot manual vs auto centroids per pulse and phase

WT_Auto_Present <- WT_Auto_Clouds_Removed_2 %>%
  select(-Cloud_Name, -Time)

Present_Data <- WT_Manual_Line %>%
  bind_rows(WT_Auto_Present) %>%
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",                          #ensures naming consistency with other papers
                           Phase == "Maintenance" ~ "Steady_State")) %>%
  mutate(Phase2 = Phase) %>%
  mutate(Phase2 = case_when(Phase2 == "Expansion" ~ "Exp",                          
                                                 Phase2 == "Steady_State" ~ "SS")) %>%
  mutate(Key = paste(Phase2, Cell_Line, sep="_")) 

Present_Data$Cell_Line <- factor(Present_Data$Cell_Line, 
                                 levels = c("WT_manual", "WT_auto"), ordered = TRUE) 

#put in order for the graphs to appear on the x-axis
Present_Data$Key <- factor(Present_Data$Key, 
                           levels = c("Exp_WT_manual", "SS_WT_manual",
                                      "Exp_WT_auto", "SS_WT_auto"), ordered = TRUE)



NNA_plot <- ggplot(Present_Data, aes(x = Cell_Line, y = Distance, fill = Key)) +       #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                                   #defines type of plot
  facet_wrap(~Phase, nrow = 1) +                                                          #shows box plots for selected cell lines
  scale_fill_manual(values = c("#FC9272", "#EF3B2C",                #WT_manual exp and SS colours 
                               "#BDBDBD", "#737373" )) +            #WT_auto exp and SS colours 
  coord_cartesian(ylim = c((0),(670) )) +                                                 #adjusts scale so whiskers don't touch the end graph 
  labs(title = "C) Validate C2-C1 Compile",                                  #sets name of the axes
       x = "Cell_Line",
       y = "Nearest Neighbour Distances") +                                         
  theme                                                                                   #adds the theme which defines the font, text size etc


#view the box plot
NNA_plot



################################################################################
#STEP 8:preform t Test on Data to asses if there is a statistical difference
Refernece_Line <- "WT_manual"

#is there a general difference between the dynamic datasets?
Stat_Diff <- t.test(Present_Data[which(Present_Data$Cell_Line==Refernece_Line),]$Distance,
                         Present_Data[which(Present_Data$Cell_Line==New_Line_Name),]$Distance)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Present_Data[which(Present_Data$Cell_Line==Refernece_Line),]$Distance,
                         Present_Data[which(Present_Data$Cell_Line==New_Line_Name),]$Distance, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Present_Data[which(Present_Data$Cell_Line==Refernece_Line),]$Distance,
                         Present_Data[which(Present_Data$Cell_Line==New_Line_Name),]$Distance, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
dynamic_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Compare_Both_Phases", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_1 <- "Overall, the automated compile has statistically higher distance in automated than the manual compile"
  } else {
    outcome_1 <- "Overall, the new automated compile has statistically lower distance in Dynamic Dataset than the manual compile"
  }} else {
    outcome_1 <- "Overall, there is no statistical difference in distance in Dynamic Dataset between the automated compile and the manual compile"
  }

#is there a general difference between the pulse 1 data?
Expansion <- filter(Present_Data, Phase == "Expansion")
Stat_Diff <- t.test(Expansion[which(Expansion$Cell_Line==Refernece_Line),]$Distance,
                         Expansion[which(Expansion$Cell_Line==New_Line_Name),]$Distance)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Expansion[which(Expansion$Cell_Line==Refernece_Line),]$Distance,
                         Expansion[which(Expansion$Cell_Line==New_Line_Name),]$Distance, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Expansion[which(Expansion$Cell_Line==Refernece_Line),]$Distance,
                         Expansion[which(Expansion$Cell_Line==New_Line_Name),]$Distance, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
Expansion_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Compare_Expansion_WT", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_2 <- "Overall, the automated compile has statistically higher distance in Expansion than the manual compile"
  } else {
    outcome_2 <- "Overall, the new automated compile has statistically lower distance in Expansion than the manual compile"
  }} else {
    outcome_2 <- "Overall, there is no statistical difference in distance in Expansion between the automated compile and the manual compile"
  }

#is there a general difference between the pulse 1 data?
steady_State <- filter(Present_Data, Phase == "Steady_State")
Stat_Diff <- t.test(steady_State[which(steady_State$Cell_Line==Refernece_Line),]$Distance,
                         steady_State[which(steady_State$Cell_Line==New_Line_Name),]$Distance)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(steady_State[which(steady_State$Cell_Line==Refernece_Line),]$Distance,
                         steady_State[which(steady_State$Cell_Line==New_Line_Name),]$Distance, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(steady_State[which(steady_State$Cell_Line==Refernece_Line),]$Distance,
                         steady_State[which(steady_State$Cell_Line==New_Line_Name),]$Distance, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than th   WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
steady_State_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Compare_steady_State_WT", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_3 <- "Overall, the automated compile has statistically higher distance in steady_State than the manual compile"
  } else {
    outcome_3 <- "Overall, the new automated compile has statistically lower distance in steady_State than the manual compile"
  }} else {
    outcome_3 <- "Overall, there is no statistical difference in distance in steady_State between the automated compile and the manual compile"
  }


#save p values
p_values <- dynamic_p_values %>%
  bind_rows(Expansion_p_values,
            steady_State_p_values)

setwd(paste(File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, "nascent_Xist_dynamics", sep="\\"))

write_csv(p_values, paste(New_Line_Name, "Auto_vs_Manual_C2-C1_Compile_t-test_p_values_Clouds_Removed_2.csv", sep="_"))
write_csv(WT_Auto_Clouds_Removed_2, "Updated_WT_auto_NNA_2to1_Compile.csv")

#return outcome of the stats test
print(outcome_1)
print(outcome_2)
print(outcome_3)
