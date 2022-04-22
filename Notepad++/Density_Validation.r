################################################################################
#MANUAL VS AUTOMATIC C1-C1/C2-C2 COMPILE VALIDATION
################################################################################
#for WT data only use nascent dynamics data set to calculate density!
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
#STEP 1: upload all centroid data to pick a WT sample

setwd(paste(File_Path, Cell_Type, "Density", "All_Cell_Lines_Merged", sep="\\"))

#choose which cell lines to compare too
WT_Original <- read_csv("Cleaned_Other_Cell_Lines_Density_Compile.csv") %>%
  filter(Cell_Line == "WT") %>%
  mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_original"))

WT_Manual_Line <- read_csv("Manual_Compile_Density_Data.csv")

################################################################################
#STEP 2: upload WT automated compile to inspect - tibble includes cloud names

setwd(paste(File_Path, Cell_Type, "Density", New_Line_Name,sep="\\"))

#data relating to new cell line - also contains Cloud Names to allow identification of specific clouds 
WT_Auto_Line <- read_csv(paste(New_Line_Name, "Density_Compile.csv", sep="_")) %>%
  filter(Data_Set == "Dynamic")

#outcomes/observations
#the wT_Auto_Line contains 21,269 observations compared to 23,030 observation in the WT_Manual_Line (manual compiles)
#  --> as the automated has less observations it suggests that the manual compile contains data from clouds which weren't included in the C1/2-centroids compile
#  --> repeat automated compile including all cloud data

#load data including datapoints from previously excluded clouds
WT_Auto_Line_ALL <- read_csv(paste(New_Line_Name, "Density_Compile_ALL_DATA.csv", sep="_")) %>%
  filter(Data_Set == "Dynamic") 

#outcomes/observations
#this data still contains the less observations than the for the original/manual compile 
#  --> difference must be due to something else
#and the number observations is the same when the excluded clouds aren't used --> shows its valid to exclude the clouds
#could it be that the density data contains data from Turnover at 60min?


WT_Auto_Plus_T60 <- read_csv(paste(New_Line_Name, "Density_Compile.csv", sep="_")) %>%
  filter(Data_Set == "Dynamic" | Data_Set == "Turnover") %>%
  filter(Time %in%
           c("10", "20", "30", "40", "50", "60"))
  
#outcomes/observations
#this data still contains the more observations than the for the original/manual compile 
#  --> including this data has made the pulse 2 at steady-state more similar
#  --> however, including this data causes there to be s sig. df. between the manual/auto compile
#  --> this sig. dif. is probably due to data being included in the WT_Auto_Plus_T60 that was not included in the WT_Manual_Line
#  --> further analysis showed turnover 60min data was not included


#remove any potential duplications
WT_Manual_Line_deduped <- WT_Manual_Line %>% 
  distinct()

#outcomes/observations
#There are 23,030 observations in the WT_Manual_Line and 21,269 observations in the WT_Auto_Line 
#  --> suggests the WT_Manual_Line has an extra 1,761 observations
#  --> it is unclear where these extra observations come from
#  --> maybe these observations come from duplications in the manual compile
#removing duplicates still leaves 22,977 observations, which is still 1,708 more observations than the WT_Auto_Line
#  --> increase in observations is not due to duplication in the WT_Manual_Line 
#investigate what is different between the WT_Auto_Line and WT_Manual_Line


################################################################################
#STEP 3: What is the same and different about the manual (WT_Manual_Line) and automated (WT_Auto_Line) compiles?

#finds which rows are the same in the manual and auto compile
#  --> does create duplicates for when the same Distance ar recorded in the manual for the same, Time, Phase, Pulse and Dataset
same_to_4dp <- inner_join(mutate(WT_Auto_Line, Distance.auto = Distance,
                          Distance.auto.4dp = round(Distance, digits = 4)),
                   mutate(WT_Manual_Line, Distance.manual = Distance,
                          Distance.manual.4dp = round(Distance, digits = 4)),
                   by= c("Phase",
                         "Pulse",
                         "Distance.auto.4dp" = "Distance.manual.4dp")) %>% 
  distinct(Data_Set.x,
           Phase,
           Pulse,
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
                                  "Pulse",
                                  "Distance.auto.1dp" = "Distance.manual.1dp")) %>% 
                      distinct(Data_Set.x,
                               Phase,
                               Pulse,
                               Time,
                               Cloud_Name,
                               Distance.manual, .keep_all=TRUE) %>%
                      select(-Distance.x, -Distance.auto.4dp, -Distance.y, -Distance.manual.4dp) %>%
                      mutate(Difference = Distance.auto - Distance.manual) %>%    #calculates the differnce in original distance values (not rounded)
                      mutate(Difference_squ = Difference*Difference)              #need to remove negatives

#outcomes/observations
#out of the 23,030 observations in the WT_Manual_Line , 19,863 are the same as the WT_Auto_Line
#there are 1,403 observations that are only present in WT_Auto_Line
#  --> likely represents data from clouds which were excluded during the manual compile
#
#there are 928 observations that are only present in the WT_Manual_Line --> why only 928 and not 1,708??
#  --> however, 237 of these observations are found in the WT_Auto_Line when the data is at 1dp rather than 4dp
#  --> this could be due to differences in the rounded values stored in the manual vs automated compile
#  --> by comparing the differences in the raw distance values only 7 datpoints differ by 0.01 while 230 differ by at least 1e-6
#  --> this suggests that these 230 data points are the same in the WT_Auto_Line and WT_Manual_Line
#  --> the other 7 observations have no obvious explanation why they are present in the WT_Manual_Line
#
#therefore, 20,093 observations are the same in the WT_Auto_Line and WT_Manual_Line (19,863 + 230)
#  --> with only 1,173 data points unique in the WT_Auto_Line (1,403 - 230)
#      and 693 observations unique in the WT_Manual_Line (828 + 230)
#investigate which cloud data is only in the WT_Auto_Line


################################################################################
#STEP 3:compare data present in manual vs auto compile

#shows no data points for each phase/pulse for manual line
WT_Manual_Data_Points_Per_Phase <- WT_Manual_Line %>%
  group_by(Data_Set, Phase, Pulse) %>%
  mutate(No_Data_Points.manual = n()) %>%
  select(-Distance) %>% 
  distinct() %>%
  ungroup()

#shows no data points for each phase/pulse for auto line
WT_Auto_Data_Points_Per_Phase <- WT_Auto_Line %>%
  select(-Time, -Cloud_Name) %>% 
  group_by(Data_Set, Phase, Pulse) %>%
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
#  - Initiation pulse 1 only differs by 41 data points
#  - Initiation pulse 2 differs by 1,392 data points --> suggests automated compile is missing data
#  - Maintenance pulse 1 only differs by -275 data points 
#  - Maintenance pulse 2 differs by 603 data points --> suggests automated compile is missing data
#  --> together suggests that the automated compile does not contain all the data it is meant too
#try to identify which data it is missing

#shows no. data points for each cloud for auto line
WT_Auto_Data_Points_Per_Cloud <- WT_Auto_Line %>%
  group_by(Data_Set, Phase, Pulse, Time, Cloud_Name) %>%
  mutate(No_Data_Points = n()) %>%
  ungroup() %>% 
  select(-Distance, -Pulse) %>% 
  distinct()

#show number of clouds per time point --> doesn't work give n = 2,000 etc
WT_Auto_Clouds_Per_Time <- WT_Auto_Line %>%
  select(-Distance, -Pulse) %>% 
  group_by(Data_Set, Phase, Time) %>%
  ungroup() %>% 
  # select(-Cloud_Name) %>% 
  # group_by(Data_Set, Phase, Time) %>%
  mutate(No_Clouds = n()) %>% 
  distinct()


################################################################################
#STEP 4: identify which clouds were excluded in the WT_Manual_Line compile
#there are 1,173 data points only found in the WT_Auto_Line
#  --> assume this difference is due to specific clouds being excluded from the WT_Manual_Line compile
#  --> identify which cloud data is only present in WT_Auto_Line and remove these clouds

#find cloud names and number of data points that are used in the WT_Auto_Line compile
name_Auto_Line <- WT_Auto_Line %>%
  group_by(Data_Set, Phase, Pulse, Time, Cloud_Name) %>%
  mutate(No_Data_Points = n()) %>%
  select(-Distance) %>%
  distinct() 


#find cloud names which are included in both WT_Auto_Line and WT_Manual_Line
name_same <- same_to_4dp %>%
  group_by(Data_Set.x, Phase, Pulse, Time, Cloud_Name) %>%
  mutate(No_Data_Points = n()) %>%
  select(Cell_Line.x, Data_Set.x, Phase, Pulse, Time, Cloud_Name, No_Data_Points) %>%
  distinct() %>%
  rename(Cell_Line = Cell_Line.x,
         Data_Set = Data_Set.x)

#identify which clouds are only found in name_only_in_auto - will show which clouds were not included in manual compile
excluded_clouds <- anti_join(name_Auto_Line,
                             name_same,
                             by= c("Cell_Line",
                                   "Data_Set",
                                   "Phase",
                                   "Pulse",
                                   "Time",
                                   "Cloud_Name")) %>%
  ungroup() %>% 
  distinct()

#outcomes/observations
#only 267 clouds are used for the manual compile analysis --> 290 clouds have been included in the auto compile
#there are 23 clouds which are not present in the manual compile
#  --> 08_Cloud-1 in Dynamic Maintenance Pulse 1 at 40min
#  --> 10_Cloud-1 in Dynamic Maintenance Pulse 1 at 10min
#  --> the other 21 clouds are all from Dynamic Initiation Pulse 2 at 40min --> suggests this data has been missed during manual compile
#As no cloud data was missing in both the pulse 1 and pulse 2 it suggests that no additional clouds were excluded from the manual compile
#
#however, there there is still no explanation for why the WT_Manual_Line has an additional 1,761 observations
#  --> suggests the automated compile has data missing (most likely in pulse 2)
#  --> re-download the original raw data files for the WT and repeat the compile (maybe some files weren't downloaded)


#from inspecting what was only present in the manual and auto it shows the majority of differences are with pulse 2
#as Pulse 1 is near identical in the manual and automated compiles this validates the computational method
#  --> however, pulse 2 data does seem slightly different - could be due to manual having 1,700 more datapoints

################################################################################
#STEP 5: Plot manual vs auto centroids per pulse and phase

WT_Auto_Present <- WT_Auto_Line %>%
  select(-Cloud_Name, -Time)

Present_Data <- WT_Manual_Line %>%
  bind_rows(WT_Auto_Present) %>%
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",                          #ensures naming consistency with other papers
                           Phase == "Maintenance" ~ "Steady_State")) %>%
  mutate(Phase2 = Phase) %>%
  mutate(Phase2 = case_when(Phase2 == "Expansion" ~ "Exp",                          
                                                 Phase2 == "Steady_State" ~ "SS")) %>%
  mutate(Key = paste(Phase2, Cell_Line, Pulse, sep="_")) 

Present_Data$Cell_Line <- factor(Present_Data$Cell_Line, 
                                 levels = c("WT_manual", "WT_auto"), ordered = TRUE) 

#put in order for the graphs to appear on the x-axis
Present_Data$Key <- factor(Present_Data$Key, 
                           levels = c("Exp_WT_manual_Pulse_1", "Exp_WT_manual_Pulse_2", "SS_WT_manual_Pulse_1", "SS_WT_manual_Pulse_2",
                                      "Exp_WT_auto_Pulse_1", "Exp_WT_auto_Pulse_2", "SS_WT_auto_Pulse_1", "SS_WT_auto_Pulse_2"), ordered = TRUE)



Density_plot <- ggplot(Present_Data, aes(x = Phase, y = Distance, fill = Key)) +       #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                                   #defines type of plot
  facet_wrap(~Pulse, nrow = 1) +                                                          #shows box plots for selected cell lines
  scale_fill_manual(values = c("#FC9272", "#FC9272", "#EF3B2C", "#EF3B2C",                #WT_manual exp and SS colours 
                               "#BDBDBD", "#BDBDBD", "#737373", "#737373" )) +            #WT_auto exp and SS colours 
  coord_cartesian(ylim = c((365),(2350) )) +                                                 #adjusts scale so whiskers don't touch the end graph 
  labs(title = "B) Validate C1-C1 and C2-C2 Compile",                                  #sets name of the axes
       x = "Phase",
       y = "Median Distances Per Pulse") +                                         
  theme                                                                                   #adds the theme which defines the font, text size etc


#view the box plot
Density_plot



################################################################################
#STEP 7:preform t Test on Data to asses if there is a statistical difference
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
  transmute(Comparision = paste("DYNAMIC_Compare_Both_Pulses", New_Line_Name, sep='_'),
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
Pulse_1 <- filter(Present_Data, Pulse == "Pulse_1")
Stat_Diff <- t.test(Pulse_1[which(Pulse_1$Cell_Line==Refernece_Line),]$Distance,
                         Pulse_1[which(Pulse_1$Cell_Line==New_Line_Name),]$Distance)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Pulse_1[which(Pulse_1$Cell_Line==Refernece_Line),]$Distance,
                         Pulse_1[which(Pulse_1$Cell_Line==New_Line_Name),]$Distance, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Pulse_1[which(Pulse_1$Cell_Line==Refernece_Line),]$Distance,
                         Pulse_1[which(Pulse_1$Cell_Line==New_Line_Name),]$Distance, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
pulse_1_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Compare_Pulse_1", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_2 <- "The automated compile has statistically higher distance in pulse 1 than the manual compile"
  } else {
    outcome_2 <- "The new automated compile has statistically lower distance in pulse 1 than the manual compile"
  }} else {
    outcome_2 <- "There is no statistical difference in distance in pulse 1 between the automated compile and the manual compile"
  }

#is there a general difference between the pulse 1 data?
Pulse_2 <- filter(Present_Data, Pulse == "Pulse_2")
Stat_Diff <- t.test(Pulse_2[which(Pulse_2$Cell_Line==Refernece_Line),]$Distance,
                         Pulse_2[which(Pulse_2$Cell_Line==New_Line_Name),]$Distance)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Pulse_2[which(Pulse_2$Cell_Line==Refernece_Line),]$Distance,
                         Pulse_2[which(Pulse_2$Cell_Line==New_Line_Name),]$Distance, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Pulse_2[which(Pulse_2$Cell_Line==Refernece_Line),]$Distance,
                         Pulse_2[which(Pulse_2$Cell_Line==New_Line_Name),]$Distance, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
pulse_2_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Compare_Pulse_2", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_3 <- "The automated compile has statistically higher distance in pulse 2 than the manual compile"
  } else {
    outcome_3 <- "The new automated compile has statistically lower distance in pulse 2 than the manual compile"
  }} else {
    outcome_3 <- "There is no statistical difference in distance in pulse 2 between the automated compile and the manual compile"
  }


#save p values
p_values <- dynamic_p_values %>%
  bind_rows(pulse_1_p_values,
            pulse_2_p_values)

setwd(paste(File_Path, Cell_Type, "Density", New_Line_Name, sep="\\"))

write_csv(p_values, paste(New_Line_Name, "Auto_vs_Manual_C1-C1_C2-C2_Compile_t-test_p_values.csv", sep="_"))
write_csv(WT_Auto_Plus_T60, "Updated_WT_Auto")

#return outcome of the stats test
print(outcome_1)
print(outcome_2)
print(outcome_3)
