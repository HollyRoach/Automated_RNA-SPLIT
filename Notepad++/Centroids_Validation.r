################################################################################
#MANUAL VS AUTOMATIC C1/C2-centroids COMPILE VALIDATION
################################################################################
#
#load libraries
library(extrafont)
font_import() #need to do once at beginning before the rest of the script can run
loadfonts(device = "win") #loads windows fonts to use
library(tidyverse)
library(stats)
library(RColorBrewer)
library(stringr)

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

setwd(paste(File_Path, Cell_Type, "Turnover", "All_Cell_Lines_Merged", sep="\\"))

#choose which cell lines to compare too
WT_Manual_Line <- read_csv("Original_WT_Centroids_Compile_Includes_Turnover_Pulse_2.csv")

################################################################################
#STEP 2: upload WT automated compile to inspect - tibble includes cloud names

setwd(paste(File_Path, Cell_Type, "Turnover", New_Line_Name,sep="\\"))

#data relating to new cell line - also contains Cloud Names to allow identification of specific clouds 
WT_Auto_Line <- read_csv(paste(New_Line_Name, "Centroids_Compile_ALL_DATA.csv", sep="_"))

# #remove Turnover Pulse 2 data as manual compile does not contain this data
# No_Turnover_P2_Auto <- WT_Auto_Line %>%
#   filter(Data_Set == "Dynamic" | Data_Set == "Turnover" & Pulse == "Pulse_1")
# 
# Only_Turnover_P2_Auto <- WT_Auto_Line %>%
#   filter(Data_Set == "Turnover" & Pulse == "Pulse_2")

#outcomes/observations
#the wT_Auto_Line contains 1342 observations compared to 1336 observation in the WT_Manual_Line (manual compiles)

################################################################################
#STEP 3: What is the same and different about the manual (WT_Manual_Line) and automated (WT_Auto_Line) compiles?

#finds which rows are the same in the manual and auto compile
#  --> does create duplicates for when the same no_centroids ar recorded in the manual for the same, Time, Phase, Pulse and Dataset
same <- inner_join(mutate(WT_Auto_Line, No_Centroids.auto = No_Centroids),
                   mutate(WT_Manual_Line, No_Centroids.manual = No_Centroids),
                   by= c("Data_Set",
                         "Phase",
                         "Pulse",
                         "Time",
                         "No_Centroids" = "No_Centroids"))

#creates table of the duplicates created by inner_join - leaves only rows which have the same observations in manual and auto
same_dups <- same[duplicated(same, by=c("Data_Set",
                             "Phase",
                             "Pulse",
                             "Time",
                             "Cloud_Name",)), ]

#removes the duplicates created by inner_join - leaves only rows which have the same observations in manual and auto
same_deduped <- same  %>% 
  distinct(Data_Set,
           Phase,
           Pulse,
           Time,
           Cloud_Name, .keep_all=TRUE)

#finds which observations are only present in the automated compile
only_in_auto <- anti_join(WT_Auto_Line, WT_Manual_Line, by= c("Data_Set",
                                               "Phase",
                                               "Pulse",
                                               "Time",
                                               "No_Centroids" = "No_Centroids"))

#finds which observations are only present in the automated compile
only_in_manual <- anti_join(WT_Manual_Line, WT_Auto_Line, by= c("Data_Set",
                                                "Phase",
                                                "Pulse",
                                                "Time",
                                                "No_Centroids")) %>%
  transmute(Cell_Line = Cell_Line,                                    #puts in same order to allow easier comparison between tables
            Data_Set = Data_Set,
            Phase = Phase,
            Pulse = Pulse,
            Time = Time,
            No_Centroids = No_Centroids,
            Normalised_No_Centroids = Normalised_No_Centroids)

#outcomes/observations
#out of the 1342 observations in the WT_Auto_Line, 1329 observations are the same/found in the WT_Manual_Line 
#there are 13 observations that are only present in the WT_Auto_Line
#there are 8 observations that are only present in the WT_Manual_Line


################################################################################
#STEP 4: what doe these difference mean and how to proceed

# Difference is only +/- 1 centroid for 3 clouds therefore are likely to be referring to the same cloud - difference could be due to human error during manual compile
# -	Dynamic Initiation Pulse 1 at 10min = 161 vs 162 (01_Cloud-1)
# -	Dynamic Initiation Pulse 1 at 40min = 218 vs 217 (05_Cloud-1)
# -	Turnover Maintenance Pulse 1 at 0min = 4 vs 3 (10_Cloud-1)
# ---> will keep data in both manual and automated compiles

# Data only present in automated compile suggests these clouds were removed during manual compile - represents 10 datapoints, corresponding to 8 different clouds
# -	Dynamic Initiation Pulse 1 at 50min = 0 (03_Cloud-1)
# -	Turnover Initiation Pulse 1 at 120min = 265 (20_Cloud-1)
# -	Turnover Initiation Pulse 2 at 120min = 30 (05_Cloud-2)
# -	Turnover Initiation Pulse 2 at 120min = 30 (20_Cloud-1)
# -	Turnover Maintenance Pulse 1 at 140min = 30 (04_Cloud-1)
# -	Turnover Maintenance Pulse 1 at 200min = 312 (19_Cloud-1) 
# -	Turnover Maintenance Pulse 1 at 200min = 503 (20_Cloud-1)
# -	Turnover Maintenance Pulse 2 at 200min = 13 (04_Cloud-1)
# -	Turnover Maintenance Pulse 2 at 200min = 198 (19_Cloud-1)
# -	Turnover Maintenance Pulse 2 at 200min = 60 (20_Cloud-1)
# ---> only remove clouds which were present in pulse 1 and 2 as this is a strong indicator they were removed before manual compile 
#      and not used in any following data analysis
#      Turnover Initiation Pulse 1 at 120min = 265 (20_Cloud-1)
#      Turnover Initiation Pulse 2 at 120min = 30 (20_Cloud-1) --> remove watershed algorithm results for this cloud
#      Turnover Maintenance Pulse 1 at 200min = 312 (19_Cloud-1)
#      Turnover Maintenance Pulse 2 at 200min = 198 (19_Cloud-1) --> remove watershed algorithm results for this cloud
#      Turnover Maintenance Pulse 1 at 200min = 503 (20_Cloud-1)
#      Turnover Maintenance Pulse 2 at 200min = 60 (20_Cloud-1) --> remove watershed algorithm results for this cloud
# ---> clouds which only appear in either pulse 1 or pulse 2 might have been missed during manual compiling rather than being purposely excluded from data analysis 

# However, there are still different 5 observations in the manual compile which cannot be accounted for 
# -	Dynamic Maintenance Pulse 1 at 50min = 0 doesn't exist in automated compile
# -	Turnover Initiation Pulse 2 at 100min = 0 doesn't exist in automated compile
# -	Turnover Maintenance Pulse 1 at 140min = 2 doesn't exist in automated compile
# -	Turnover Maintenance Pulse 1 at 140min = 1 doesn't exist in automated compile
# -	Turnover Maintenance Pulse 2 at 200min = 32 doesn't exist in automated compile
# ---> will keep these datapoints in the manual compile as there is no evidence to suggest they need removing



################################################################################
#STEP 5: check whether the datasets become more similar by removing cloud data from the WT_Auto_Line that wasn't in the WT_Manual_Line


#uses auto compiled data which doesn't include data for the specified clouds 
WT_Auto_Clouds_Removed <- read_csv("Updated_WT_Auto_Centroid_Compile.csv")


#check the differences between the 2 data sets
update_same <- inner_join(mutate(WT_Auto_Clouds_Removed, No_Centroids.auto = No_Centroids),
                   mutate(WT_Manual_Line, No_Centroids.manual = No_Centroids),
                   by= c("Data_Set",
                         "Phase",
                         "Pulse",
                         "Time",
                         "No_Centroids" = "No_Centroids"))

update_same_dups <- same[duplicated(update_same, by=c("Data_Set",
                                                      "Phase",
                                                      "Pulse",
                                                      "Time",
                                                      "Cloud_Name",)), ]

update_same_deduped <- update_same  %>% 
  distinct(Data_Set,
           Phase,
           Pulse,
           Time,
           Cloud_Name, .keep_all=TRUE)


update_only_in_auto <- anti_join(WT_Auto_Clouds_Removed, WT_Manual_Line, by= c("Data_Set",
                                                                              "Phase",
                                                                              "Pulse",
                                                                              "Time",
                                                                              "No_Centroids" = "No_Centroids"))


update_only_in_manual <- anti_join(WT_Manual_Line, WT_Auto_Clouds_Removed, by= c("Data_Set",
                                                                                 "Phase",
                                                                                 "Pulse",
                                                                                 "Time",
                                                                                 "No_Centroids")) %>%
  transmute(Cell_Line = Cell_Line,                                    #puts in same order to allow easier comparison between tables
            Data_Set = Data_Set,
            Phase = Phase,
            Pulse = Pulse,
            Time = Time,
            No_Centroids = No_Centroids,
            Normalised_No_Centroids = Normalised_No_Centroids)


#outcomes/observations
#out of the 1328 observations in the WT_Auto_Line, 1325 observations are the same/found in the WT_Manual_Line 
#there are 3 observations that are only present in the WT_Auto_Line (update_only_in_auto)
#  ---> these observations equate to the 3 previously identified clouds which only have differ by +/- 1 centroid in the manual vs automated compile
#there are 9 observations that are only present in the WT_Manual_Line
#  ---> 3 of these equate to the 3 previously identified clouds which only have differ by +/- 1 centroid in the manual vs automated compile
#  ---> 5 of these correspond to the data which was previously identified in the manual compile 
#      (these 5 difference are still unaccounted for but were expected to be present in the update_only_in_manual table)
#        -	Dynamic Maintenance Pulse 1 at 50min = 0 doesn't exist in automated compile
#        -	Turnover Initiation Pulse 2 at 100min = 0 doesn't exist in automated compile
#        -	Turnover Maintenance Pulse 1 at 140min = 2 doesn't exist in automated compile
#        -	Turnover Maintenance Pulse 1 at 140min = 1 doesn't exist in automated compile
#        -	Turnover Maintenance Pulse 2 at 200min = 32 doesn't exist in automated compile
#  ---> this leaves only 1 unexpected difference: Turnover Maintenance Pulse 2 at 140min = 2 
#       it it possible that this could be due to the removal of 04_Cloud-1 from the	Turnover Maintenance at 140min in the automated pipeline
#       which was previously identified for removal - it is possible that by mistake the pulse 2 data was left in while the pulse 1 data was removed during the manual compile
#as there are no was increase in differences in the update_only_in_manual table it suggests that removing the specific cloud data from the auto compile is justified

#OVERALL THERE ARE ONLY 5 (out of >1300 observations) UNACCOUNTABLE DIFFERENCES BETWEEN AUTO AND MANUAL COMPILE 
# --> this VALIDATES this automated compile for the C1/C2-centroids.csv files worked 


################################################################################
#STEP 6: Plot manual vs auto centroids per pulse and phase

WT_Auto_Present <- WT_Auto_Clouds_Removed %>%
  select(-Cloud_Name)

Present_Data <- WT_Manual_Line %>%
  mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_manual")) %>%
  bind_rows(WT_Auto_Present) %>%
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",                          #ensures naming consistency with other papers
                           Phase == "Maintenance" ~ "Steady_State")) %>%
  mutate(Phase2 = Phase) %>%
  mutate(Phase2 = case_when(Phase2 == "Expansion" ~ "Exp",                          
                                                 Phase2 == "Steady_State" ~ "SS"))


#plot how the number of centroids in DYANIC DATASET
Present_Dynamic <- filter(Present_Data, Data_Set == "Dynamic") %>%
  mutate(Key = paste(Phase2, Cell_Line, Pulse, sep="_")) 

Present_Dynamic$Cell_Line <- factor(Present_Dynamic$Cell_Line, 
                                 levels = c("WT_manual", "WT_auto"), ordered = TRUE) #use names defined in line 51

#put in order for the graphs to appear on the x-axis
Present_Dynamic$Key <- factor(Present_Dynamic$Key, 
                           levels = c("Exp_WT_manual_Pulse_1", "Exp_WT_manual_Pulse_2", "SS_WT_manual_Pulse_1", "SS_WT_manual_Pulse_2",
                                      "Exp_WT_auto_Pulse_1", "Exp_WT_auto_Pulse_2", "SS_WT_auto_Pulse_1", "SS_WT_auto_Pulse_2"), ordered = TRUE)



Dynamic_plot <- ggplot(Present_Dynamic, aes(x = Phase, y = No_Centroids, fill = Key)) +       #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                                   #defines type of plot
  facet_wrap(~Pulse, nrow = 1) +                                                          #shows box plots for selected cell lines
  scale_fill_manual(values = c("#FC9272", "#FC9272", "#EF3B2C", "#EF3B2C",                #WT_manual exp and SS colours 
                               "#BDBDBD", "#BDBDBD", "#737373", "#737373" )) +            #WT_auto exp and SS colours 
  coord_cartesian(ylim = c((0),(230) )) +                                                 #adjusts scale so whiskers don't touch the end graph 
  labs(title = "A) Validate C1 and C2-Centroids Compile for Dynamic Dataset",                                  #sets name of the axes
       x = "Phase",
       y = "No. Centroids Per Pulse") +                                         
  theme                                                                                   #adds the theme which defines the font, text size etc


#view the box plot
Dynamic_plot


#plot how the number of centroids in TURNOVER DATASET
Present_Turnover <- filter(Present_Data, Data_Set == "Turnover") %>%
  mutate(Key = paste(Phase2, Cell_Line, Pulse, sep="_")) 

Present_Turnover$Cell_Line <- factor(Present_Turnover$Cell_Line, 
                                    levels = c("WT_manual", "WT_auto"), ordered = TRUE) #use names defined in line 51

#put in order for the graphs to appear on the x-axis
Present_Turnover$Key <- factor(Present_Turnover$Key, 
                              levels = c("Exp_WT_manual_Pulse_1", "Exp_WT_manual_Pulse_2", "SS_WT_manual_Pulse_1", "SS_WT_manual_Pulse_2",
                                         "Exp_WT_auto_Pulse_1", "Exp_WT_auto_Pulse_2", "SS_WT_auto_Pulse_1", "SS_WT_auto_Pulse_2"), ordered = TRUE)


Turnover_plot <- ggplot(Present_Turnover,aes(x = Phase, y = No_Centroids, fill = Key)) +  #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                                   #defines type of plot
  facet_wrap(~Pulse, nrow = 1) +                                                          #shows box plots for selected cell lines
  scale_fill_manual(values = c("#FC9272", "#FC9272", "#EF3B2C", "#EF3B2C",                #WT_manual exp and SS colours 
                               "#BDBDBD", "#BDBDBD", "#737373", "#737373" )) +            #WT_auto exp and SS colours 
  coord_cartesian(ylim = c((0),(155) )) +                                                 #adjusts scale so whiskers don't touch the end graph 
  labs(title = "B) Validate C1 and C2-Centroids Compile for Turnover Dataset",            #sets name of the axes
       x = "Phase",
       y = "No. Centroids Per Pulse") +                                         
  theme                                                                                #adds the theme which defines the font, text size etc


#view the box plot
Turnover_plot

################################################################################
#STEP 7:preform Wilcox Test on Data to asses if there is a statistical difference
Refernece_Line <- "WT_manual"

#is there a general difference between the dynamic datasets?
Stat_Diff <- t.test(Present_Dynamic[which(Present_Dynamic$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Dynamic[which(Present_Dynamic$Cell_Line==New_Line_Name),]$No_Centroids)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Present_Dynamic[which(Present_Dynamic$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Dynamic[which(Present_Dynamic$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Present_Dynamic[which(Present_Dynamic$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Dynamic[which(Present_Dynamic$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
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
    outcome_1 <- "Overall, the automated compile has statistically higher centroid count in Dynamic Dataset than the manual compile"
  } else {
    outcome_1 <- "Overall, the new automated compile has statistically lower centroid count in Dynamic Dataset than the manual compile"
  }} else {
    outcome_1 <- "Overall, there is no statistical difference in centroid count in Dynamic Dataset between the automated compile and the manual compile"
  }

#is there a general difference between the pulse 1 data?
Pulse_1 <- filter(Present_Dynamic, Pulse == "Pulse_1")
Stat_Diff <- t.test(Pulse_1[which(Pulse_1$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_1[which(Pulse_1$Cell_Line==New_Line_Name),]$No_Centroids)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Pulse_1[which(Pulse_1$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_1[which(Pulse_1$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Pulse_1[which(Pulse_1$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_1[which(Pulse_1$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
D_pulse_1_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("DYNAMIC_Compare_Pulse_1_WT", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_2 <- "The automated compile has statistically higher No_Centroids in pulse 1, Turnover data, than the manual compile"
  } else {
    outcome_2 <- "The new automated compile has statistically lower No_Centroids in pulse 1, Turnover data, than the manual compile"
  }} else {
    outcome_2 <- "There is no statistical difference in No_Centroids in pulse 1, Turnover data, between the automated compile and the manual compile"
  }

#is there a general difference between the pulse 1 data?
Pulse_2 <- filter(Present_Dynamic, Pulse == "Pulse_2")
Stat_Diff <- t.test(Pulse_2[which(Pulse_2$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_2[which(Pulse_2$Cell_Line==New_Line_Name),]$No_Centroids)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Pulse_2[which(Pulse_2$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_2[which(Pulse_2$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Pulse_2[which(Pulse_2$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_2[which(Pulse_2$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
D_pulse_2_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("DYNAMIC_Compare_Pulse_2_WT", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_3 <- "The automated compile has statistically higher No_Centroids in pulse 2, dyanmic data, than the manual compile"
  } else {
    outcome_3 <- "The new automated compile has statistically lower No_Centroids in pulse 2, dyanmic data, than the manual compile"
  }} else {
    outcome_3 <- "There is no statistical difference in No_Centroids in pulse 2, dyanmic data, between the automated compile and the manual compile"
  }


#is there a general difference between the turnover datasets?
Stat_Diff <- t.test(Present_Turnover[which(Present_Turnover$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Turnover[which(Present_Turnover$Cell_Line==New_Line_Name),]$No_Centroids)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Present_Turnover[which(Present_Turnover$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Turnover[which(Present_Turnover$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Present_Turnover[which(Present_Turnover$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Turnover[which(Present_Turnover$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
turnover_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("TURNOVER_Compare_Both_Pulses", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_4 <- "Overall, the automated compile has statistically higher centroid count in Turnover Dataset than the manual compile"
  } else {
    outcome_4 <- "Overall, the new automated compile has statistically lower centroid count in Turnover Dataset  than the manual compile"
  }} else {
    outcome_4 <- "Overall, there is no statistical difference in centroid count in Turnover Dataset between the automated compile and the manual compile"
  }

#is there a general difference between the pulse 1 data?
Pulse_1 <- filter(Present_Turnover, Pulse == "Pulse_1")
Stat_Diff <- t.test(Pulse_1[which(Pulse_1$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_1[which(Pulse_1$Cell_Line==New_Line_Name),]$No_Centroids)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Pulse_1[which(Pulse_1$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_1[which(Pulse_1$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Pulse_1[which(Pulse_1$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_1[which(Pulse_1$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
T_pulse_1_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Turnover_Compare_Pulse_1_WT", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_5 <- "The automated compile has statistically higher No_Centroids in pulse 1, dyanmic data, than the manual compile"
  } else {
    outcome_5 <- "The new automated compile has statistically lower No_Centroids in pulse 1, dyanmic data, than the manual compile"
  }} else {
    outcome_5 <- "There is no statistical difference in No_Centroids in pulse 1, dyanmic data, between the automated compile and the manual compile"
  }

#is there a general difference between the pulse 1 data?
Pulse_2 <- filter(Present_Turnover, Pulse == "Pulse_2")
Stat_Diff <- t.test(Pulse_2[which(Pulse_2$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_2[which(Pulse_2$Cell_Line==New_Line_Name),]$No_Centroids)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Pulse_2[which(Pulse_2$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_2[which(Pulse_2$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Pulse_2[which(Pulse_2$Cell_Line==Refernece_Line),]$No_Centroids,
                         Pulse_2[which(Pulse_2$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
T_pulse_2_p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Turnover_Compare_Pulse_2_WT", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_6 <- "The automated compile has statistically higher No_Centroids in pulse 2, Turnover data, than the manual compile"
  } else {
    outcome_6 <- "The new automated compile has statistically lower No_Centroids in pulse 2, Turnover data, than the manual compile"
  }} else {
    outcome_6 <- "There is no statistical difference in No_Centroids in pulse 2, Turnover data, between the automated compile and the manual compile"
  }

#save p values
p_values <- dynamic_p_values %>%
  bind_rows(turnover_p_values,
            D_pulse_1_p_values,
            D_pulse_2_p_values,
            T_pulse_1_p_values,
            T_pulse_2_p_values)

setwd(paste(File_Path, Cell_Type, "Turnover", New_Line_Name, sep="\\"))

write_csv(p_values, paste(New_Line_Name, "Auto_vs_Manual_Count_t-test_p_values.csv", sep="_"))

#retrun outcome of the stats test
print(outcome_1)
print(outcome_2)
print(outcome_3)
print(outcome_4)
print(outcome_5)
print(outcome_6)


################################################################################
#TOTAL MOLECULE cOUNT ANALSYIS 
#calculate total number of centroids per cloud - WT manual compile
#separate into pulse 1 and pulse 2 data

#calculate total molecules for dynamic manual WT
Manual_D_Pulse_1 <- WT_Manual %>%
  filter(Data_Set == "Dynamic" & Pulse == "Pulse_1") %>%
  rename(P1_No_Centroids = No_Centroids)

Manual_D_Pulse_2 <- WT_Manual %>%
  filter(Data_Set == "Dynamic" & Pulse == "Pulse_2") %>%
  rename(P2_No_Centroids = No_Centroids)

#puts pulse 1 and pulse 2 columns in the same tibble --> uses these coulmns to find total count
Manual_D_Total_Count <- Manual_D_Pulse_1 %>%
  mutate(P2_No_Centroids = Manual_D_Pulse_2$P2_No_Centroids,
         Total_Centroids = Manual_D_Pulse_1$P1_No_Centroids + Manual_D_Pulse_2$P2_No_Centroids) %>%
  select(-Normalised_No_Centroids)

#calculate total molecules for Turnover manual WT
Manual_T_Pulse_1 <- WT_Manual %>%
  filter(Data_Set == "Turnover" & Pulse == "Pulse_1") %>%
  rename(P1_No_Centroids = No_Centroids)

Manual_T_Pulse_2 <- WT_Manual %>%
   filter(Data_Set == "Turnover" & Pulse == "Pulse_2") %>%
   rename(P2_No_Centroids = No_Centroids)

#puts pulse 1 and pulse 2 columns in the same tibble --> uses these coulmns to find total count
Manual_T_Total_Count <- Manual_T_Pulse_1 %>%
   mutate(P2_No_Centroids = Manual_T_Pulse_2$P2_No_Centroids,
          Total_Centroids = Manual_T_Pulse_1$P1_No_Centroids + Manual_T_Pulse_2$P2_No_Centroids) %>%
   select(-Normalised_No_Centroids)

#combine dynamic and turnover data for manual compile
Manual_Total_Count <- Manual_D_Total_Count %>%
  bind_rows(Manual_T_Total_Count) %>%
  mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_manual"))

################################################################################
#STEP b: TOTAL MOLECULE COUNT ANALSYIS 
##calculate total number of centroids per cloud - WT auto compile
#separate into pulse 1 and pulse 2 data

#calculate total molecules for dynamic Auto WT
Auto_D_Pulse_1 <- WT_Auto %>%
  filter(Data_Set == "Dynamic" & Pulse == "Pulse_1") %>%
  rename(P1_No_Centroids = No_Centroids)

Auto_D_Pulse_2 <- WT_Auto %>%
  filter(Data_Set == "Dynamic" & Pulse == "Pulse_2") %>%
  rename(P2_No_Centroids = No_Centroids)

#puts pulse 1 and pulse 2 columns in the same tibble --> uses these coulmns to find total count
Auto_D_Total_Count <- Auto_D_Pulse_1 %>%
  mutate(P2_No_Centroids = Auto_D_Pulse_2$P2_No_Centroids,
         Total_Centroids = Auto_D_Pulse_1$P1_No_Centroids + Auto_D_Pulse_2$P2_No_Centroids) %>%
  select(-Normalised_No_Centroids)

#calculate total molecules for Turnover Auto WT
Auto_T_Pulse_1 <- WT_Auto %>%
  filter(Data_Set == "Turnover" & Pulse == "Pulse_1") %>%
  rename(P1_No_Centroids = No_Centroids)

#cannot calculate as Auto compile does not contain pulse 2 turnover data 
Auto_T_Pulse_2 <- WT_Auto %>%
   filter(Data_Set == "Turnover" & Pulse == "Pulse_2") %>%
   rename(P2_No_Centroids = No_Centroids)

#puts pulse 1 and pulse 2 columns in the same tibble --> uses these coulmns to find total count
 Auto_T_Total_Count <- Auto_T_Pulse_1 %>%
   mutate(P2_No_Centroids = Auto_T_Pulse_2$P2_No_Centroids,
          Total_Centroids = Auto_T_Pulse_1$P1_No_Centroids + Auto_T_Pulse_2$P2_No_Centroids) %>%
   select(-Normalised_No_Centroids)

#combine dynamic and turnover data for Auto compile
Auto_Total_Count <- Auto_D_Total_Count %>%
  bind_rows(Auto_T_Total_Count)


################################################################################
#STEP c: TOTAL MOLECULE COUNT ANALSYIS 

# #plot original vs manual vs auto compile
# setwd(paste(File_Path, Cell_Type, "Molecule_Count", "All_Cell_Lines_Merged", sep="\\"))
# Original_Total_Count <- read_csv("Cleaned_Other_Cell_Lines_Molecule_Count_Compile.csv") %>%
#   filter(Cell_Line == "WT") %>%
#   mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_original"))
  
  
#combine auto and manual data to present
Present_Count <- Auto_Total_Count %>%
  select(-Cloud_Name) %>%
  bind_rows(Manual_Total_Count) %>% 
  select(Cell_Line, Phase, Total_Centroids) %>%
  #bind_rows(Original_Total_Count) %>%
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",                          #ensures naming consistency with other papers
                           Phase == "Maintenance" ~ "Steady_State")) %>%                 
  mutate(Key = paste(Phase, Cell_Line, sep="_"))                                         #creates key variable to base colour coding off


#order data to be presented in box plot 
Present_Count$Cell_Line <- factor(Present_Count$Cell_Line, 
                                 levels = c(#"WT_original",
                                            "WT_manual", "WT_auto"), ordered = TRUE) #use names defined in line 51

#put in order for the graphs to appear on the x-axis
Present_Count$Key <- factor(Present_Count$Key, 
                           levels = c(#"Expansion_WT_original", "Steady_State_WT_original",
                                      "Expansion_WT_manual", "Steady_State_WT_manual",
                                      "Expansion_WT_auto", "Steady_State_WT_auto"), ordered = TRUE)


#create theme for plots - defines font, text size etc
theme <- theme(plot.title = element_text(family = "Calibri", face = "bold", size = (13)),
               legend.title = element_text(family = "Calibri", size = (11)), 
               legend.text = element_text(family = "Calibri", size = (11)), 
               axis.title = element_text(family = "Calibri", face = "bold", size = (11)),
               axis.text = element_text(family = "Calibri", face = "bold", size = (9)))


#make box plot
count_plot <- ggplot(Present_Count, aes(x = Phase, y = Total_Centroids, fill = Key)) +       #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                                   #defines type of plot
  facet_wrap(~Cell_Line, nrow = 1) +                                                      #shows box plots for selected cell lines
  scale_fill_manual(values = c(#"#FA9FB5", "#DD3497",                                      #WT_original exp and SS colours 
                               "#FC9272", "#EF3B2C",                                      #WT_manual exp and SS colours 
                               "#BDBDBD", "#737373" )) +                                  #WT_auto exp and SS colours 
  coord_cartesian(ylim = c((0),(290) )) +                                             #adjusts scale so whiskers don't touch the end graph 
  labs(title = "A) Validate Total Molecule Count Analysis",                                                     #sets name of the axes
       x = "Phase",
       y = "Total Xist Count [no. centroids]") +                                         
  theme                                                                                   #adds the theme which defines the font, text size etc


#view the box plot
count_plot

################################################################################
#STEP d: TOTAL MOLECULE COUNT ANALSYIS 
#preform t-test Test on Data to asses if there is a statistical difference
Refernece_Line <- "WT_manual"

#is there a general difference between the dynamic datasets?
Stat_Diff <- t.test(Present_Count[which(Present_Count$Cell_Line==Refernece_Line),]$Total_Centroids,
                         Present_Count[which(Present_Count$Cell_Line==New_Line_Name),]$Total_Centroids)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- t.test(Present_Count[which(Present_Count$Cell_Line==Refernece_Line),]$Total_Centroids,
                         Present_Count[which(Present_Count$Cell_Line==New_Line_Name),]$Total_Centroids, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- t.test(Present_Count[which(Present_Count$Cell_Line==Refernece_Line),]$Total_Centroids,
                         Present_Count[which(Present_Count$Cell_Line==New_Line_Name),]$Total_Centroids, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("TOTAL_Compare_Both_Phases_WT_vs", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < 0.05) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_1 <- "Overall, the automated compile has statistically higher centroid count than the manual compile"
  } else {
    outcome_1 <- "Overall, the new automated compile has statistically lower centroid count than the manual compile"
  }} else {
    outcome_1 <- "Overall, there is no statistical difference in centroid count between the automated compile and the manual compile"
  }


#save p values

setwd(paste(File_Path, Cell_Type, "Molecule_Count", New_Line_Name, sep="\\"))

write_csv(p_values, paste(New_Line_Name, "Auto_vs_Manual_Count_t-test_p_values.csv", sep="_"))

#return outcome of the stats test
print(outcome_1)


################################################################################
#TRANSCRIPTION DYNAMICS ANALSYIS 
#clean Transcription Dynamics
#CHECK - I think Lisa has included 0min and 60min from turnover data - need to repeat when I get access to the turnover pulse 2

Manual_Transcription <- WT_Manual %>%
  filter(Data_Set == "Dynamic" | Data_Set == "Turnover") %>% 
  filter(Pulse == "Pulse_2") %>%                             #only want to include pulse 2 (this relates to newly synthesised Xist RNPs)
  filter(Time %in%                                           #selects only time 20, 40 and 60min to present
           c("10", "20", "30", "40", "50", "60")) %>%
  mutate(Time = case_when(Time == 10 ~ 20,                    #new data needs to grouped by times - 20, 40, 60
                          Time == 20 ~ 20,
                          Time == 30 ~ 40,
                          Time == 40 ~ 40,
                          Time == 50 ~ 60,
                          Time == 60 ~ 60)) %>%
  select(Cell_Line = Cell_Line,                              #only need to keep cell line, phase, no_centroid columns
         Phase = Phase,
         Time = Time,
         No_Centroids = No_Centroids)

Auto_Transcription <- WT_Auto %>%
  filter(Data_Set == "Dynamic" | Data_Set == "Turnover") %>% 
  filter(Pulse == "Pulse_2") %>%                             #only want to include pulse 2 (this relates to newly synthesised Xist RNPs)
  filter(Time %in%                                           #selects only time 20, 40 and 60min to present
           c( "10", "20", "30", "40", "50", "60")) %>%
  mutate(Time = case_when(Time == 10 ~ 20,                   #new data needs to grouped by times - 20, 40, 60
                          Time == 20 ~ 20,
                          Time == 30 ~ 40,
                          Time == 40 ~ 40,
                          Time == 50 ~ 60,
                          Time == 60 ~ 60)) %>%
  select(Cell_Line = Cell_Line,                              #only need to keep cell line, phase, no_centroid columns
         Phase = Phase,
         Time = Time,
         No_Centroids = No_Centroids)

# #load original transcription dynamics data
# setwd(paste(File_Path, Cell_Type, "Transcription_Dynamics", "All_Cell_Lines_Merged", sep="\\"))
# Original_Transcription <- read_csv("Cleaned_Other_Cell_Lines_Transcription_Dynamics_Compile.csv") %>%
#   mutate(Cell_Line = "WT_original") %>%
#   filter(Time %in%                                           #selects only time 20, 40 and 60min to present
#            c( "0", "10", "20", "30", "40", "50", "60"))

#load automated data before specific datapoints are removed
setwd(paste(File_Path, Cell_Type, "Turnover", New_Line_Name, sep="\\"))

Auto_Original <- read_csv("WT_auto_Centroids_Compile_ALL_DATA.csv") %>%
  filter(Data_Set == "Dynamic" | Data_Set == "Turnover") %>% #only want to include nascent_Xist_dynamics data
  filter(Pulse == "Pulse_2") %>%                             #only want to include pulse 2 (this relates to newly synthesised Xist RNPs)
  filter(Time %in%                                           #selects only time 20, 40 and 60min to present
           c( #"0",
              "10", "20", "30", "40", "50", "60")) %>%
  mutate(Time = case_when(Time == 0 ~ 20,                    #new data needs to grouped by times - 20, 40, 60
                          Time == 10 ~ 20,       
                          Time == 20 ~ 20,
                          Time == 30 ~ 40,
                          Time == 40 ~ 40,
                          Time == 50 ~ 60,
                          Time == 60 ~ 60)) %>%
  select(Cell_Line = Cell_Line,                              #only need to keep cell line, phase, no_centroid columns
         Phase = Phase,
         Time = Time,
         No_Centroids = No_Centroids)

################################################################################
#STEP b: TRANSCRIPTION DYNAMICS ANALSYIS 
#plot manual vs auto compile

Present_Transcription  <- Manual_Transcription %>%
  mutate(Cell_Line = case_when(Cell_Line == "WT" ~ "WT_manual"))    %>%
  bind_rows(Auto_Original,
            #Original_Transcription
            ) %>%
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",
                           Phase == "Maintenance" ~ "Steady_State")) %>%
  mutate(Key = paste(Phase, Cell_Line, sep="_")) %>%                                  #creates variable to base colour coding in plot
  filter(Time %in%                                                                    #selects only time 20, 40 and 60min to present
           c("20","40", "60"))

#order data to be presented in box plot 
Present_Transcription$Cell_Line <- factor(Present_Transcription$Cell_Line, 
                                    levels = c(#"WT_original",
                                               "WT_manual", "WT_auto"), ordered = TRUE) 

#put in order for the graphs to appear on the x-axis
Present_Transcription$Key <- factor(Present_Transcription$Key, 
                                levels = c(#"Expansion_WT_original", "Steady_State_WT_original",
                                           "Expansion_WT_manual", "Steady_State_WT_manual",
                                           "Expansion_WT_auto", "Steady_State_WT_auto"), ordered = TRUE)

Present_Transcription$Time  <- factor(Present_Transcription$Time,
                             levels = c("20","40", "60"), ordered = TRUE)

#create theme for plots - defines font, text size etc
theme <- theme(plot.title = element_text(family = "Calibri", face = "bold", size = (13)),
               legend.title = element_text(family = "Calibri", size = (11)), 
               legend.text = element_text(family = "Calibri", size = (11)), 
               axis.title = element_text(family = "Calibri", face = "bold", size = (11)),
               axis.text = element_text(family = "Calibri", face = "bold", size = (9)))

#calculate the y scale for the plot
Min_y = min(Present_Transcription$No_Centroids, na.rm = TRUE)
Max_y = max(Present_Transcription$No_Centroids, na.rm = TRUE)  #adjust value according to graph output as doesn't contain outliers

#make box plot 
transcription_plot <- ggplot(Present_Transcription, aes(x = Time, y = No_Centroids, fill = Key)) + #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                         #defines type of plot
  facet_wrap(~Cell_Line, nrow = 1) +                                                #shows box plots seperated by phase
  scale_fill_manual(values = c(#"#FA9FB5", "#DD3497",                                      #WT_original exp and SS colours 
                               "#FC9272", "#EF3B2C",                                      #WT_manual exp and SS colours 
                               "#BDBDBD", "#737373" )) +                                  #WT_auto exp and SS colours 
  coord_cartesian(ylim = c(Min_y, 90)) +                                        #adjusts so whiskers don't touch the end graph 
  labs(title = "C) Validate Transcription Dynamics Analysis",                              #sets name of the axes
       x = "Time [min]",
       y = "Newly Synthesized Xist Count [no. centroids]") +           
  theme                                                                         #adds the theme which defines the font, text size etc

#view the box plot
transcription_plot

################################################################################
#c: TRANSCRIPTION DYNAMICS ANALSYIS 
#preform t-test Test on Data to asses if there is a statistical difference
Refernece_Line <- "WT_manual"

#is there a general difference between the dynamic datasets?
Stat_Diff_T <- t.test(Present_Transcription[which(Present_Transcription$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Transcription[which(Present_Transcription$Cell_Line==New_Line_Name),]$No_Centroids)
Stat_Diff_T $p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less_T  <- t.test(Present_Transcription[which(Present_Transcription$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Transcription[which(Present_Transcription$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Diff_Less_T $p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More_T  <- t.test(Present_Transcription[which(Present_Transcription$Cell_Line==Refernece_Line),]$No_Centroids,
                         Present_Transcription[which(Present_Transcription$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
Diff_More_T $p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
d <- Stat_Diff_T $p.value
e <- Diff_Less_T $p.value
f <- Diff_More_T $p.value
p_values <- tibble(d, e, f) %>%
  transmute(Comparision = paste("TOTAL_Compare_Both_Phases", New_Line_Name, sep='_'),
            General_Stat_Diff = d,
            Diff_Statistically_Less = e,
            Diff_Statistically_More = f) 

if (Stat_Diff_T$p.value < 0.05) {
  if (Diff_Less_T$p.value < Diff_More_T$p.value){
    outcome_2 <- "Overall, the automated compile has statistically higher transcription rate than the manual compile"
  } else {
    outcome_2 <- "Overall, the new automated compile has statistically lower transcription rate than the manual compile"
  }} else {
    outcome_2 <- "Overall, there is no statistical difference in transcription rate between the automated compile and the manual compile"
  }


#save p values

setwd(paste(File_Path, Cell_Type, "Transcription_Dynamics", New_Line_Name, sep="\\"))

write_csv(p_values, paste(New_Line_Name, "Auto_vs_Manual_Transcription_Dynamics_t-test_p_values.csv", sep="_"))

#return outcome of the stats test
print(outcome_1)
print(outcome_2)