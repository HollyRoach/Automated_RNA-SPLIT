#############################################################################################################################################
#SCRIPT CREATES VIOLIN PLOT FOR NNA AND PREFORMS WILCOX STAT TEST
#please email Holly Roach at hmroach@hotmail.co.uk if you have any questions
#################################################################################

# #load libraries 
# library(tidyverse)
# library(RColorBrewer)
# library(stringr)
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
# 
# ################################################################################
# #USER INPUT REQUIRED - check lines 22-74
# 
# #define type of cells used in experiment
# Cell_Type <- "mESCs"                     #either "mESCs" or "NPCs"
# 
# #######################
# #SELECT DATA TO PRESENT
# 
# #set name of new cell line, and up to 3 others to present in the plot
# #all names should be spelled the same as file names within the directory
# #can leave names empty if don't want to present them all
# New_Line_Name <- "Test"       
# 
# Reference_Line_Name <- "WT"        #Reference_Line_Name will be used to perform stats test with

# #define which dataset you want to present 
# New_Line_Data <- "Turnover" #should be "Turnover" or "Dynamic"
# Ref_Line_Data <- "Dynamic" #should be "Turnover" or "Dynamic"
# 
# #define the threshold used to determine a significant result during wilcoxon test
# p <- 0.05

#######################
#DEFINE ASTHETICS OF PLOT

#may need to adjust scale for y-axis on line 

# #define colours for each cell line
# Colour_1_EdU <- "#FDE0DD"     #sets colour for EdU control of Reference_Line_Name
# Colour_1_exp <- "#FA9FB5"    #sets colour for expansion phase Reference_Line_Name
# Colour_1_ss <- "#DD3497"     #sets colour for steady state phase Reference_Line_Name
# Colour_1_Random <- "#FCC5C0" #sets colour for Random control of Reference_Line_Name
# 
# Colour_2_EdU <- "#ECE2F0"     #sets colour for EdU control of new cell line
# Colour_2_exp <- "#A6BDDB"    #sets colour for expansion phase new cell line
# Colour_2_ss <- "#3690C0"     #sets colour for steady state phase new cell line
# Colour_2_Random <- "#D0D1E6" #sets colour for Random control of new cell line



# #set name of how cell lines should be presented in the plot - ensures name consistency with other papers
# #if name contains delta/triangle symbol use - "SPEN^"~Delta*"RRM"
# Name_1 <- bquote("WT")                   #sets name for Reference_Line_Name
# Name_2 <- bquote("METTL3_FKBP12"^"F36V") #sets name for new cell line

#set name of plot titles
Title <- "Assessment of Xist RNP Coupling"
x_axis <- ""
y_axis <- "Distance to Nearest Neighbour [nm]"

#################################################################################
#STEP 1: load NNA data for plotting

#define file path to where "Pulse_Chase_Analysis" is located - this is where the compiled data is stored
File_Path <- Output_File_Path

#set file path to the location of the all Density files 
Input_Path <- paste(File_Path, Cell_Type, "Nearest_Neighbour", "All_Cell_Lines", sep="/")

#stores name of all Nearest_Neighbour_Compile.csv files , in the directory, into a vector
Files <- list.files(path = Input_Path, pattern = "Nearest_Neighbour_Compile.csv", full.names  = TRUE)   

#creates a list containing the Nearest_Neighbour_Compile.csv files
File_List <- lapply(Files, read_csv, col_types="cccncn")

#concatenates all the individual Nearest_Neighbour_Compile.csv files into 1 tibble
All_NNA_Data <- bind_rows(File_List)


#checks if new cell line data is stored in the dataframe containing all NNA data
if (!any(All_NNA_Data$Cell_Line == New_Line_Name)) {
  stop("This new cell line does not exist in Main Dataframe (All_NNA_Data)
       - make sure to run Cloud_Volume_Compilation_Manipulation.R ")
}

#################################################################################
#STEP 2: SELECT WHICH DATA TO PRESENT IN THE PLOT
#ultimately will present 2 sets of violin plots

#select new cell line data to present
Present_New <- All_NNA_Data %>% 
  filter(Cell_Line == New_Line_Name) %>% 
  filter(Data_Set == New_Line_Data)

#select new cell line data to present
Present_Ref <- All_NNA_Data %>% 
  filter(Cell_Line == Reference_Line_Name) %>% 
  filter(Data_Set == Ref_Line_Data)

#combine ref and new line data in order to present together
Present_Data <- bind_rows(Present_New, Present_Ref) %>%               
  mutate(Phase = case_when(Phase == "EdU" ~ "EdU",                        #ensures naming consistency with other papers                      
                           Phase == "Initiation" ~ "Expansion",
                           Phase == "Maintenance" ~ "Steady_State",
                           Phase == "random_sample" ~ "Random",)) %>%
  mutate(Key = paste(Phase, Cell_Line, sep="_"))                          #creates variable to base colour coding in plot


#order data to be present in violin plot
Present_Data$Key <- factor(Present_Data$Key, 
                             levels = c(paste("EdU", Reference_Line_Name, sep="_"),
                                        paste("Expansion", Reference_Line_Name, sep="_"),
                                        paste("Steady_State", Reference_Line_Name, sep="_"),
                                        paste("Random", Reference_Line_Name, sep="_"), 
                                        paste("EdU", New_Line_Name, sep="_"),
                                        paste("Expansion", New_Line_Name, sep="_"),
                                        paste("Steady_State", New_Line_Name, sep="_"),
                                        paste("Random", New_Line_Name, sep="_"), ordered = TRUE)) 

Present_Data$Cell_Line <- factor(Present_Data$Cell_Line, 
                                 levels = c(Reference_Line_Name, New_Line_Name), ordered = TRUE,
                                 labels = c(Name_1, Name_2))

Present_Data$Phase <- factor(Present_Data$Phase, 
                           levels = c("EdU", "Expansion", "Steady_State", "Random"), ordered = TRUE) 


#######################################################################################
#STEP 3: plot violin plots which compare the NNA Data for the chosen cell line

#sets name of labels for each cell line
Labels = c(Cell_Line = Name_1,
           Cell_Line = Name_2)

#create violin plot
NNA_plot <- ggplot(Present_Data, aes(x = Phase, y = Distance, fill = Key)) +    #assigns data to be plotted, sets colour based on cell line
  geom_violin(aes(fill = Key)) +                                                #sets data to be plotted as a violin plot - displays median as a line (draw_quantiles = c(0.5)
  facet_wrap(~Cell_Line, nrow = 1, labeller = label_parsed) +                   #separates cell lines, and sets labels
  coord_cartesian(ylim=c(0, 1500)) +                                            #adjusts scale of the y axis 
  scale_fill_manual(values = c(Colour_1_EdU, Colour_1_exp,                      #defines the colour of each boxplot
                               Colour_1_ss, Colour_1_Random,                        
                               Colour_2_EdU, Colour_2_exp,
                               Colour_2_ss, Colour_2_Random)) +                  
  geom_boxplot(width=0.1, outlier.shape = NA, aes(fill = Key)) +                #adds box plot over the violin plots
  labs(title = Title,                                                           #sets name of the axes
       x = x_axis,
       y = y_axis) 

#view violin plot
print(NNA_plot)

################################################################################
#STEP 4: need to manually export the plots


################################################################################
#STEP 5: preform Wilcox Test on Data to asses if there is a statistical difference
#can only directly compare results to internal control as NNA distances are largely effected by internal measurements
#if coupling has occurred:
#  - there will be significantly higher distances than the technical (EdU) control
#  - there will be significantly lower distances than the random control
################################################################################

#selects data for new cell line to calculate p values
Wilcox <- filter(All_NNA_Data, Cell_Line == New_Line_Name)

###################################
#IS THERE A STATISTICAL DIFFERENCE BETWEEN THE Expansion AND INTERNAL EdU CONTROL?

#is there a general difference in distances recorded between the data sets?
Edu_Exp_Stat_Diff <- wilcox.test(Wilcox[which(Wilcox$Phase=="EdU"),]$Distance,
                         Wilcox[which(Wilcox$Phase=="Initiation"),]$Distance)
Edu_Exp_Stat_Diff$p.value #p = 1 suggests a statistical difference / p nearly 0 suggests no statistical difference


#if there is a difference in distance, is it statistically less? 
Edu_Exp_Diff_Less <- wilcox.test(Wilcox[which(Wilcox$Phase=="EdU"),]$Distance,
                 Wilcox[which(Wilcox$Phase=="Initiation"),]$Distance, alternative = "less")
Edu_Exp_Diff_Less$p.value #p = 1 suggests new cell line lower distances than EdU 

#if there is a difference in distance, is it statistically Greater?
Edu_Exp_Diff_Greater <- wilcox.test(Wilcox[which(Wilcox$Phase=="EdU"),]$Distance,
                 Wilcox[which(Wilcox$Phase=="Initiation"),]$Distance, alternative = "greater")
Edu_Exp_Diff_Greater$p.value #p = 1 suggests new cell line has higher distances than EdU

#record  EdU p values in table
Edu_Exp_Diff <- Edu_Exp_Stat_Diff$p.value
Edu_Exp_Less <- Edu_Exp_Diff_Less$p.value
Edu_Exp_Greater <- Edu_Exp_Diff_Greater$p.value

Edu_Exp_p_values <- tibble(Edu_Exp_Diff, Edu_Exp_Less, Edu_Exp_Greater) %>%
  transmute(Comaprision = "EdU_vs_Exp_Data",
            General_Stat_Diff = Edu_Exp_Diff,
            Diff_Statistically_Less = Edu_Exp_Less,
            Diff_Statistically_Greater = Edu_Exp_Greater) 

#says whether there is a stat difference between the new data and EdU Control
if (Edu_Exp_Stat_Diff$p.value < p) {
  if (Edu_Exp_Diff_Less$p.value < Edu_Exp_Diff_Greater$p.value){
    outcome_1 <- "The new cell line has statistically higher distances in Expansion phase than the EdU Control"
    outcome_1a <- TRUE
  }  else {
    outcome_1 <- "The new cell line has statistically lower distances in Expansion phase than the EdU Control"
    outcome_1a <- FALSE
  }}else {
    outcome_1 <- "There is no statistical difference in distances between the new cell line in Expansion phase and the EdU Control"
    outcome_1a <- FALSE
  }

###################################
#IS THERE A STATISTICAL DIFFERENCE BETWEEN THE Expansion DATA AND INTERNAL RANDOM CONTROL?

#is there a general difference in distances recorded between the data sets?
Random_Exp_Stat_Diff <- wilcox.test(Wilcox[which(Wilcox$Phase=="random_sample"),]$Distance,
                             Wilcox[which(Wilcox$Phase=="Initiation"),]$Distance)
Random_Exp_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference in distance, is it statistically less? 
Random_Exp_Diff_Less <- wilcox.test(Wilcox[which(Wilcox$Phase=="random_sample"),]$Distance,
                             Wilcox[which(Wilcox$Phase=="Initiation"),]$Distance, alternative = "less")
Random_Exp_Diff_Less$p.value #p = 1 suggests new cell line lower distances than Random 

#if there is a difference in distance, is it statistically Greater?
Random_Exp_Diff_Greater <- wilcox.test(Wilcox[which(Wilcox$Phase=="random_sample"),]$Distance,
                             Wilcox[which(Wilcox$Phase=="Initiation"),]$Distance, alternative = "greater")
Random_Exp_Diff_Greater$p.value #p = 1 suggests new cell line has higher distances than Random

#record  Random p values in table
Random_Exp_Diff <- Random_Exp_Stat_Diff$p.value
Random_Exp_Less <- Random_Exp_Diff_Less$p.value
Random_Exp_Greater <- Random_Exp_Diff_Greater$p.value

Random_Exp_p_values <- tibble(Random_Exp_Diff, Random_Exp_Less, Random_Exp_Greater) %>%
  transmute(Comaprision = "Random_vs_Expansion_Data",
            General_Stat_Diff = Random_Exp_Diff,
            Diff_Statistically_Less = Random_Exp_Less,
            Diff_Statistically_Greater = Random_Exp_Greater) 

#says whether there is a stat difference between the new expansion data and Random Control
if (Random_Exp_Stat_Diff$p.value < p) {
  if (Random_Exp_Diff_Less$p.value < Random_Exp_Diff_Greater$p.value){
    outcome_2 <- "The new cell line has statistically higher distances in Expansion phase than the Random Control"
    outcome_2a <- FALSE
  }  else {
    outcome_2 <- "The new cell line has statistically lower distances in Expansion phase than the Random Control"
    outcome_2a <- TRUE
  }}else {
    outcome_2 <- "There is no statistical difference in distances between the new cell line in Expansion phase and the Random Control"
    outcome_2a <- FALSE
  }

###################################
#IS THERE A STATISTICAL DIFFERENCE BETWEEN THE Steady_State AND INTERNAL EdU CONTROL?

#is there a general difference in distances recorded between the data sets?
Edu_SS_Stat_Diff <- wilcox.test(Wilcox[which(Wilcox$Phase=="EdU"),]$Distance,
                                 Wilcox[which(Wilcox$Phase=="Maintenance"),]$Distance)
Edu_SS_Stat_Diff$p.value #p = 1 suggests a statistical difference / p nearly 0 suggests no statistical difference


#if there is a difference in distance, is it statistically less? 
Edu_SS_Diff_Less <- wilcox.test(Wilcox[which(Wilcox$Phase=="EdU"),]$Distance,
                                 Wilcox[which(Wilcox$Phase=="Maintenance"),]$Distance, alternative = "less")
Edu_SS_Diff_Less$p.value #p = 1 suggests new cell line lower distances than EdU 

#if there is a difference in distance, is it statistically Greater?
Edu_SS_Diff_Greater <- wilcox.test(Wilcox[which(Wilcox$Phase=="EdU"),]$Distance,
                                    Wilcox[which(Wilcox$Phase=="Maintenance"),]$Distance, alternative = "greater")
Edu_SS_Diff_Greater$p.value #p = 1 suggests new cell line has higher distances than EdU

#record  EdU p values in table
Edu_SS_Diff <- Edu_SS_Stat_Diff$p.value
Edu_SS_Less <- Edu_SS_Diff_Less$p.value
Edu_SS_Greater <- Edu_SS_Diff_Greater$p.value

Edu_SS_p_values <- tibble(Edu_SS_Diff, Edu_SS_Less, Edu_SS_Greater) %>%
  transmute(Comaprision = "EdU_vs_SS_Data",
            General_Stat_Diff = Edu_SS_Diff,
            Diff_Statistically_Less = Edu_SS_Less,
            Diff_Statistically_Greater = Edu_SS_Greater) 

#says whether there is a stat difference between the new data and EdU Control
if (Edu_SS_Stat_Diff$p.value < p) {
  if (Edu_SS_Diff_Less$p.value < Edu_SS_Diff_Greater$p.value){
    outcome_3 <- "The new cell line has statistically higher distances in Steady_State phase than the EdU Control"
    outcome_3a <- TRUE
  }  else {
    outcome_3 <- "The new cell line has statistically lower distances in Steady_State phase than the EdU Control"
    outcome_3a <- FALSE
  }}else {
    outcome_3 <- "There is no statistical difference in distances between the new cell line in Steady_State phase and the EdU Control"
    outcome_3a <- FALSE
  }

###################################
#IS THERE A STATISTICAL DIFFERENCE BETWEEN THE Steady_State DATA AND INTERNAL RANDOM CONTROL?

#is there a general difference in distances recorded between the data sets?
Random_SS_Stat_Diff <- wilcox.test(Wilcox[which(Wilcox$Phase=="random_sample"),]$Distance,
                                    Wilcox[which(Wilcox$Phase=="Maintenance"),]$Distance)
Random_SS_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference in distance, is it statistically less? 
Random_SS_Diff_Less <- wilcox.test(Wilcox[which(Wilcox$Phase=="random_sample"),]$Distance,
                                    Wilcox[which(Wilcox$Phase=="Maintenance"),]$Distance, alternative = "less")
Random_SS_Diff_Less$p.value #p = 1 suggests new cell line lower distances than Random 

#if there is a difference in distance, is it statistically Greater?
Random_SS_Diff_Greater <- wilcox.test(Wilcox[which(Wilcox$Phase=="random_sample"),]$Distance,
                                       Wilcox[which(Wilcox$Phase=="Maintenance"),]$Distance, alternative = "greater")
Random_SS_Diff_Greater$p.value #p = 1 suggests new cell line has higher distances than Random

#record  Random p values in table
Random_SS_Diff <- Random_SS_Stat_Diff$p.value
Random_SS_Less <- Random_SS_Diff_Less$p.value
Random_SS_Greater <- Random_SS_Diff_Greater$p.value

Random_SS_p_values <- tibble(Random_SS_Diff, Random_SS_Less, Random_SS_Greater) %>%
  transmute(Comaprision = "Random_vs_Steady_State_Data",
            General_Stat_Diff = Random_SS_Diff,
            Diff_Statistically_Less = Random_SS_Less,
            Diff_Statistically_Greater = Random_SS_Greater) 

#says whether there is a stat difference between the new Steady_State data and Random Control
if (Random_SS_Stat_Diff$p.value < p) {
  if (Random_SS_Diff_Less$p.value < Random_SS_Diff_Greater$p.value){
    outcome_4 <- "The new cell line has statistically higher distances in Steady_State phase than the Random Control"
    outcome_4a <- FALSE
  }  else {
    outcome_4 <- "The new cell line has statistically lower distances in Steady_State phase than the Random Control"
    outcome_4a <- TRUE
  }}else {
    outcome_4 <- "There is no statistical difference in distances between the new cell line in Steady_State phase and the Random Control"
    outcome_4a <- FALSE
  }


###################################
#IS THERE A STATISTICAL DIFFERENCE BETWEEN THE STEADY-STATE AND Steady_State PHASES?

#is there a general difference in distances recorded between the data sets?
Phase_Stat_Diff <- wilcox.test(Wilcox[which(Wilcox$Phase=="Maintenance"),]$Distance,
                                Wilcox[which(Wilcox$Phase=="Initiation"),]$Distance)
Phase_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference in distance, is it statistically less? 
Phase_Diff_Less <- wilcox.test(Wilcox[which(Wilcox$Phase=="Maintenance"),]$Distance,
                                Wilcox[which(Wilcox$Phase=="Initiation"),]$Distance, alternative = "less")
Phase_Diff_Less$p.value #p = 1 suggests new cell line lower distances than Phase 

#if there is a difference in distance, is it statistically Greater?
Phase_Diff_Greater <- wilcox.test(Wilcox[which(Wilcox$Phase=="Maintenance"),]$Distance,
                                Wilcox[which(Wilcox$Phase=="Initiation"),]$Distance, alternative = "greater")
Phase_Diff_Greater$p.value #p = 1 suggests new cell line has higher distances than Phase

#record  Phase p values in table
Phase_Diff <- Phase_Stat_Diff$p.value
Phase_Less <- Phase_Diff_Less$p.value
Phase_Greater <- Phase_Diff_Greater$p.value

Phase_p_values <- tibble(Phase_Diff, Phase_Less, Phase_Greater) %>%
  transmute(Comaprision = "Steady_State_vs_Expansion_Data",
            General_Stat_Diff = Phase_Diff,
            Diff_Statistically_Less = Phase_Less,
            Diff_Statistically_Greater = Phase_Greater) 


#says whether there is a stat difference between the expansion and steady state phases
if (Phase_Stat_Diff$p.value < p) {
  if (Phase_Diff_Less$p.value < Phase_Diff_Greater$p.value){
    outcome_5 <- "The Expansion Phase has statistically higher distances than the Steady-State"
  }  else {
    outcome_5 <- "The Expansion Phase has statistically lower distances than the Steady-State"
  }}else {
    outcome_5 <- "There is no statistical difference in distances between the Expansion Phase and Steady-State"
  }

#record all p values for each test and save the results
NNA_p_Values <- Edu_Exp_p_values %>%
  bind_rows(Edu_SS_p_values,
            Random_Exp_p_values,
            Random_SS_p_values,
            Phase_p_values)

#save p values
Save_Path <- paste(File_Path, Cell_Type, "Nearest_Neighbour", New_Line_Name, sep="/")

write_csv(NNA_p_Values, paste(Save_Path, paste(New_Line_Name, "NNA_Wilcox_Test_p_values.csv", sep=""), sep="/"))


###################################
#return outcomes of the Wilcox Test

#overall states whether coupling has been observed in the new cell line
if(outcome_1a & outcome_2a) {
  outcome_6 <- "Suggests that the new cell line has conversved coupling in the expansion phase"
} else {
  outcome_6 <- "Suggests that the new cell line disrupts couplingin the expansion phase"
}


if(outcome_3a & outcome_4a) {
  outcome_7 <- "Suggests that the new cell line has conversved coupling in the steady-state phase"
} else {
  outcome_7 <- "Suggests that the new cell line disrupts couplingin the steady-state phase"
}

print(outcome_5)
print(outcome_6)
print(outcome_7)