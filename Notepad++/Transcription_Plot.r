#############################################################################################################################################
#SCRIPT CREATES BOX PLOT FOR Transcription Dynamics DATA AND PREFORMS WILCOX STAT TEST
#############################################################################################################################################
#load libraries 

#run these two lines on their own first
library(extrafont)
#font_import() 

#finish running rest of the script from this line
loadfonts(device = "win") #loads windows fonts to use
library(tidyverse)
library(RColorBrewer)
library(stringr)


################################################################################
#USER INPUT REQUIRED - check lines 21-71

#define type of cells used in experiment
Cell_Type <- "mESCs"                     #either "mESCs" or "NPCs"

#######################
#SELECT DATA TO PRESENT

#set name of new cell line, and up to 3 others to present in the plot
#all names should be spelled the same as file names within the directory
#can leave names empty if don't want to present them all
New_Line_Name <- "Mettl3_dTAG"       

Refernece_Line <- "WT"        #reference_line will be used to perform stats test with

Cell_Line_3 <- "Ciz1_KO"

Cell_Line_4 <- "SPEN_RRM_del"  

#define the threshold used to determine a significant result during wilcoxon test
p <- 0.05


#######################
#DEFINE ASTHETICS OF PLOT

#may need to adjust scale for y-axis on line 

#select which time points to present
Timepoints <- c("20", "40", "60")

#define colours for each cell line - leave blank if not presenting all 4 cell lines
Colour_1_exp <- "#FA9FB5" #sets colour for expansion phase new reference_line
Colour_1_ss <- "#DD3497"  #sets colour for steady state phase new reference_line

Colour_2_exp <- "#A6BDDB" #sets colour for expansion phase new cell line
Colour_2_ss <- "#3690C0"  #sets colour for steady state phase new cell line

Colour_3_exp <- "#FEB24C" #sets colour for expansion phase new cell_line_3
Colour_3_ss <- "#FC4E2A"  #sets colour for steady state phase new cell_line_3

Colour_4_exp <- "#ADDD8E" #sets colour for expansion phase new cell_line_3
Colour_4_ss <- "#41AB5D"  #sets colour for steady state phase new cell_line_3



#set name of how cell lines should be presented in the plot - ensures name consistency with other papers
#if name contains delta/triangle symbol use - "SPEN^"~Delta*"RRM"
Name_1 <- bquote("METTL3_FKBP12"^"F36V") #sets name for new cell line
Name_2 <- bquote("WT")                   #sets name for reference_line
Name_3 <- bquote("Ciz1_KO")              #sets name for cell_line_3
Name_4 <- bquote("SPEN"^~Delta*"RRM")    #sets name for cell_line_4

#set name of plot titles
Title <- "Number of New Xist RNPs Over Time"
x_axis <- "Time [min]"
y_axis <- "Newly Synthesized Xist Count [no. centroids]"


#################################################################################
#STEP 1: load Transcription_Dynamics data for plotting

#define file path to where "Pulse_Chase_Analysis" is located - this is where the compiled data is stored
File_Path <- choose.dir(default = "", caption = "Select Pulse_Chase_Analysis folder, where compiled data is stored")

#open Transcription_Dynamics data for all existing cell lines
if (file.exists(paste(File_Path, Cell_Type, "Transcription_Dynamics", "All_Cell_Lines_Merged", "New_All_Cell_Lines_Transcription_Dynamics_Compile.csv", sep="/"))) { 
  All_Transcription_Dynamics_Data <- read_csv(paste(File_Path, Cell_Type, "Transcription_Dynamics", "All_Cell_Lines_Merged", "New_All_Cell_Lines_Transcription_Dynamics_Compile.csv", sep="/"))
} else{
  stop("New_All_Cell_Lines_Transcription_Dynamics_Compile.csv file does not exist in directory
       - check File_Path input")
}


#checks if new cell line data is stored in the dataframe containing all Transcription_Dynamics data
if (any(All_Transcription_Dynamics_Data$Cell_Line == New_Line_Name)) {
  print("New cell line data is found within Main Dataframe")
} else {
  stop("This new cell line does not exist in Main Dataframe (All_Transcription_Dynamics_Data)
       - make sure to run C1_C2-Centroids_Compilation_Manipulation.R ")
}

#################################################################################
#STEP 2: SELECT WHICH DATA TO PRESENT IN THE PLOT
#ultimately will present 2 sets of box plots

#chose which data to present in box plot
Present_Data <- All_Transcription_Dynamics_Data %>%
  filter(Cell_Line %in%
           c(Refernece_Line, New_Line_Name, Cell_Line_3, Cell_Line_4) &  #select which cell lines to present
           Phase %in%
           c("Initiation", "Maintenance") &                              #select which Phases to present
           Time %in%                                                     #select which Phases to present
           Timepoints) %>%                            
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",          #ensures naming consistency with other papers
                           Phase == "Maintenance" ~ "Steady_State")) %>%                 
  mutate(Key = paste(Phase, Cell_Line, sep="_"))                         #creates key variable to base colour coding off


#order data to be presented in box plot 
Present_Data$Cell_Line <- factor(Present_Data$Cell_Line, 
                                 levels = c(Refernece_Line, New_Line_Name, Cell_Line_3, Cell_Line_4), ordered = TRUE) 

#put in order for the graphs to appear on the x-axis
Present_Data$Key <- factor(Present_Data$Key, 
                           levels = c(paste("Expansion", Refernece_Line, sep="_"), paste("Steady_State", Refernece_Line, sep="_"),
                                      paste("Expansion", New_Line_Name, sep="_"), paste("Steady_State", New_Line_Name, sep="_"),
                                      paste("Expansion", Cell_Line_3, sep="_"), paste("Steady_State", Cell_Line_3, sep="_"),
                                      paste("Expansion", Cell_Line_4, sep="_"), paste("Steady_State", Cell_Line_4, sep="_"), ordered = TRUE))

Present_Data$Time  <- factor(Present_Data$Time,
                                 levels = Timepoints, ordered = TRUE)


#################################################################################
#STEP 4: PLOT THE DATA

#create theme for plots - defines font, text size etc
theme <- theme(plot.title = element_text(family = "Calibri", face = "bold", size = (20)),
               legend.title = element_text(family = "Calibri", size = (16)), 
               legend.text = element_text(family = "Calibri", size = (14)), 
               axis.title = element_text(family = "Calibri", face = "bold", size = (16)),
               axis.text = element_text(family = "Calibri", face = "bold", size = (12)),
               strip.text.x = element_text(family = "Calibri", size = (16)))

#sets name of labels for each cell line
Labels = c(Cell_Line = Name_1,
           Cell_Line = Name_2,
           Cell_Line = Name_3,
           Cell_Line = Name_4)

#calculate the y scale for the plot
Min_y = min(Present_Data$No_Centroids, na.rm = TRUE)
Max_y = max(Present_Data$No_Centroids, na.rm = TRUE)  #adjust value according to graph output as doesn't contain outliers

#make box plot 

box_plot <- ggplot(Present_Data, aes(x= Time, y = No_Centroids, fill = Key)) + #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                        #defines type of plot
  facet_wrap(~Phase, nrow = 1) +                                               #shows box plots for selected cell lines
  scale_fill_manual(values = c(Colour_1_exp, Colour_1_ss,                      #defines the colour of each boxplot 
                               Colour_2_exp, Colour_2_ss,                                               
                               Colour_3_exp, Colour_3_ss,                                                  
                               Colour_4_exp, Colour_4_ss)) +                                                     
  coord_cartesian(ylim = c((Min_y),(90) )) +                                   #adjusts scale so whiskers don't touch the end graph 
  labs(title = Title,                                                          #sets name of the axes
       x = x_axis,
       y = y_axis) +                                          
  theme  

#view the box plot
box_plot

################################################################################
#STEP 4: need to manually export the plots

################################################################################
#STEP 5: preform Wilcox Test on Data to asses if there is a statistical difference

#is there a difference between the 2 data sets?
#is there a general difference between the data sets?
Stat_Diff <- wilcox.test(All_Transcription_Dynamics_Data[which(All_Transcription_Dynamics_Data$Cell_Line==Refernece_Line),]$No_Centroids,
                         All_Transcription_Dynamics_Data[which(All_Transcription_Dynamics_Data$Cell_Line==New_Line_Name),]$No_Centroids)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- wilcox.test(All_Transcription_Dynamics_Data[which(All_Transcription_Dynamics_Data$Cell_Line==Refernece_Line),]$No_Centroids,
                         All_Transcription_Dynamics_Data[which(All_Transcription_Dynamics_Data$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less transcription than the WT 

#if there is a difference, is it statistically more
Diff_More <- wilcox.test(All_Transcription_Dynamics_Data[which(All_Transcription_Dynamics_Data$Cell_Line==Refernece_Line),]$No_Centroids,
                         All_Transcription_Dynamics_Data[which(All_Transcription_Dynamics_Data$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more transcription than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Compare_Both_Phases", Refernece_Line , "vs", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < p) {
  if (Diff_Less$p.value < Diff_More$p.value ){
    outcome_1 <- paste("Overall the new cell line has statistically higher trancription rate than the", Refernece_Line, sep=" ")
  } else {
    outcome_1 <- paste("Overall the new cell line has statistically lower transcription rate than the", Refernece_Line, sep=" ")
  }} else {
    outcome_1 <- paste("Overall there is no statistical difference in transcription rate between the new cell line and the", Refernece_Line, sep=" ")
  }


################################################################################
#STEP 6: preform Wilcox Test on Data to asses if there is a statistical difference between expansion phases in WT and new cell line?

#select data only for the expansion phase
Expansion_Only <- filter(All_Transcription_Dynamics_Data, Phase == "Initiation")

#is there a difference between the 2 data sets?
#is there a general difference between the data sets?
Exp_Stat_Diff <- wilcox.test(Expansion_Only[which(Expansion_Only$Cell_Line==Refernece_Line),]$No_Centroids,
                         Expansion_Only[which(Expansion_Only$Cell_Line==New_Line_Name),]$No_Centroids)
Exp_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Exp_Diff_Less <- wilcox.test(Expansion_Only[which(Expansion_Only$Cell_Line==Refernece_Line),]$No_Centroids,
                         Expansion_Only[which(Expansion_Only$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
Exp_Diff_Less$p.value #p = 1 suggests new cell line less transcription than the WT 

#if there is a difference, is it statistically more
Exp_Diff_More <- wilcox.test(Expansion_Only[which(Expansion_Only$Cell_Line==Refernece_Line),]$No_Centroids,
                         Expansion_Only[which(Expansion_Only$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
Exp_Diff_More$p.value #p = 1 suggests new cell line has more transcription than the WT

#record p values in table
e <- Exp_Stat_Diff$p.value
f <- Exp_Diff_Less$p.value
g <- Exp_Diff_More$p.value
Exp_p_values <- tibble(e, f, g) %>%
  transmute(Comparision = paste("Compare_Both_Phases", Refernece_Line , "vs", New_Line_Name, sep='_'),
            General_Stat_Diff = e,
            Diff_Statistically_Less = f,
            Diff_Statistically_More = g) 

if (Exp_Stat_Diff$p.value < p) {
  if (Exp_Diff_Less$p.value < Exp_Diff_More$p.value ){
    outcome_2 <- paste("The new cell line has statistically higher trancription rate during expansion phase than the", Refernece_Line, sep=" ")
  } else {
    outcome_2 <- paste("The new cell line has statistically lower transcription rate during expansion phase than the", Refernece_Line, sep=" ")
  }} else {
    outcome_2 <- paste("There is no statistical difference in transcription rate during expansion phase between the new cell line and the", Refernece_Line, sep=" ")
  }

################################################################################
#STEP 7: preform Wilcox Test on Data to asses if there is a statistical difference between steady-state phases in WT and new cell line?

#select data only for the expansion phase
Steady_State_Only <- filter(All_Transcription_Dynamics_Data, Phase == "Maintenance")

#is there a difference between the 2 data sets?
#is there a general difference between the data sets?
SS_Stat_Diff <- wilcox.test(Steady_State_Only[which(Steady_State_Only$Cell_Line==Refernece_Line),]$No_Centroids,
                             Steady_State_Only[which(Steady_State_Only$Cell_Line==New_Line_Name),]$No_Centroids)
SS_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
SS_Diff_Less <- wilcox.test(Steady_State_Only[which(Steady_State_Only$Cell_Line==Refernece_Line),]$No_Centroids,
                             Steady_State_Only[which(Steady_State_Only$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "less")
SS_Diff_Less$p.value #p = 1 suggests new cell line less transcription than the WT 

#if there is a difference, is it statistically more
SS_Diff_More <- wilcox.test(Steady_State_Only[which(Steady_State_Only$Cell_Line==Refernece_Line),]$No_Centroids,
                             Steady_State_Only[which(Steady_State_Only$Cell_Line==New_Line_Name),]$No_Centroids, alternative = "greater")
SS_Diff_More$p.value #p = 1 suggests new cell line has more transcription than the WT

#record p values in table
h <- SS_Stat_Diff$p.value
i <- SS_Diff_Less$p.value
j <- SS_Diff_More$p.value
SS_p_values <- tibble(h, i, j) %>%
  transmute(Comparision = paste("Compare_Both_Phases", Refernece_Line , "vs", New_Line_Name, sep='_'),
            General_Stat_Diff = h,
            Diff_Statistically_Less = i,
            Diff_Statistically_More = j) 

if (SS_Stat_Diff$p.value < p) {
  if (SS_Diff_Less$p.value < SS_Diff_More$p.value ){
    outcome_3 <- paste("The new cell line has statistically higher trancription rate during steady-state phase than the", Refernece_Line, sep=" ")
  } else {
    outcome_3 <- paste("The new cell line has statistically lower transcription rate during steady-state phase than the", Refernece_Line, sep=" ")
  }} else {
    outcome_3 <- paste("There is no statistical difference in transcription rate during steady-state phase between the new cell line and the", Refernece_Line, sep=" ")
  }


#################################
#STEP 8: preform Wilcox Test on Data to asses if there is a statistical difference between expansion and SS phases in WT and new cell line?

#select data only for the expansion phase
New_Line_Name_Only <- filter(All_Transcription_Dynamics_Data, Cell_Line == New_Line_Name)

#is there a difference between the 2 data sets?
#is there a general difference between the data sets?
New_Line_Name_Stat_Diff <- wilcox.test(New_Line_Name_Only[which(New_Line_Name_Only$Phase=="Initiation"),]$No_Centroids,
                            New_Line_Name_Only[which(New_Line_Name_Only$Phase=="Maintenance"),]$No_Centroids)
New_Line_Name_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
New_Line_Name_Diff_Less <- wilcox.test(New_Line_Name_Only[which(New_Line_Name_Only$Phase=="Initiation"),]$No_Centroids,
                            New_Line_Name_Only[which(New_Line_Name_Only$Phase=="Maintenance"),]$No_Centroids, alternative = "less")
New_Line_Name_Diff_Less$p.value #p = 1 suggests new cell line less transcription than the New_Line

#if there is a difference, is it statistically more
New_Line_Name_Diff_More <- wilcox.test(New_Line_Name_Only[which(New_Line_Name_Only$Phase=="Initiation"),]$No_Centroids,
                            New_Line_Name_Only[which(New_Line_Name_Only$Phase=="Maintenance"),]$No_Centroids, alternative = "greater")
New_Line_Name_Diff_More$p.value #p = 1 suggests new cell line has more transcription than the New Line

#record p values in table
n <- New_Line_Name_Stat_Diff$p.value
o <- New_Line_Name_Diff_Less$p.value
q <- New_Line_Name_Diff_More$p.value
New_Line_Name_p_values <- tibble(n, o, q) %>%
  transmute(Comparision = paste(New_Line_Name, "Internal_Exp_vs_SS", sep="_"),
            General_Stat_Diff = n,
            Diff_Statistically_Less = o,
            Diff_Statistically_More = q) 

#compare all p values together
p_values <- p_values %>%
  bind_rows(Exp_p_values,
            SS_p_values,
            New_Line_Name_p_values) 

if (New_Line_Name_Stat_Diff$p.value < p) {
  if (New_Line_Name_Diff_Less$p.value < New_Line_Name_Diff_More$p.value ){
    outcome_4 <- "The new cell line has statistically higher trancription rate during steady-state phase than the expansion phase"
  } else {
    outcome_4 <- "The new cell line has statistically lower transcription rate during steady-state phase than the expansion phase"
  }} else {
    outcome_4 <- "There is no statistical difference in transcription rate during steady-state phase and the expansion phase in the new cell line"
  }

#save p values
Save_Path <- paste(File_Path, Cell_Type, "Transcription_Dynamics", New_Line_Name, sep="/")

write_csv(p_values, paste(Save_Path, paste(New_Line_Name, "_Transcription_Dynamics_p_values.csv", sep=""), sep="/"))


#provide outcomes of Wilcoxon tests

if (Exp_Stat_Diff$p.value < p) {
  print(outcome_1)
  print(outcome_2)
  print(outcome_3)
  print(outcome_4)
} else {
  print(outcome_1)
  }