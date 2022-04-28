#################################################################################
#SRCIPT CREATES BOX PLOT FOR CLOUD VOLUME DATA AND PREFORMS WILCOX STAT TEST
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
# New_Line_Name <- "Mettl3_dTAG"       
# 
# Reference_Line_Name <- "WT"        #Reference_Line_Name_Name will be used to perform stats test with
# 
# Cell_Line_3 <- "Ciz1_KO"
# 
# Cell_Line_4 <- "SPEN_RRM_del"  
# 
# #define the threshold used to determine a significant result during wilcoxon test
# p <- 0.05
# 
# 
# #######################
# #DEFINE ASTHETICS OF PLOT
# 
# #may need to adjust scale for y-axis on line 
# 
# #define colours for each cell line - leave blank if not presenting all 4 cell lines
# Colour_1_exp <- "#FA9FB5" #sets colour for expansion phase reference line
# Colour_1_ss <- "#DD3497"  #sets colour for steady state phase reference line
# 
# Colour_2_exp <- "#A6BDDB" #sets colour for expansion phase new cell line
# Colour_2_ss <- "#3690C0"  #sets colour for steady state phase new cell line
# 
# Colour_3_exp <- "#FEB24C" #sets colour for expansion phase new cell_line_3
# Colour_3_ss <- "#FC4E2A"  #sets colour for steady state phase new cell_line_3
# 
# Colour_4_exp <- "#ADDD8E" #sets colour for expansion phase new cell_line_3
# Colour_4_ss <- "#41AB5D"  #sets colour for steady state phase new cell_line_3
# 
# 
# #set name of how cell lines should be presented in the plot - ensures name consistency with other papers
# #if name contains delta/triangle symbol use - "SPEN^"~Delta*"RRM"
# Name_1 <- bquote("METTL3_FKBP12"^"F36V") #sets name for new cell line
# Name_2 <- bquote("WT")                   #sets name for Reference_Line_Name_Name
# Name_3 <- bquote("Ciz1_KO")              #sets name for cell_line_3
# Name_4 <- bquote("SPEN"^~Delta*"RRM")    #sets name for cell_line_4

#set name of plot titles
Title <- "Xist Cluster Volume"
x_axis <- "Phase"
y_axis <- bquote("Xist"~"Territory"~"Volume"~"["*mu*"m"^"3"*"]")


#################################################################################
#STEP 1: load Cloud_Volume data for plotting

#define file path to where "Pulse_Chase_Analysis" is located - this is where the compiled data is stored
File_Path <- Output_File_Path

#set file path to the location of the all cloud volume files 
Input_Path <- paste(File_Path, Cell_Type, "Cloud_Volume", "All_Cell_Lines", sep="/")

#stores name of all Cloud_Volume_Compile.csv files , in the directory, into a vector
Files <- list.files(path = Input_Path, pattern = "Cloud_Volume_Compile.csv", full.names  = TRUE)   

#creates a list containing the Cloud_Volume_Compile.csv files
File_List <- lapply(Files, read_csv, col_types="cccncn")

#concatenates all the individual Cloud_Volume_Compile.csv files into 1 tibble
All_Cloud_Volume_Data <- bind_rows(File_List)


#checks if new cell line data is stored in the dataframe containing all Cloud_Volume data
if (!any(All_Cloud_Volume_Data$Cell_Line == New_Line_Name)) {
  stop("This new cell line does not exist in Main Dataframe (All_Cloud_Volume_Data)
       - make sure to run Cloud_Volume_Compilation_Manipulation.R ")
}

#################################################################################
#STEP 2: SELECT WHICH DATA TO PRESENT IN THE PLOT
#ultimately will present 2 sets of box plots

#chose which data to present in box plot
Present_Data <- All_Cloud_Volume_Data %>%
  filter(Cell_Line %in%
           c(Reference_Line_Name, New_Line_Name, Cell_Line_3, Cell_Line_4) &  #select which cell lines to present
        Phase %in%
          c("Initiation", "Maintenance")) %>%                            #select which Phases to present
  mutate(Phase = case_when(Phase == "Initiation" ~ "Expansion",          #ensures naming consistency with other papers
                           Phase == "Maintenance" ~ "Steady_State")) %>%                 
  mutate(Key = paste(Phase, Cell_Line, sep="_"))                         #creates key variable to base colour coding off


#order data to be presented in box plot 
Present_Data$Cell_Line <- factor(Present_Data$Cell_Line, 
                                 levels = c(Reference_Line_Name, New_Line_Name, Cell_Line_3, Cell_Line_4), ordered = TRUE,
                                 labels = c(Name_2, Name_1, Name_3, Name_4))  

#put in order for the graphs to appear on the x-axis
Present_Data$Key <- factor(Present_Data$Key, 
                                 levels = c(paste("Expansion", Reference_Line_Name, sep="_"), paste("Steady_State", Reference_Line_Name, sep="_"),
                                            paste("Expansion", New_Line_Name, sep="_"), paste("Steady_State", New_Line_Name, sep="_"),
                                            paste("Expansion", Cell_Line_3, sep="_"), paste("Steady_State", Cell_Line_3, sep="_"),
                                            paste("Expansion", Cell_Line_4, sep="_"), paste("Steady_State", Cell_Line_4, sep="_"), ordered = TRUE))


#################################################################################
#STEP 3: PLOT THE DATA

#sets name of labels for each cell line
Labels = c(Cell_Line = Name_1,
           Cell_Line = Name_2,
           Cell_Line = Name_3,
           Cell_Line = Name_4)

#calculate the y scale for the plot
Min_y = min(Present_Data$Cloud_Volume, na.rm = TRUE) 
Max_y = max(Present_Data$Cloud_Volume, na.rm = TRUE)  

#make box plot
box_plot <- ggplot(Present_Data, aes(x= Phase, y = Cloud_Volume, fill = Key)) + #defines data to present
  geom_boxplot(outlier.shape = NA) +                                            #defines type of plot
  facet_wrap(~Cell_Line, nrow = 1, labeller = label_parsed) +                   #separates cell lines, and sets labels
  scale_fill_manual(values = c(Colour_1_exp, Colour_1_ss,                       #defines the colour of each boxplot 
                               Colour_2_exp, Colour_2_ss,                                               
                               Colour_3_exp, Colour_3_ss,                                                  
                               Colour_4_exp, Colour_4_ss)) +                                                     
  coord_cartesian(ylim = c((Min_y-10),(1200) )) +                               #adjusts scale so whiskers don't touch the end graph 
  labs(title = Title,                                                           #sets name of the axes
       x = x_axis,
       y = y_axis) 

#view the box plot
print(box_plot)


################################################################################
#STEP 4: need to manually export the plots

################################################################################
#STEP 5: preform Wilcox Test on Data to asses if there is a statistical difference

#is there a general difference between the data sets?
Stat_Diff <- wilcox.test(All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==Reference_Line_Name),]$Cloud_Volume,
                         All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==New_Line_Name),]$Cloud_Volume)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- wilcox.test(All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==Reference_Line_Name),]$Cloud_Volume,
                         All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- wilcox.test(All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==Reference_Line_Name),]$Cloud_Volume,
                         All_Cloud_Volume_Data[which(All_Cloud_Volume_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "greater")
Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT

#record p values in table
a <- Stat_Diff$p.value
b <- Diff_Less$p.value
c <- Diff_More$p.value
p_values <- tibble(a, b, c) %>%
  transmute(Comparision = paste("Compare_Both_Phases", Reference_Line_Name , "vs", New_Line_Name, sep='_'),
            General_Stat_Diff = a,
            Diff_Statistically_Less = b,
            Diff_Statistically_More = c) 

if (Stat_Diff$p.value < p) {
  if (Diff_Less$p.value < Diff_More$p.value){
    outcome_1 <- paste("Overall, the new cell line (", New_Line_Name,") has statistically bigger clouds than the", Reference_Line_Name, sep=" ")
  } else {
    outcome_1 <- paste("Overall, the new cell line (", New_Line_Name,") has statistically smaller clouds than the", Reference_Line_Name, sep=" ")
  }} else {
    outcome_1 <- paste("Overall, there is no statistical difference in cloud size between the new cell line (", New_Line_Name,") and the", Reference_Line_Name, sep=" ")
  }

print(outcome_1)

#is there a general difference during the expansion phase between the data sets?
Expansion_Data <- filter(All_Cloud_Volume_Data, Phase == "Initiation")

#only perform stats test if data exists
#create empty tibble for in case it doesn't exist
Exp_p_values <- tibble()
if (New_Line_Name %in% Expansion_Data$Cell_Line) {

  Exp_Stat_Diff <- wilcox.test(Expansion_Data[which(Expansion_Data$Cell_Line==Reference_Line_Name),]$Cloud_Volume,
                               Expansion_Data[which(Expansion_Data$Cell_Line==New_Line_Name),]$Cloud_Volume)
  Exp_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference
  
  
  #if there is a difference, is it statistically less 
  Exp_Diff_Less <- wilcox.test(Expansion_Data[which(Expansion_Data$Cell_Line==Reference_Line_Name),]$Cloud_Volume,
                               Expansion_Data[which(Expansion_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "less")
  Exp_Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 
  
  #if there is a difference, is it statistically more>
  Exp_Diff_More <- wilcox.test(Expansion_Data[which(Expansion_Data$Cell_Line==Reference_Line_Name),]$Cloud_Volume,
                               Expansion_Data[which(Expansion_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "greater")
  Exp_Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT
  
  #record p values in table
  e <- Exp_Stat_Diff$p.value
  f <- Exp_Diff_Less$p.value
  g <- Exp_Diff_More$p.value
  Exp_p_values <- tibble(e, f, g) %>%
    transmute(Comparision = paste("Compare_Exp_Phase", Reference_Line_Name , "vs", New_Line_Name, sep='_'),
              General_Stat_Diff = e,
              Diff_Statistically_Less = f,
              Diff_Statistically_More = g) 
  
  if (Exp_Stat_Diff$p.value < p) {
    if (Exp_Diff_Less$p.value < Exp_Diff_More$p.value){
      outcome_2 <- paste("The new cell line (", New_Line_Name,") has statistically bigger clouds, during the expansion phase, than the", Reference_Line_Name, sep=" ")
    } else {
      outcome_2 <- paste("The new cell line (", New_Line_Name,") has statistically smaller clouds, during the expansion phase, than the", Reference_Line_Name, sep=" ")
    }} else {
      outcome_2 <- paste("There is no statistical difference in cloud size, during the expansion phase, between the new cell line (", New_Line_Name,") and the", Reference_Line_Name, sep=" ")
    }
  
  print(outcome_2)
}

#is there a general difference during the Steady_State phase between the data sets?
Steady_State_Data <- filter(All_Cloud_Volume_Data, Phase == "Maintenance")

#only perform stats test if data exists
#create empty tibble for in case it doesn't exist
SS_p_values <- tibble()
if (New_Line_Name %in% Steady_State_Data$Cell_Line) {

  SS_Stat_Diff <- wilcox.test(Steady_State_Data[which(Steady_State_Data$Cell_Line==Reference_Line_Name),]$Cloud_Volume,
                              Steady_State_Data[which(Steady_State_Data$Cell_Line==New_Line_Name),]$Cloud_Volume)
  SS_Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference
  
  
  #if there is a difference, is it statistically less 
  SS_Diff_Less <- wilcox.test(Steady_State_Data[which(Steady_State_Data$Cell_Line==Reference_Line_Name),]$Cloud_Volume,
                              Steady_State_Data[which(Steady_State_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "less")
  SS_Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 
  
  #if there is a difference, is it statistically more>
  SS_Diff_More <- wilcox.test(Steady_State_Data[which(Steady_State_Data$Cell_Line==Reference_Line_Name),]$Cloud_Volume,
                              Steady_State_Data[which(Steady_State_Data$Cell_Line==New_Line_Name),]$Cloud_Volume, alternative = "greater")
  SS_Diff_More$p.value #p = 1 suggests new cell line has more coupling than the WT
  
  #record p values in table
  h <- SS_Stat_Diff$p.value
  i <- SS_Diff_Less$p.value
  j <- SS_Diff_More$p.value
  SS_p_values <- tibble(h, i, j) %>%
    transmute(Comparision = paste("Compare_SS_Phase", Reference_Line_Name , "vs", New_Line_Name, sep='_'),
              General_Stat_Diff = h,
              Diff_Statistically_Less = i,
              Diff_Statistically_More = j) 
  
  if (SS_Stat_Diff$p.value < p) {
    if (SS_Diff_Less$p.value < SS_Diff_More$p.value){
      outcome_3 <- paste("The new cell line (", New_Line_Name,") has statistically bigger clouds, during the Steady_State phase, than the", Reference_Line_Name, sep=" ")
    } else {
      outcome_3 <- paste("The new cell line (", New_Line_Name,") has statistically smaller clouds, during the Steady_State phase, than the", Reference_Line_Name, sep=" ")
    }} else {
      outcome_3 <- paste("There is no statistical difference in cloud size, during the Steady_State phase, between the new cell line (", New_Line_Name,") and the", Reference_Line_Name, sep=" ")
    }
  
  print(outcome_3)
}

#combine all p values
p_values <- p_values %>%
  bind_rows(Exp_p_values,
            SS_p_values)

#save p values
Save_Path <- paste(File_Path, Cell_Type, "Cloud_Volume", New_Line_Name, sep="/")

write_csv(p_values, paste(Save_Path, paste(New_Line_Name, "_Cloud_Volume_Wilcox_p_values.csv", sep=""), sep="/"))



