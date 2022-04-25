#############################################################################################################################################
#SCRIPT CREATES BOX PLOT FOR DENSITY DATA AND PREFORMS WILCOX STAT TEST
#please email Holly Roach at hmroach@hotmail.co.uk if you have any questions
#############################################################################################################################################

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
# Colour_1 <- "#F768A1" #sets colour for reference line
# 
# Colour_2 <- "#67A9CF" #sets colour for new cell line
# 
# Colour_3 <- "#FD8D3C"  #sets colour for cell_line_3
# 
# Colour_4 <- "#78C679"  #sets colour for cell_line_4
# 
# 
# #set name of how cell lines should be presented in the plot - ensures name consistency with other papers
# #if name contains delta/triangle symbol use - "SPEN^"~Delta*"RRM"
# Name_1 <- bquote("METTL3_FKBP12"^"F36V") #sets name for new cell line
# Name_2 <- bquote("WT")                   #sets name for Reference_Line_Name
# Name_3 <- bquote("Ciz1_KO")              #sets name for cell_line_3
# Name_4 <- bquote("SPEN"^~Delta*"RRM")    #sets name for cell_line_4

#set name of plot titles
Title <- "Xist RNP Density"
x_axis <- "Cell_Line"
y_axis <- "Distance between Xist RNPs [nm]"

  
#################################################################################
#STEP 1: load density data for plotting

#define file path to where "Pulse_Chase_Analysis" is located - this is where the compiled data is stored
File_Path <- Output_File_Path

#set file path to the location of the all Density files 
Input_Path <- paste(File_Path, Cell_Type, "Density", "All_Cell_Lines", sep="/")

#stores name of all Density_Compile.csv files , in the directory, into a vector
Files <- list.files(path = Input_Path, pattern = "Density_Compile.csv", full.names  = TRUE)   

#creates a list containing the Density_Compile.csv files
File_List <- lapply(Files, read_csv, col_types="ccccncn")

#concatenates all the individual Density_Compile.csv files into 1 tibble
All_Density_Data <- bind_rows(File_List)


#checks if new cell line data is stored in the dataframe containing all density data
if (!(any(All_Density_Data$Cell_Line == New_Line_Name))) {
  stop("This new cell line does not exist in Main Dataframe (Data_Combined)
       - make sure to run Cloud_Volume_Data_Clean.R")
}

#################################################################################
#STEP 2: SELECT WHICH DATA TO PRESENT IN THE PLOT
#ultimately will present 2 sets of box plots

#chose which data to present in box plot
Present_Data <- All_Density_Data %>%
  filter(Cell_Line %in%
           c(Reference_Line_Name, New_Line_Name, Cell_Line_3, Cell_Line_4)) #select which cell lines to present


#order data to be presented in box plot 
Present_Data$Cell_Line <- factor(Present_Data$Cell_Line, 
                                 levels = c(Reference_Line_Name, New_Line_Name, Cell_Line_3, Cell_Line_4), ordered = TRUE,
                                 labels = c(Name_2, Name_1, Name_3, Name_4))  


#################################################################################
#STEP 3: PLOT THE DATA

#sets name of labels for each cell line
Labels = c(Name_2,
           Name_1,
           Name_3,
           Name_4)

#calculate the y scale for the plot
Min_y = min(Present_Data$Distance, na.rm = TRUE) 
Max_y = max(Present_Data$Distance, na.rm = TRUE) 

#make box plot
box_plot <- ggplot(Present_Data, aes(x= Cell_Line, y = Distance, fill = Cell_Line)) +    #defines data to present
  geom_boxplot(outlier.shape = NA,  ) +                                                  #defines type of plot
  scale_x_discrete(labels = Labels) +                                                    #adds labels to plot
  scale_fill_manual(values = c(Colour_1,                                                 #defines the colour of each boxplot 
                               Colour_2,                                                 
                               Colour_3,                                                  
                               Colour_4),
                    labels = Labels) +                                                     
  coord_cartesian(ylim = c((Min_y),(2500) )) +                                           #adjusts scale so whiskers don't touch the end graph 
  labs(title = Title,                                                                    #sets name of the axes
       x = x_axis,
       y = y_axis)    

#view the box plot
print(box_plot)


################################################################################
#STEP 4: need to manually export the plots

################################################################################
#STEP 5: preform Wilcox Test on Data to asses if there is a statistical difference

#is there a general difference between the data sets?
Stat_Diff <- wilcox.test(All_Density_Data[which(All_Density_Data$Cell_Line==Reference_Line_Name),]$Distance,
                         All_Density_Data[which(All_Density_Data$Cell_Line==New_Line_Name),]$Distance)
Stat_Diff$p.value #p < 0 suggests a statistical difference / p > 0 suggests not statistical difference


#if there is a difference, is it statistically less 
Diff_Less <- wilcox.test(All_Density_Data[which(All_Density_Data$Cell_Line==Reference_Line_Name),]$Distance,
                         All_Density_Data[which(All_Density_Data$Cell_Line==New_Line_Name),]$Distance, alternative = "less")
Diff_Less$p.value #p = 1 suggests new cell line less coupling than the WT 

#if there is a difference, is it statistically more>
Diff_More <- wilcox.test(All_Density_Data[which(All_Density_Data$Cell_Line==Reference_Line_Name),]$Distance,
                         All_Density_Data[which(All_Density_Data$Cell_Line==New_Line_Name),]$Distance, alternative = "greater")
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
    outcome_1 <- paste("The new cell line (", New_Line_Name,") has statistically greater distances between Xist RNPs (lower density) than the", Reference_Line_Name, sep=" ")
  } else {
    outcome_1 <- paste("The new cell line (", New_Line_Name,") has statistically shorter distances between Xist RNPs (higher density) than the", Reference_Line_Name, sep=" ")
  }} else {
    outcome_1 <- paste("There is no statistical difference in Xist RNP distances (density) between the new cell line (", New_Line_Name,") and the", Reference_Line_Name, sep=" ")
  }

#save p values
Save_Path <- paste(File_Path, Cell_Type, "Density", New_Line_Name, sep="/")

write_csv(p_values, paste(Save_Path, paste(New_Line_Name, "_Density_Wilcox_p_values.csv", sep=""), sep="/"))


#returns statements of significance
print(outcome_1)

