//manually Therhold all C2/3 MCF-1 files before running script

dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory ");
list = getFileList(dir1);     //gets list of files in dir1
setBatchMode(false);
for (i=0; i<list.length; i++) 
{
    C1_file1 = dir1 + list[i]; // defining MCF filepath
    
    if (startsWith(C1_file1, dir1+"C1") && endsWith(C1_file1, "MCF.tif")) 
    {
    	C1_Original = list[i]; //defines C1 files
	    C1_Threshold = replace(C1_Original, "MCF", "MCF-1");
	    C1_file2 = dir1 + C1_Threshold;  
	    
	    C2_file1 = replace(C1_file1, "C1", "C2"); //defines C2 Files
	    C2_Original = replace(C1_Original, "C1", "C2"); 
	    C2_Threshold = replace(C1_Threshold, "C1", "C2");
	    C2_file2 = replace(C1_file2, "C1", "C2"); 
	       
	    C3_file1 = replace(C2_file1, "C2", "C3"); //defines C3 Files
	    C3_Original = replace(C2_Original, "C2", "C3"); 
	    C3_Threshold = replace(C2_Threshold, "C2", "C3");
	    C3_file2 = replace(C2_file2, "C2", "C3");

	    // opens C1 MCN files
		run("Bio-Formats Importer", "open="+C1_file1+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); 
		//open C1 MCN-1 files
		//run("Bio-Formats Importer", "open="+C1_file2+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); 

		// opens C2 MCN files
		run("Bio-Formats Importer", "open="+C2_file1+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); 
		//open C2 MCN-1 files
		run("Bio-Formats Importer", "open="+C2_file2+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); 

		// opens C3 MCN files
		run("Bio-Formats Importer", "open="+C3_file1+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); 
		//open C3 MCN-1 files
		run("Bio-Formats Importer", "open="+C3_file2+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); 

 
    	//creates Result of C2 image which has background removed
		selectWindow(C2_Threshold); 
		run("Divide...", "value=255.000 stack");
		imageCalculator("Multiply create stack", C2_Threshold, C2_Original);
		saveAs("TIFF", dir2 + "Results of "+C2_Threshold);

		//creates Result of C3 image which has background removed
		selectWindow(C3_Threshold); 
		run("Divide...", "value=255.000 stack");
		imageCalculator("Multiply create stack", C3_Threshold, C3_Original);
		saveAs("TIFF", dir2 + "Results of "+C3_Threshold);

		//merge channels
		run("Merge Channels...", "c1="+C1_Original+" c2=[Results of "+C2_Threshold+"] c3=[Results of "+C3_Threshold+"] create");		
		//run("Merge Channels...", "c1=C1-20201014_RNA_SPLIT_Mettl3dTAG_Initiation_0_01_FUS_SIR_THR_MCF.tif c2=[Result of C2-20201014_RNA_SPLIT_Mettl3dTAG_Initiation_0_01_FUS_SIR_THR_MCF-1.tif] c3=[Result of C3-20201014_RNA_SPLIT_Mettl3dTAG_Initiation_0_01_FUS_SIR_THR_MCF-1.tif] create");
		
		selectWindow(substring(C1_Original, 3)); 
		run("Grays"); //ensures channel 1 is gray not red
		saveAs("TIFF", dir2 + getTitle());

		//close all windows
		close("*");
			
    }

    /*else {
    	print("FAILED :(");
    }*/
} 

print("Compiled Images with C2/3 Background Removed");
print("Next run Chromagnon on resulting SIR_THR_MCF files"); 