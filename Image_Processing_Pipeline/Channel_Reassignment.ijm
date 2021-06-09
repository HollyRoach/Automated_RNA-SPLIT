//before running script make sure channels have been aligned using chromagon

dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory ");
list = getFileList(dir1);     //gets list of files in dir1
setBatchMode(false);
for (i=0; i<list.length; i++) 
{
    ALN_file1 = dir1 + list[i]; // defining ALN_full filepath

    if (endsWith(ALN_file1, "ALN_full.tif")) 
    {	
    	 // opens ALN files
    	 open(ALN_file1); 

		ALN = getTitle();
    	 
    	 // spilts ALN file into C1, C2, C3
		run("Split Channels");
   
		//define images
    	C1_ALN = "C1-"+ALN;
    	C2_ALN = "C2-"+ALN;
    	C3_ALN = "C3-"+ALN;
		
		//merge channels: C1= old C3, C2= old C2, remove C1
	   run("Merge Channels...", "c1="+C3_ALN+" c2="+C2_ALN+" create");
	   selectWindow("Composite");
	   ALN_Save_Name = replace(ALN_file1, "_ALN_full", "_ALN_Reassignment");
	   saveAs("TIFF", ALN_Save_Name);

		//close all windows
		close("*");
    }

}

print("Completed Channel Reassignment");
print("Next Crop Clouds");
		
	    