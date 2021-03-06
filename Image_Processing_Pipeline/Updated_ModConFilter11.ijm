//## Copyright (C) 2016 Ezequiel Miron <eze.miron@bioch.ox.ac.uk>
//##
//## This program is free software: you can redistribute it and/or modify
//## it under the terms of the GNU General Public License as published by
//## the Free Software Foundation, either version 3 of the License, or
//## (at your option) any later version.
//##
//## This program is distributed in the hope that it will be useful,
//## but WITHOUT ANY WARRANTY; without even the implied warranty of
//## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//## GNU General Public License for more details.
//##
//## You should have received a copy of the GNU General Public License
//## along with this program.  If not, see <http://www.gnu.org/licenses/>.

dir1 = getDirectory("Choose Source Directory ");
dir2 = getDirectory("Choose Destination Directory ");
list = getFileList(dir1);     //gets list of files in dir1
setBatchMode(false);
for (i=0; i<list.length; i++) 
{
    file1 = dir1 + list[i]; // defining THR filepath
    
    if (endsWith(file1, "THR.tif")) 
    {
    	th = list[i];
	    MC = replace(th, "SIR_THR", "MCN");
	    file2 = dir1 + MC;  // defining MCN filepath

	    // opems THR file
		run("Bio-Formats Importer", "open="+file1+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); 
		//open MCN file
		run("Bio-Formats Importer", "open="+file2+" autoscale color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT"); 

		//defines parameters
		sigma = 0.8;
		min_thresh = 4.0;
		
		thresh = filename(th);
		selectWindow(th);
		rename(thresh);
		MCNR = filename(MC);
		selectWindow(MC);
		rename(MCNR);
		
		/*print(min_thresh);*/ 
		MCF(thresh, MCNR, min_thresh);
		run("Grays");
		contrast_stack(thresh+"_MCF");

		run("Split Channels"); //seperates channels
		selectWindow("C1-"+thresh+"_MCF"); //saves seperate channels
		saveAs("TIFF", dir2 + "C1-"+thresh+"_MCF");
		selectWindow("C2-"+thresh+"_MCF");
		run("Green"); //ensures C2 image signal is green 
		run("Duplicate...", "duplicate"); //duplicate C2 for use in Step 2c macro
		saveAs("TIFF", dir2 + "C2-"+thresh+"_MCF");
		saveAs("TIFF", dir2 + "C2-"+thresh+"_MCF-1");
		selectWindow("C3-"+thresh+"_MCF");
		run("Blue"); //ensures C3 image signal is blue
		run("Duplicate...", "duplicate"); //duplicate C2 for use in Step 2c macro
		saveAs("TIFF", dir2 + "C3-"+thresh+"_MCF");
		saveAs("TIFF", dir2 + "C3-"+thresh+"_MCF-1");
		close("*");
		}
}

print("Completed Running ModConFilter11 For Data Set with correct 255 vs 65355");
print("Completed Spltting Channels for Images");
print("Completed Duplicating C2 and C3 Channels");
print("Saved all Relevent Images");
print("Next mannuslly Threshold C2/3 MCF-1 images")

function MCF (THR, MCN, min_thresh){
	selectWindow(MCN);
	Stack.getDimensions(width, height, channels, slices, frames);
	run("Scale...", "x=2 y=2 z=1.0 width=" + width + " height=" + height + "depth=" + slices + "interpolation=Bicubic average process create title="+MCN+"Scale");
    selectWindow(MCN+"Scale");
    getMinAndMax(min,max);
    setBatchMode(true);
    for (i=1; i<=nSlices; i++) {
       setSlice(i);
       getMinAndMax(min, max);
       setThreshold(min_thresh, max);
       run("Create Mask");
       run("Copy");
       selectWindow(MCN+"Scale");
       run("Paste");
       setSlice(1);
       resetThreshold;}
    setBatchMode(false);
    run("16-bit");
    run("Divide...", "value=255 stack"); //65535
	getMinAndMax(min, max);
	if (max>40000){		//ensures britgthness is not ridculously high
		run("Multiply...", "value=255 stack");
		run("Divide...", "value=65535 stack");
	}  
    rename("Mask "+THR);
    imageCalculator("Multiply create stack", "Mask "+THR, THR);
    selectWindow("Result of Mask "+THR);
    rename(THR+"_MCF");
    run("Gaussian Blur...", "sigma="+sigma+" stack");
}
	
function filename (name){
	dot = indexOf(name, ".");
    if (dot>=0) name = substring(name, 0, dot);
    else if (dot == NaN) name = name;
    return name;
}

function substack(thresh, MCNR, array){
  for (i=0; i<channels; i++){
    if (array[i] == 1){
	  selectWindow(thresh);
	  Stack.getDimensions(width, height, channels, slices, frames);
	  run("Make Substack...", "channels="+i+" slices=1-"+slices);
      resetMinAndMax;
      rename("C"+i+"-"+thresh);
      selectWindow(MCNR);
      run("Make Substack...", "channels="+i+" slices=1-"+slices);
      resetMinAndMax;
      rename("C"+i+"-"+MCNR);
}}

function contrast_stack(stack){
  Stack.getDimensions(width, height, channels, slices, frames);
  setSlice(round(slices/2));
  run("Enhance Contrast", "saturated=0.35");
}
  