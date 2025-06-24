/*	Macro to measure and display fluorescence lifetimes or intensitues of individual cells in time-lapse experiments.
 * 	Also used for detecting hits in a screen.
 * 	
 * 	Input:
 * 	► 2-channel .tif files representing two lifetime-components measured with TCSPC (e.g. files exported from Stellaris/SP8)
 * 	  Optionally a third channel with a nuclear marker			
 * 	► .fli files from the Lambert Instruments Frequency-Domain FLIM microscope
 * 	► Leica TauContrast, TauSeparation timelapse images
 * 	► 2-channel .tif files with ratio imaging
 * 	► single-channel timelapse images with intensities
 *
 *	Brief workflow: 			  
 *	► Segment cells (create labelmap) using Cellpose, and optionally segment nuclei with StarDist
 *	► If a nuclear marker is used, nuclear labels are assigned to cellular labels. Single nuclei without a detected cell are also assigned a cellular label
 *	► Measure the intensity-weighted lifetime of all labels
 *	► Display the lifetime traces in a graph
 *	► Display the lifetime traces in a kymograph-like image
 *	► Determine stimulation and calibration time points from the second derivative of the average lifetime trace of all labels
 *	► Sort the kymograph on cellular response
 *	
 *	 * Requires the following update sites:
 * - CLIJ
 * - CLIJ2
 * - CSBDeep
 * - ImageScience
 * - IJ-PB Plugins
 * - PTBIOP
 * - SCF MPI CBG
 * - StarDist
 * You also need a working Cellpose Python environment, and the 'Turbo' LUT (https://github.com/cleterrier/ChrisLUTs/blob/master/Turbo.lut)
 * For FD-FLIM analysis you need the plugin by Rolf Harkes (https://github.com/rharkes/fdFLIM/)
 * 
 *	Author: Bram van den Broek, The Netherlands Cancer Institute, b.vd.broek@nki.nl
 *	
 *	version 1.4:
 *	- Added flow threshold parameter for Cellpose
 *	- Added sensitivity for detecting stimulation and calibration
 *	- Added possibility to restrict the image for Cellpose to the first n frames
 *	- Small bug fixes
 *	
 *	version 1.5:
 *	- Added feature to overlay lifetime with the intensity movie instead of the still intensity image
 *	- Keep ROI manager open to allow saving the ROIs at the end. Due to this, converting the labelmap to ROIs is now a bit slower.
 *	- More small bug fixes
 *	
 *	version 1.7:
 *	- Reinstated the nuclei detection part
 *	
 *	Version 2.0:
 *	- United with the script that was used for screening, and improved that part on many aspects, notably:
 *	- Cell-nucleus linking faster and more robust, now using CLIJ2_argMaximumZProjection()
 *	- Added hit selection criteria
 *	- Create a top N hits selection
 *	- Plotting the hits and tiles
 *	
 *	Version 2.1:
 *	- NaN pixels in the weighted lifetime do no longer cause missing pixels in the overlay
 *	
 *	Version 2.2:
 *	- Added support for Ratio Imaging (most figures and tables are still called 'lifetime')
 *	
 *	Version 2.2:
 *	- The positions in the hit list are now correct when no nuclei channel is present in the .lif file
 *	- Improved .rgn file generation, with UUID.
 *	  
 *	Version 2.5:
 *	- Screening: plot the lifetime traces of the hits in a separate plot
 *	- Screening: Color the ROIs according to the traces
 *	
 *	Version 2.6:
 *	- Screening: Hit finding on mean response can now be set to the first n seconds
 *	- Screening: Added fitting traces (not complete yet - need a lot of selection criteria)
 *	
 *	Version 2.7:
 *	- Screening: Hit finding on mean response can now be set to a flexible time window
 *	- Added frame interval to RGB overlay image
 *	  
 *	Version 2.8:
 *	- Added possibility to smooth the traces (and for now, the data as well)
 *	- Other small improvements/changes
 *	
 *	Version 2.9:
 *	- Screening: Possibility to determine hits based on the max difference to the *average* baseline of all cells 
 *	
 *	Version 2.95:
 *	- Screening: Possibility to remove the last frame (sometimes (partially) empty) 
 *	- Various small bug fixes concerning finding hits
 *	- renamed 'lifetime traces image' to 'kymograph'
 *	
 *	Version 2.96-2.97:
 *	- Screening: Save empty graphs when no hits are found (necessary for pooling results in 'Screen_inspect_output_traces.ijm'
 *	- Moved 'remove last frame' option to the script parameters. (fixing a bug where removeLastFrame_boolean did not exist when not running a screen)
 *	
 *	Version 3.0:
 *	- Possibility to re-analyze screening data from saved data (no segmentation - much faster)
 *	- Weighted lifetime image is now also saved (Why not before? Beats me!)
 *	
 *	Version 3.1:
 *	- Segmentation can now be performed on a custom range of frames
 *	
 *	Version 3.2:
 *	- Critical bug fix where the RGB overlay colors did not match the LUT
 *	- Incorporated different Lookup Tables (for the visually challenged)
 *	- Added parameter for the brightness of the RGB overlay
 *	- Fancyfied starting dialog with HTML code
 *	
 *	Version 3.3:
 *	- save Stage Positions in a text file
 *	
 *	Version 3.4:
 *	- Greatly speeded up conversion from labelmap to ROI Manager, using CLIJ2
 *	
 *	Version 3.5:
 *	- Save a table with the cell coordinates for every tile
 *	- Screening: added possibility for absolute maximum lifetime
 *	
 *	Version 3.6:
 *	- cell coordinates are now saved in a single table containing all tiles
 *	
 *	Version 3.7:
 *	- Added possibility to analyze intensity-only data
 *	
 *	Version 3.8:
 *	- Added a parameter for the axis font size of the output plots 
 *	
 *	Version 3.9:
 *	- Finally made traces fitting available
 *	
 *	Version 4.0:
 *	- Added TauContrast support (in 'microscopy')
 *	- Fixed some bugs in the trace fitting
 *	
 *	Version 4.01:
 *	- Added .lif file support for TauContrast (but only single series files)
 *	
 *	Version 4.03:
 *	- Fixed bugs concerning non-timelapse experiments
 *	
 *	Version 4.05:
 *	- Added confidence interval (3 sigma around the mean) in the traces plot
 *	
 *	Version 4.1:
 *	- Added TauContrast and TauSeparation channels in the script parameters
 *	- FDFLIM settings are now hidden, because we don't use them.
 *	
 *	Version 4.2:
 *	- Added xy drift correction
 *	- Fixed typo in Cellpose command, causing 'additional flags' message
 *	- TCSPC fitted / TauSeparation files are now opened using Bio-Formats (because in the latter case it is a .lif file)
 *	- NaNs in the weighted lifetime image are now set to 0, removing artifacts caused by drift / movements (presumably)
 *	
 *	Version 4.4:
 *	- Added Fast FLIM option on exported .tif files containing intensity + lifetime channels, with lifetime values in picoseconds - for Ron Hoebe
 *	- Added multi-series support for .lif files
 *	
 *	Version 4.5:
 *	- histogramBins is now a variable (but not in the script parameters)
 *	- Fixed a bug where the macro would crash if smoothing the RGB image/movie was set to 0 for non-timelapse experiments
 *	
 *	Version 4.6:
 *	- Fixed crashing of the macro when no cells are detected
 *	- Added the Cellpose model as parameter
 *	- Calibration bar position is now a parameter
 *	
 *	Version 4.7:
 *	- Rolled back NaNs to zeros from version 4.2, because it causes artifacts in the RGB overlay when smoothing.
 *	- Finally fixed the number of decimals in the script parameters
 *	- The macro now also outputs a statistics table containing area, intensity (of the projection!), coordinates etc. of the measured cells.
 *	- Scatter plots are created of lifetime vs area and lifetime vs intensity (of the projection!)
 *	- For non-timelapse images also scatter plots lifetime-area and lifetime-intensity are generated.
 *	- Changed default saturatedPixels variable (was set to 3). This changes the cell segmentation.
 *	
 *	Version 4.8:
 *	- Fixed drift correction and smoothing artifact issue by (for smoothing) replacing zeros for NaNs and then padding them with the mean of surrounding pixels
 *	- Separated XY and Time smoothing
 *	- Created parameter for creating scatterplots
 *	- Fixed cell mismatch in scatterplots - colors now match with colors of the traces
 *	 
 *	Version 4.92:
 *	- Rank vector is now also displayed and saved (necessary for selecting traces in the kymograph)
 *	- Showing gidlines in the plot(s) is now a parameter
 *	- Fixed orientation of sorted kymograph (horizontal flip) 
 *	- Fixed formatting of doubles in script parameters
 *	- 4.92b: temporarily reinstated FDFLIM options
 *	
 *	Version 4.93:
 *	- Added possibility to measure intensity in an additional channel
 *	
 *	Version 4.94:
 *	- Fixed a regression issue where filenames were 'doubled' in the output
 *	
 *	Version 4.95:
 *	- fixed some issues regarding curve fitting
 *	
 *	Version 4.96:
 *	- Added the 'first segment, then analyze' option (not a parameter yet). To speed up segmentation all images are loaded and put into a hyperstack. Cellpose runs only once and saves the labelmaps.
 *	  The implementation comes with a lot of if..else statements and should be rewritten at some point.
 *	  Also, lifetimes and intensity images are calculated twice, as well as drift correction; when activted the processing can actually be slower.
 *	  
 *	Version 4.97 + 4.98:
 *	- Added timelapse lifetime histograms
 *	- bug fixes regarding the 'first segment, then analyze' option.
 *	
 *	Version 4.99:
 *	- Included a 'Density plot', a 2D-histogram of 1D-histograms vs time.
 *	
 *	Version 5.01:
 *	- Fixed mistake concerning bidirectional correction
 *	- (Temporarily) fixed pixel size issue by not using Bio-Formats.
 *	- Hits can now be detected on non-timelapse images.
 *	- Hit coordinate files (.rgn) are saved in chucks of 'MaxNrHits'.
 *	
 *	Version 5.02:
 *	- Bidirectional phase mismatch correction does not require images with 2^n x 2^n pixels any more.
 *	
 *	Version 5.03:
 *	- Added possibility to sorts hits on lifetime.
 *	
 *	Version 5.05:
 *	- Fixed hit positions for merged images.
 *	
 *	Version 5.06:
 *	- Use 'NKI Labelmap to ROI Manager' groovy script in stead of IJ1 Macro function (speeds up this step a bit).
 *	
 *	Version 5.07:
 *	- Reverted change of 5.06.
 *	- Included the new 'Lifetime' LUT.
 *	
 *	Version 5.10:
 *	- Fixed bug when measuring additional channel (Exiting batch mode for Cellpose now displays the orignal image)
 *	- A timelapse scatterplot of additional channel intensity vs lifetime is created.
 *	
 *	Version 5.11:
 *	- Reinstated FDFLIM options
 *	- Fixed bug where multiseries fli files were opened multiple times
 *	- Introduced a 1-minute waiting step in case Cellpose fails (This happens with FDFLIM images; for unknown reasons the temp image is not saved unless you wait some time).
 *	
 *	Version 5.12:
 *	- Fixed small bug regarding drift correction
 *	- Changed appearance of RGB overlay script parameter (now radiobuttons)
 *	- Added some comments
 *	
 *	Version 5.20:
 *	- Many small changes and improvement in the screening part; More options, GUI changes, easier calculation, fool-proofing
 *	
 *	Version 5.21:
 *	- Added option to use custom Cellpose model
 *	
 *	Version 5.22 & 5.23:
 *	- Screening: added additional option to restrict hits to traces with a baseline close to the average baseline
 *	- Fixed bug that caused macro to crash on Taucontrast .tif files (image displaying issue)
 *	- Moved nrTilesX and nrTilesY down to the screening part of the script parameter dialog
 *	
 *	Version 5.24:
 *	- TauSeparation .lif files are opened using BioFormats importer (not showing the popup window). TCSPC .tif files are opened with the standard opener.
 *	
 *	Version 5.30:
 *	- Screening: added option to find hits based on the mean intensity in the additional channel + some extras
 *	
 *	Version 5.40:
 *	- Added possibility to load a labelmap from disk
 *	- Added option to skip RGB overlay visualization in order to speed up the processing
 *	- Added option to speed up ROI handling (using BIOP plugin)
 *	- Fixed bug where drift corretion crashed when the data has only 1 frame
 *	
 *	Version 5.50:
 *	- Screening: tile layout is found automatically for images processed with 'stitch_tiles.ijm'.
 *	- Hit positions should now always be correct
 *	- Other small improvements
 *	
 *	Version 5.6:
 *	- Various screening improvements (e.g. calculating averages without NaNs)
 *	- lifetime and intensity measurements are now performed using CLIJ2 (GPU). This makes the processing much faster for experiments with many cells
 *	
 *	Version 5.8:
 *	- Added histograms and kymographs of hits only, and of non-hits (still to be finished)
 *	
 *	Version 5.9:
 *	- Screening parameters are now saved and loaded to/from disk
 *	
 *	Version 5.91:
 *	- Last version that creates tables with cells as columns (becomes very slow for many cells)
 *	
 *	Version 6.00:
 *	- Lifetime table is created transposed, only if the (not displayed) parameter 'generateTables' below is set to true.
 *	  The table is not necessary any more, because all graphs are created from the kymographs.
 *	- Other small improvements
 *	
 *	Version 6.01:
 *	- Enabled negative overlap detection from metadata of the stitched image
 *	
 *	Version 6.03:
 *	- Activated FD-FLIM parameters
 *	- Added a minimum cell intensity selection parameter
 *	
 *	Version 6.10:
 *	- New faster method for measuring intensities without generating a huge table, using ROI Manager
 *	- Generate 'normalized kymograph' and traces graph, where all baselines are moved to the average baseline
 *	- Switched some parameters in the GUI (e.g. Default Frame Interval)
 *	- Added text overlay to the traces graph with the smoothing factor in pixels 
 *	
 *	Version 6.11:
 *	- Included option for screening on the rise time of the response
 *	
 *	Version 6.12:
 *	- Small bug fixes concerning the changes in 6.11
 *	
 *	Version 6.13:
 *	- Bidirectional phase mismatch correction now works on the timelapse instead of only the intensity image
 *	
 *	Version 6.15:
 *	- Moved screening parameters to the screening dialog
 *	- Added feature to measure baseline only (no stimulation/calibration) on all frames
 *	
 *	Version 6.16:
 *	- Fixed hit finding in additional channel. Also works without tables.
 *	- Inverted the parameter to create the RGB overlay image
 *	
 *	Version 6.17:
 *	- Added possibility for screening with only baseline and calibration. To activate, manually enter the stimulation and calibration frames, where the stimulation frame is 0.
 *	- Added a maximum circularity parameter to exclude very round cells
 *	- Copied some screening parameters back to the main script parameter dialog (necessary for sorting kymographs when not screening).
 *	
 *	Version 6.19:
 *	- 6.18: Some bug fixes introduced in version 6.17
 *	- 6.19: Added a time column (at the end) in the lifetime table 
 *	
 *	Version 6.20:
 *	- Possibility to manually 'segment' cells or regions. To do so, select 'manual segmentation' in the Cellpose model field
 *	
 *	Version 6.21:
 *	- Improved GUI with headings
 *	
 *	Version 6.31:
 *	- Adapted for new Cellpose Fiji wrapper
 *	- Possibility to generate random hits in the screening (convenient for testing purposes)
 *	- Added a boolean for loading the labelmap
 *	
 *	Version 6.32:
 *	- Implemented support for Lambert Instruments FD-FLIM recordings with SPAD Toggel camera.
 *	- Added a choice of using lifetime from phase or modulation
 *	
 *	Version 6.40:
 *	- Automatic Cellpose env_path and env_type detection
 *	
 *	Version 6.41:
 *	- Added baseline_calibration_difference_number_high parameter
 *	
 *	Version 6.48:
 *	- Screening: hit selections criteria can now be combined with AND or OR logic
 *	- Fixed a small bug where positions would be wrong in re-analysis when no RGB image was created.
 *	- Added 'fraction of max response' as hit searching criterium
 *	
 *	Version 6.50:
 *	- Synchonization between script parameters and dialog parameters (stimulationSensitivity and manualStimCalFrames)
 *	- Hit detection region is now visualized on the hit only plot
 *	
 *	Version 6.50:
 *	- Fixed small mistake (update of 6.41 didn't work properly).
 *
 *	Version 6.60:
 *	- !! Changed the calculation of the intensity image for TCSPC images !! 
 *	  This was previously done by adding the two components, but since the components represent amplitudes of the two lifetimes 
 *	  it actually should be weighted by the lifetimes. This is now implemented correctly, but could lead to slightly different results than in earlier versions.
 *	  
 *	Version 6.61:
 *	- Fixed a bug where the macro would crash when in re-analysis no hits were found
 *	
 *	Version 6.64:
 *	- Added possibility to segment cells on a user-defined channel
 *	- Fixed bug where the macro would crash in baseline only mode
 *	
 *	Version 6.65
 *	- Corrected selection of the baseline, response and calibration parts in the kymograph (crashed when the image size was smaller than the number of frames) 
 *	- Fixed kymograph and density plot of non-hits
 *	
 *	To do:
 *	- Separate segmentation from measurements; create possibility to run single, or both.
 *		- separate macros? macro tools?
 *	- Recombine with StarDist / Cellpose nuclei
 *	- Check (and fix) behaviour with nuclei as second channel
 *	- Check and update hit finding scheme
 *	- possibility to change min and max lifetime / other visualization options
 *	- 'Tracking' by using 3D Cellpose / Trackmate / BIOP Trackmate on ROIs (check earlier macro for this)?
 *	- Hide parameters in other dialogs, store using preferences or hidden script parameters
 *	
 *	- Clean up code
 *	
 */

version = 6.65;

#@ String	file_and_image_message	(value="<html><p style='font-size:12px; color:#3366cc; font-weight:bold'>Input file settings</p></html>", visibility="MESSAGE") 
#@ String	microscope				(label = "Input file microscopy type", choices={"Confocal TCSPC / TauSeparation", "Fast FLIM", "TauContrast", "Frequency Domain FLIM", "Ratio Imaging", "Intensity only"}, style="listBox")
#@ File[]	file_list				(label = "Input files (exported .tif, .lif or .fli files)", style="File")
#@ File		output					(label = "Output folder", style = "directory") 
#@ Boolean	segmentFirstAnalyzeLater(label="Segment all images before analyzing (faster, but dimensions must be the same)", value=true)

#@ String 	confocal_message		(value="<html><p style='font-size:12px; color:#3366cc; font-weight:bold'>FLIM and image settings</p></html>", visibility="MESSAGE")
#@ Integer	intensityChannel		(label="TCSPC fitted / TauContrast / TauSeparation (first) channel", value=1, description="TauContrast saves with an additional (hidden) intensity channel. Choose the LAS X channel that shows the TauContrast lifetime.")
#@ Double	tau1					(label="Lifetime component 1 (ns)", value=0.6, description="For fitted TCSPC or TauSeparation data", style="format:0.00")
#@ Double	tau2					(label="Lifetime component 2 (ns)", value=3.4, description="For fitted TCSPC or TauSeparation data", style="format:0.00")
#@ Double	default_frameInterval	(label="Default frame interval (s) (if not found)", value = 5, style="format:#.0")
#@ Boolean	correctBidirect_boolean	(label="Correct bidirectional phase mismatch",value=false)
#@ Boolean	correctDrift_boolean	(label="Correct xy drift (using intensity)",value=false)
#@ String	registration_against	(label = "Registration against [Drift correction]", choices={"First frame","Last frame","Previous frame"}, style="listBox")
#@ Boolean	edge_detect				(label="Edge-detect before registration? [Drift correction]", value=false)
#@ Integer	additionalChannel		(label="Measure intensity in additional channel (-1 if N/A)", value=-1, min=-1)
#@ Boolean	removeLastFrame_boolean	(label="Remove last frame (sometimes (partially) empty)", value=false)
#@ Integer	nucleiChannel			(label="Nuclei channel in the .lif file (-1 if not present)", value=-1)

//#@ String	FDFLIM_message			(value="<html><p style='font-size:12px; color:#3366cc; font-weight:bold'>Frequency Domain FLIM settings</p></html>", visibility="MESSAGE")
//#@ File		reference				(visibility=INVISIBLE, label = "Reference file", style = "file", required=false)
//#@ Integer	phases					(label="Number of phases", value = 12)
//#@ Integer	freq					(label="Frequency (MHz)", value = 40)
//#@ Double	tau_ref					(label="Lifetime of the reference", value = 3.93)
//#@ String	tauPhiOrTauMod			(label = "Get lifetime from", choices={"phase", "modulation"}, style="radioButtonHorizontal")
//#@ Boolean	FDFLIM_SPAD_camera		(label="SPAD camera?", value = false)

#@ String	segmentation_message	(value="<html><p style='font-size:12px; color:#3366cc; font-weight:bold'>Cell segmentation settings</p></html>", visibility="MESSAGE")
#@ Integer	segmentationChannel		(label="Segment cells on channel (-1 if no dedicated channel is present)", value=-1)
#@ Boolean	load_labelmap_boolean	(label="Load labelmap from disk instead of Cellpose segmentation?", value=false)
#@ File		labelmapPath			(label = "Labelmap path", style="File", required=false, description="Load segmented labels from disk instead of running Cellpose segmentation.")
#@ String	CellposeModel			(label = "Cellpose model", choices={"cyto2_cp3","cyto3","nuclei","custom", "manual segmentation"}, style="listBox", value="cyto3")
#@ File		cellposeModelPath		(label = "Custom Cellpose model path", style="File", required=false)
#@ Boolean	equalize_contrast_cp	(label="Enhance contrast before cell segmentation (Gamma correction)", value=false)
#@ Integer	CellposeStartFrame		(label="Start frame for segmentation (-1 for the full timelapse)", value=-1)
#@ Integer	CellposeEndFrame		(label="End frame for segmentation (-1 for the full timelapse)", value=-1)
//#@ Integer	upSampleFactor		(label="Image upscale factor (for 10X images)", value=1, visibility="INVISIBLE")
#@ Integer	CellposeDiameter		(label="Cell diameter (0 for automatic)", value=0)
#@ Double	CellposeFlowThreshold	(label="Cell flow error threshold (Cellpose default=0.4, higher -> more cells)", style="scroll bar", value=1, min=0, max=3, stepSize=0.1)
#@ Double	CellposeProbability		(label="Cell probability threshold (default=0, higher -> smaller cell area)", style="scroll bar", value=0, min=-6, max=6, stepSize=0.25)
#@ Integer	minCellSize				(label="Minimum cell size (pixels)", value=50, min=0)
#@ Double	maxCircularity			(label="Maximum cell circularity [0-1]", value=1.0, min=0.0, max=1.0)
#@ Integer	minCellBrightness		(label="Minimum cell intensity (gray values)", value=0, min=0)
#@ Double	StarDistProbThreshold	(label="Nuclei Probability threshold [StarDist]", style="scroll bar", value=0.4, min=0.0, max=1.0, stepSize=0.1)

#@ String	display_message			(value="<html><p style='font-size:12px; color:#3366cc; font-weight:bold'>Output and visualization settings</p></html>", visibility="MESSAGE")
#@ String	lut						(label = "Lookup table", choices={"Lifetime", "Turbo", "Fire", "mpl-viridis", "mpl-plasma", "mpl-viridis", "phase", "glow", "Grays"}, style="listBox", value="Turbo")
#@ Double	minLifetime				(label="Min. displayed lifetime(ns)", value=2.0, style="format:0.0")
#@ Double	maxLifetime				(label="Max. displayed lifetime(ns)", value=3.4, style="format:0.0")
#@ Double	smoothRadiusTraces		(label="Smooth traces (graphs, scatterplots, hit detection) with radius", value=0.0, min=0.0, style="format:0.0")
#@ Double	smoothRadiusOverlayXY	(label="Smooth lifetime in overlay movie with radius (x and y)", value=1.0, min=0.0, style="format:0.0")
#@ Double	smoothRadiusOverlayTime	(label="Smooth lifetime in overlay movie with radius (time)", value=1.0, min=0.0, style="format:0.0")
#@ String	overlayMovie			(label = "Create RGB lifetime overlay on", choices={"intensity image", "intensity movie"}, style="radioButtonHorizontal")
#@ Double	RGB_brightness			(label="Brightness of RGB overlay [0-5]", value=3, min=0, max=5, style="format:0.0")
#@ String	calibrationBarPosition	(label = "Calibration bar position", choices={"Upper Right","Upper Left","Lower Right","Lower Left"}, style="listBox", value="Upper Right")
#@ Boolean	displayGrid				(label="Display grid lines in plot", value=true)
#@ Integer	axisFontSize			(label="Output plot axis font size", value=18)
#@ Boolean	createHistogram			(label="Create timelapse histogram", value=true)
#@ Boolean	createScatterPlots		(label="Create scatterplots (lifetime vs intensity/area)", value=true)
#@ Boolean	createRGBOverlay		(label="Create RGB overlay", value=false)
#@ Boolean	generateTables			(label="Generate lifetime table", value=true)

#@ String	other_options_message	(value="<html><p style='font-size:12px; color:#3366cc; font-weight:bold'>Screening settings</p></html>", visibility="MESSAGE")
#@ Boolean	runScreen				(label="Activate screening: find hits and write the coordinates to a file", value=false)
//#@ File	lif_file				(label = "Input .lif file (for stage coordinates and nuclei)", style="File", default="-")
#@ Double	stimulationSensitivity	(label="Sensitivity in detecting stimulation/calibration frames (#stdDevs of 2nd derivative)", value=1, min=0.0, style="format:0.0")
#@ String	manualStimCalFrames		(label="Manual stimulation and calibration frames with (comma-separated - empty for automatic)", value="")
//#@ Boolean	speedup				(label="Faster ROI processing", value=false)
#@ Boolean	debugMode				(label="Debug mode", value=false)

print("\\Clear");
if(additionalChannel == 0) additionalChannel = -1;
speedup = false; //Faster ROI processing using BIOP labels to ROIs plugin, but messes up composite ROIs. Therefore set to false by default. 

risetime_fraction = 0.75;
useColoredHitOverlay = false;
histogramBins = 50;		//Number of bins in the histograms
saturatedPixels = 0.35;	//Contrast settings before applying the LUT for Cellpose - Some saturation gives better results if there are also dim cells
RGB_brightness = pow(10, (0.5*(RGB_brightness-2)) );
sigma_CC = 4;			//Blur Cross Correlation with sigma (pixels)
show_CC = false;		//Show Cross Correlation image
subpixel = false;		//If enabled, the shift positions are calculated using the center of mass of (a crop of) the cross correlation images, instead of the maximum pixel. This allows for subpixel accuracy, but it is not always better.
displayConfidenceInterval = false;	//Display red lines in the traces plot depicting the confidence interval
confidence_interval_sigma = 2;
createNormalizedKymograph = false;
upSampleFactor = 1;
lifetimeChannel = intensityChannel + 1;
reanalyze_boolean = false;

//Cellpose parameters
List.setCommands;
if (List.get("Cellpose ...")!="") newCellposeWrapper = true;
else if (List.get("Cellpose Advanced")!="") newCellposeWrapper = false;
else exit("ERROR: Cellpose wrapper not found!");
env_path = getPref("Packages.ch.epfl.biop.wrappers.cellpose.ij2commands.Cellpose", "env_path");
env_type = getPref("Packages.ch.epfl.biop.wrappers.cellpose.ij2commands.Cellpose", "env_type");
if(env_type == "<null>") env_type = "conda";	//Default is conda, but returns <null>
print("Cellpose environment ("+env_type+") found at "+env_path);

output = output + File.separator;
if(!File.exists(output)) {
	print("Output directory "+output+" does not exist. Creating it.");
	File.makeDirectory(output);
}

upSampleFactor = 1;
var pixelWidth;
var tileSizeX;
var tileSizeY;
var stimulationFrame = -1;
var calibrationFrame = -1;
var singleFrame = false;
var nrHits;
var baselineOnly = false;
var calibrationOnly = false;
maxNrHits = 1200;		//Maximum number of hits in a single .rgn file

hitList = "Hit List";
var responseTimeStart = 0;				//Required for drawing hit finding region on graph
var hit_avg_response_time_window = 0;	//Required for drawing hit finding region on graph
var baselineLength = 0;

saveSettings();

run("Colors...", "foreground=white background=black");
run("Conversions...", " ");

run("CLIJ2 Macro Extensions", "cl_device=");
Ext.CLIJ2_clear();

windowList = getList("window.titles");
for(i=0; i<windowList.length; i++) {
	selectWindow(windowList[i]);
	run("Close");
}
//close_windows("Hit");
//close_windows("Cell_coordinates");
//close_windows("Lifetime_Data");
//close_windows("Cell_statistics");

//Create dialog for screening. First load previous settings from disk.
if(runScreen == true) {

	//Set these two (copy from script parameters)
	call("ij.Prefs.set", "stimulation.Sensitivity", stimulationSensitivity);
	call("ij.Prefs.set", "manual.Stim.Cal.Frames", manualStimCalFrames);
	
	//get Prefs (not ordered)
	lif_file = 										call("ij.Prefs.get", "lifFile.Path", "");
	stimulationSensitivity =						call("ij.Prefs.get", "stimulation.Sensitivity", 2.0);
	manualStimCalFrames = 							call("ij.Prefs.get", "manual.Stim.Cal.Frames", "");
	absoluteOrRelativeOption =						call("ij.Prefs.get", "absolute.Or.Relative.Option", "response lifetime difference with baseline");
	baseline_avg_baseline_difference_boolean =		call("ij.Prefs.get", "baseline.avg.baseline.difference.boolean", true);
	baseline_avg_baseline_difference_number =		call("ij.Prefs.get", "baseline.avg.baseline.difference.number", 0.2);
	baseline_calibration_difference_boolean =		call("ij.Prefs.get", "baseline.calibration.difference.boolean", true);
	baseline_calibration_difference_number_low =	call("ij.Prefs.get", "baseline.calibration.difference.number.low", 0.5);
	baseline_calibration_difference_number_high =	call("ij.Prefs.get", "baseline.calibration.difference.number.high", 1.0);
	fit_equation =									call("ij.Prefs.get", "fit.equation", "y = a + b/(1+exp(-d*(x-c)))");
	fit_traces =									call("ij.Prefs.get", "fit.traces", false);
	initialGuesses =								call("ij.Prefs.get", "initial.Guesses", "");
	hit_find_logic =								call("ij.Prefs.get", "hit.find.logic", "OR");
	hit_additional_channel_boolean =				call("ij.Prefs.get", "hit.additional.channel.boolean", false);
	hit_additional_channel_choice =					call("ij.Prefs.get", "hit.additional.channel.choice", "higher");
	hit_additional_channel_number1 =				call("ij.Prefs.get", "hit.additional.channel.number1", 10);
	hit_additional_channel_number2 =				call("ij.Prefs.get", "hit.additional.channel.number2", 65535);
	hit_additional_channel_timing =					call("ij.Prefs.get", "hit.additional.channel.timing", "Average intensity");
	hit_avg_response_boolean =						call("ij.Prefs.get", "hit.avg.response.boolean", true);
	hit_avg_response_choice =						call("ij.Prefs.get", "hit.avg.response.choice", "lower");
	hit_avg_response_number =						call("ij.Prefs.get", "hit.avg.response.number", 0.3);
	hit_avg_response_time_margin =					call("ij.Prefs.get", "hit.avg.response.time.margin", 0);
	hit_avg_response_time_window =					call("ij.Prefs.get", "hit.avg.response.time.window", 3);
	hit_baseline_boolean =							call("ij.Prefs.get", "hit.baseline.boolean", false);
	hit_baseline_choice =							call("ij.Prefs.get", "hit.baseline.choice", "higher");
	hit_baseline_number =							call("ij.Prefs.get", "hit.baseline.number", 2.8);
	hit_max_abs_response_boolean =					call("ij.Prefs.get", "hit.max.abs.response.boolean", false);
	hit_max_abs_response_choice =					call("ij.Prefs.get", "hit.max.abs.response.choice", "lower");
	hit_max_abs_response_number =					call("ij.Prefs.get", "hit.max.abs.response.number", 2.8);
	hit_max_response_boolean =						call("ij.Prefs.get", "hit.max.response.boolean", false);
	hit_max_response_choice =						call("ij.Prefs.get", "hit.max.response.choice", "lower");
	hit_max_response_number =						call("ij.Prefs.get", "hit.max.response.number", 0.2);
	hit_max_response_to_avg_baseline_boolean =		call("ij.Prefs.get", "hit.max.response.to.avg.baseline.boolean", false);
	hit_risetime_boolean =							call("ij.Prefs.get", "hit.risetime.boolean", false);
	hit_risetime_choice =							call("ij.Prefs.get", "hit.risetime.choice", "slower");
	hit_risetime_number =							call("ij.Prefs.get", "hit.risetime.number", 1);
	hit_rrr_boolean =								call("ij.Prefs.get", "hit.rrr.boolean", false);
	hit_rrr_choice =								call("ij.Prefs.get", "hit.rrr.choice", "lower");
	hit_rrr_number =								call("ij.Prefs.get", "hit.rrr.number", 0.5);
	maxNrHits =										call("ij.Prefs.get", "max.Nr.Hits", 1000);
	meanResponseTimingOption =						call("ij.Prefs.get", "mean.Response.Timing.Option", "before calibration");
	nrTilesX =										call("ij.Prefs.get", "nr.Tiles.X", 1);
	nrTilesY =										call("ij.Prefs.get", "nr.Tiles.Y", 1);
	tileOverlap =									call("ij.Prefs.get", "tile.Overlap", 0);
	reanalyze_boolean =								call("ij.Prefs.get", "reanalyze.boolean", false);
	sort_hits_column =								call("ij.Prefs.get", "sort.hits.column", "");
	sort_hits_direction = 							call("ij.Prefs.get", "sort.hits.direction", "highest first");
	topNHits =										call("ij.Prefs.get", "top.N.Hits", 1000);
	generateRandomHits_boolean =					call("ij.Prefs.get", "generate.RandomHits.boolean", false);
	nrOfRandomHits =								call("ij.Prefs.get", "nr.Of.Random.Hits", 0);

	Dialog.createNonBlocking("Settings for detecting hits in the screen.");
	Dialog.addMessage("Response time is the time between stimulation and calibration (if any). 'Response difference' is the mean amplitude from the cell's baseline.\nDetect hits as (logical OR, sorting on the last selected item).");

	Dialog.addCheckbox("Only re-analyze data (reads output files of already analyzed data from disk)", reanalyze_boolean);
	Dialog.addFile(".lif file for stage coordinates (if not found in stitched file)", lif_file);
	Dialog.addMessage("Manual stimulation/calibration frames (comma-separated):");
	Dialog.addString("     empty for automatic detection | '0' for baseline only mode", manualStimCalFrames, 25);
	Dialog.addNumber("Sensitivity in detecting stimulation & calibration frames", stimulationSensitivity, 1, 3, "(# stdDevs)");
	Dialog.addMessage("");
	
	Dialog.addMessage("Hit detection criteria:");

	Dialog.addRadioButtonGroup("    Logic", newArray("OR", "AND"), 1, 2, hit_find_logic);
	
	Dialog.setInsets(10, 20, 0);
	Dialog.addCheckbox("baseline lifetime", hit_baseline_boolean);
	Dialog.addToSameRow();
	Dialog.addChoice("", newArray("lower","higher"), hit_baseline_choice);
	Dialog.addToSameRow();
	Dialog.addNumber("than", hit_baseline_number, 2, 5, "ns");

	Dialog.setInsets(0, 20, 10);
	Dialog.addCheckbox("", hit_avg_response_boolean);
//	Dialog.addToSameRow();
	Dialog.setInsets(-35, -298, 0);
	Dialog.addChoice("", newArray("mean absolute response lifetime", "mean response lifetime difference with baseline", "fraction of max response minus baseline"), absoluteOrRelativeOption);
	Dialog.addToSameRow();
	Dialog.addNumber("in a", hit_avg_response_time_window, 0, 3, "-frames window (-1 for full range)");
	Dialog.addToSameRow();
	Dialog.addChoice("", newArray("after stimulation", "before calibration"), meanResponseTimingOption);
	Dialog.addToSameRow();
	Dialog.addNumber("with margin", hit_avg_response_time_margin, 0, 3, "frame(s)");
	Dialog.addToSameRow();
	Dialog.addChoice("", newArray("lower","higher"), hit_avg_response_choice);
	Dialog.addToSameRow();
	Dialog.addNumber("than", hit_avg_response_number, 2, 5, "ns");
	
	Dialog.setInsets(0, 20, 5);
	Dialog.addCheckbox("maximum absolute response lifetime", hit_max_abs_response_boolean);
	Dialog.addToSameRow();
	Dialog.addChoice("", newArray("lower","higher"), hit_max_abs_response_choice);
	Dialog.addToSameRow();
	Dialog.addNumber("than", hit_max_abs_response_number, 2, 5, "ns");
	
	Dialog.addCheckbox("maximum response response lifetime difference with baseline", hit_max_response_boolean);
	Dialog.addToSameRow();
	Dialog.addChoice("", newArray("lower","higher"), hit_max_response_choice);
	Dialog.addToSameRow();
	Dialog.addNumber("than", hit_max_response_number, 2, 5, "ns");
	Dialog.addToSameRow();	
	Dialog.addCheckbox("relative to average baseline of all cells", hit_max_response_to_avg_baseline_boolean);

	Dialog.addCheckbox("Rise time of the response", hit_risetime_boolean);
	Dialog.addToSameRow();
	Dialog.addChoice("", newArray("slower","faster"), hit_risetime_choice);
	Dialog.addToSameRow();
	Dialog.addNumber("than", hit_risetime_number, 0, 5, "frames");

	Dialog.addCheckbox("Rapid Response Ratio © (2nd half / 1st half of the response)", hit_rrr_boolean);
	Dialog.addToSameRow();
	Dialog.addChoice("", newArray("lower","higher"), hit_rrr_choice);
	Dialog.addToSameRow();
	Dialog.addNumber("than", hit_rrr_number, 2, 5, "");
	
	if(additionalChannel != -1) {
		Dialog.addCheckbox("Intensity in the selected additional channel ("+additionalChannel+")", hit_additional_channel_boolean);
		Dialog.addToSameRow();
		Dialog.addChoice("", newArray("lower","higher","between"), hit_additional_channel_choice);
		Dialog.addToSameRow();
		Dialog.addNumber("(than)", hit_additional_channel_number1, 1, 7, "");
		Dialog.addToSameRow();
		Dialog.addNumber("(and", hit_additional_channel_number2, 1, 7, ")");
		Dialog.setInsets(-25, 00, 0);
		Dialog.addChoice("", newArray("Average intensity", "Intensity in the first frame"), hit_additional_channel_timing);
	}
	Dialog.addMessage("AND (extra conditions):");
	Dialog.addCheckbox("baseline to calibration lifetime difference is between", baseline_calibration_difference_boolean);
	Dialog.addToSameRow();
	Dialog.addNumber("", baseline_calibration_difference_number_low, 2, 5, "ns");
	Dialog.addToSameRow();
	Dialog.addNumber("and", baseline_calibration_difference_number_high, 2, 5, "ns");
	Dialog.addCheckbox("(cell baseline - average baseline) is lower than", baseline_avg_baseline_difference_boolean);
	Dialog.addToSameRow();	
	Dialog.addNumber("", baseline_avg_baseline_difference_number, 2, 5, "ns");

	Dialog.addMessage("--------------------------------------------------------------------------------------------------------------------");

	Dialog.addNumber("Number of tiles X (overruled if found in metadata)", nrTilesX, 0, 4, "");
	Dialog.addNumber("Number of tiles Y (overruled if found in metadata)", nrTilesY, 0, 4, "");
	Dialog.addNumber("Tile overlap (%) (overruled if found in metadata)", tileOverlap, 1, 4, "");

	Dialog.addMessage("Output options:");
	Dialog.addChoice("Sort hits on criterium", newArray("mean baseline","mean response","max response","max response diff to avg baseline","rise time (frames)","rapid response ratio","mean additional channel intensity","do not sort"), sort_hits_column);
	Dialog.addChoice("Sort direction", newArray("highest first", "lowest first"), sort_hits_direction);
	Dialog.addNumber("        Generate .rgn files in chunks of ", maxNrHits, 0, 5, "hits");
	Dialog.addNumber("        Generate an additional .rgn file with top ", topNHits, 0, 5, "hits");

	Dialog.addMessage("Curve fitting options:");
	Dialog.addCheckbox("Fit traces (from stimulation to calibration, if present)", fit_traces);
	Dialog.addString("Fit equation (max 5 parameters a,b,c,d,e)", fit_equation, 30);
//	logistics curve: y = a + b/(1+exp(-d*(x-c)))
//	initial guesses for logistics curve: 2.3, 0.7, 200*frameInterval, -0.01
//	Dialog.addString("Fit equation (max 5 parameters a,b,c,d)", "y = a*exp(-(x)/(b*0.69315)) + c", 30);
	Dialog.addString("Initial parameter guesses, separated by spaces (optional)", initialGuesses, 30);

	Dialog.addMessage("--------------------------------------------------------------------------------------------------------------------");

	Dialog.addCheckbox("Do not screen. Instead generate random hits from the segmented cells (nr)", generateRandomHits_boolean);
	Dialog.addToSameRow();
	Dialog.addNumber("", nrOfRandomHits, 0, 5, "");

	Dialog.show();
	
	
	reanalyze_boolean = Dialog.getCheckbox();

	lif_file = Dialog.getString();
	lif_file = replace(lif_file, "\\", "/");
	manualStimCalFrames = Dialog.getString();
	stimulationSensitivity = Dialog.getNumber();

	hit_find_logic = Dialog.getRadioButton();
	
	hit_baseline_boolean = Dialog.getCheckbox();
	hit_baseline_choice = Dialog.getChoice();
	hit_baseline_number = Dialog.getNumber();

	hit_avg_response_boolean = Dialog.getCheckbox();
	absoluteOrRelativeOption = Dialog.getChoice();
	hit_avg_response_time_window = Dialog.getNumber();
	meanResponseTimingOption = Dialog.getChoice();
	hit_avg_response_time_margin = Dialog.getNumber();
	hit_avg_response_choice = Dialog.getChoice();
	hit_avg_response_number = Dialog.getNumber();
	
	hit_max_abs_response_boolean = Dialog.getCheckbox();
	hit_max_abs_response_choice = Dialog.getChoice();
	hit_max_abs_response_number = Dialog.getNumber();
	
	hit_max_response_boolean = Dialog.getCheckbox();
	hit_max_response_choice = Dialog.getChoice();
	hit_max_response_number = Dialog.getNumber();
	hit_max_response_to_avg_baseline_boolean = Dialog.getCheckbox();

	hit_risetime_boolean =  Dialog.getCheckbox();
	hit_risetime_choice =  Dialog.getChoice();
	hit_risetime_number =  Dialog.getNumber();

	hit_rrr_boolean = Dialog.getCheckbox();
	hit_rrr_choice = Dialog.getChoice();
	hit_rrr_number = Dialog.getNumber();

	if(additionalChannel != -1) {
		hit_additional_channel_boolean = Dialog.getCheckbox();
		hit_additional_channel_choice = Dialog.getChoice();
		hit_additional_channel_number1 = Dialog.getNumber();
		hit_additional_channel_number2 = Dialog.getNumber();
		hit_additional_channel_timing = Dialog.getChoice();
	}
	else hit_additional_channel_boolean = false;
	
	baseline_calibration_difference_boolean = Dialog.getCheckbox();
	baseline_calibration_difference_number_low = Dialog.getNumber();
	baseline_calibration_difference_number_high = Dialog.getNumber();
	baseline_avg_baseline_difference_boolean = Dialog.getCheckbox();
	baseline_avg_baseline_difference_number = Dialog.getNumber();

	nrTilesX = Dialog.getNumber();
	nrTilesY = Dialog.getNumber();
	tileOverlap = Dialog.getNumber()/100;
	
	sort_hits_column = Dialog.getChoice();
	sort_hits_direction = Dialog.getChoice();
	maxNrHits = Dialog.getNumber();
	topNHits = Dialog.getNumber();
	
	fit_traces = Dialog.getCheckbox();
	fit_equation = Dialog.getString();
	initialGuesses = Dialog.getString();
	generateRandomHits_boolean = Dialog.getCheckbox();
	nrOfRandomHits = Dialog.getNumber();


	//Set Prefs
	call("ij.Prefs.set", "lifFile.Path", lif_file);
	call("ij.Prefs.set", "stimulation.Sensitivity", stimulationSensitivity);
	call("ij.Prefs.set", "manual.Stim.Cal.Frames", manualStimCalFrames);
	call("ij.Prefs.set", "absolute.Or.Relative.Option", absoluteOrRelativeOption);
	call("ij.Prefs.set", "baseline.avg.baseline.difference.boolean", baseline_avg_baseline_difference_boolean);
	call("ij.Prefs.set", "baseline.avg.baseline.difference.number", baseline_avg_baseline_difference_number);
	call("ij.Prefs.set", "baseline.calibration.difference.boolean", baseline_calibration_difference_boolean);
	call("ij.Prefs.set", "baseline.calibration.difference.number.low", baseline_calibration_difference_number_low);
	call("ij.Prefs.set", "baseline.calibration.difference.number.high", baseline_calibration_difference_number_high);
	call("ij.Prefs.set", "fit.equation", fit_equation);
	call("ij.Prefs.set", "fit.traces", fit_traces);
	call("ij.Prefs.set", "initial.Guesses", initialGuesses);
	call("ij.Prefs.set", "hit.additional.channel.boolean", hit_additional_channel_boolean);
	call("ij.Prefs.set", "hit.additional.channel.choice", hit_additional_channel_choice);
	call("ij.Prefs.set", "hit.additional.channel.number1", hit_additional_channel_number1);
	call("ij.Prefs.set", "hit.additional.channel.number2", hit_additional_channel_number2);
	call("ij.Prefs.set", "hit.additional.channel.timing", hit_additional_channel_timing);
	call("ij.Prefs.set", "hit.find.logic", hit_find_logic);
	call("ij.Prefs.set", "hit.avg.response.boolean", hit_avg_response_boolean);
	call("ij.Prefs.set", "hit.avg.response.choice", hit_avg_response_choice);
	call("ij.Prefs.set", "hit.avg.response.number", hit_avg_response_number);
	call("ij.Prefs.set", "hit.avg.response.time.margin", hit_avg_response_time_margin);
	call("ij.Prefs.set", "hit.avg.response.time.window", hit_avg_response_time_window);
	call("ij.Prefs.set", "hit.baseline.boolean", hit_baseline_boolean);
	call("ij.Prefs.set", "hit.baseline.choice", hit_baseline_choice);
	call("ij.Prefs.set", "hit.baseline.number", hit_baseline_number);
	call("ij.Prefs.set", "hit.max.abs.response.boolean", hit_max_abs_response_boolean);
	call("ij.Prefs.set", "hit.max.abs.response.choice", hit_max_abs_response_choice);
	call("ij.Prefs.set", "hit.max.abs.response.number", hit_max_abs_response_number);
	call("ij.Prefs.set", "hit.max.response.boolean", hit_max_response_boolean);
	call("ij.Prefs.set", "hit.max.response.choice", hit_max_response_choice);
	call("ij.Prefs.set", "hit.max.response.number", hit_max_response_number);
	call("ij.Prefs.set", "hit.max.response.to.avg.baseline.boolean", hit_max_response_to_avg_baseline_boolean);
	call("ij.Prefs.set", "hit.risetime.boolean", hit_risetime_boolean);
	call("ij.Prefs.set", "hit.risetime.choice", hit_risetime_choice);
	call("ij.Prefs.set", "hit.risetime.number", hit_risetime_number);
	call("ij.Prefs.set", "hit.rrr.boolean", hit_rrr_boolean);
	call("ij.Prefs.set", "hit.rrr.choice", hit_rrr_choice);
	call("ij.Prefs.set", "hit.rrr.number", hit_rrr_number);
	call("ij.Prefs.set", "max.Nr.Hits", maxNrHits);
	call("ij.Prefs.set", "mean.Response.Timing.Option", meanResponseTimingOption);
	call("ij.Prefs.set", "nr.Tiles.X", nrTilesX);
	call("ij.Prefs.set", "nr.Tiles.Y", nrTilesY);
	call("ij.Prefs.set", "tile.Overlap", tileOverlap);
	call("ij.Prefs.set", "reanalyze.boolean", reanalyze_boolean);
	call("ij.Prefs.set", "sort.hits.column", sort_hits_column);
	call("ij.Prefs.set", "sort.hits.direction", sort_hits_direction);
	call("ij.Prefs.set", "top.N.Hits", topNHits);
	call("ij.Prefs.set", "generate.RandomHits.boolean", generateRandomHits_boolean);
	call("ij.Prefs.set", "nr.Of.Random.Hits", nrOfRandomHits);

	setScriptParameterValue("stimulationSensitivity", stimulationSensitivity);
	setScriptParameterValue("manualStimCalFrames", manualStimCalFrames);

	if(!matches(fit_equation, ".*b.*")) nrFitParameters = 1;
	else if(!matches(fit_equation, ".*c.*")) nrFitParameters = 2;
	else if(!matches(fit_equation, ".*d.*")) nrFitParameters = 3;
//	else if(!matches(fit_equation, ".*e.*")) nrFitParameters = 4;
	else nrFitParameters = 4;

	initialGuesses = split(initialGuesses, " ");
	if (initialGuesses.length != 0 && initialGuesses.length != nrFitParameters) exit("The number of inital guesses ("+initialGuesses.length+") must be equal to the number fit parameters ("+nrFitParameters+"), or zero. Exiting macro.");

	//TO DO: Print and/or save screening parameters.
	//- Cannot do this automatically for all settings, so they should be copied to a list, or an array and use String.join(array, delimiter).

}





//Skip files that are not .fli, .tif or .lif
Array.sort(file_list);
for(i=0; i<file_list.length; i++) {
	if( (!endsWith(file_list[i], ".fli")) && (!endsWith(file_list[i], ".tif")) && (!endsWith(file_list[i], ".lif"))) {
		print(file_list[i] + " is not a valid input file - skipping it. Only .lif, .tif and .fli files are supported.");
		file_list = Array.deleteIndex(file_list, i);
		i--;
	}
}

//Get stage coordinates for separate tif files (not stitched). Do this once here. For a single file, it is done later when the file is open.
if(runScreen == true && file_list.length>1) {
	//Get stage coordinates of all the tiles. N.B. Counting always starts at the first tile in the .lif file.
	oldLogWindow = getInfo("log");
	print("\\Clear");
	run("NKI get stage coordinates to log window", "file=["+lif_file+"], nrtiles="+file_list.length);	//Run a Jython script to print the stage coordinates to the log window
	logWindow = getInfo("log");
	stagePositions = split(logWindow, "\n");
	print("\\Clear");
	print(oldLogWindow);
	print("Found stage coordinates");
	tiles = file_list.length;
}
print(file_list.length + " files/tiles to analyze.");

//Loop over all files. Open intensity and lifetime from .tif file, and open nuclei from the .lif file, if present.
alreadySegmented = false;
if(file_list.length == 1) segmentFirstAnalyzeLater = false;
for (f = 0; f < file_list.length; f++) {
	if (microscope == "Confocal TCSPC / TauSeparation" || microscope == "TauContrast") {
		run("Bio-Formats Macro Extensions");
		Ext.setId(file_list[f]);
		Ext.getSeriesCount(nr_series);
	}
	else nr_series = 1;
	
	//Loop over all series
	for(s = 0; s < nr_series; s++) {
		run("Close All");
		if(isOpen("Lifetime_Data")) {
			selectWindow("Lifetime_Data");
			run("Close");
		}
		if(isOpen("Intensity_table_ch"+additionalChannel)) {
			selectWindow("Intensity_table_ch"+additionalChannel);
			run("Close");
		}
		if(isOpen("ROI Manager")) {
			selectWindow("ROI Manager");
			run("Close");
		}
		setBatchMode(false);
		run("ROI Manager...");
		setBatchMode(true);
		roiManager("reset");
		roiManager("Set Color", "grays");
		roiManager("Set Line Width", 0);
	
		setBatchMode(true);

		run("Bio-Formats Macro Extensions");			//Run extensions again, because only one macro extension can be loaded
		Ext.setId(file_list[f]);
		Ext.setSeries(s);								//Series start at 0 here
		Ext.getSeriesName(seriesName)
		seriesName = replace(seriesName,"\\/","-");		//replace slashes by dashes in the seriesName
		if(seriesName != File.getName(file_list[f])) {
			saveName = File.getNameWithoutExtension(file_list[f]) + " - " + seriesName;
		}
		else saveName = File.getNameWithoutExtension(file_list[f]);

		run("CLIJ2 Macro Extensions", "cl_device=");	//Run extensions again, because only one macro extension can be loaded
		Ext.CLIJ2_clear();	

		//open the image
		if(reanalyze_boolean == false) {
			if(segmentFirstAnalyzeLater == false || alreadySegmented == true) {
	//			!!WARNING!! Fitted TCSPC .tif images exported from LAS X are saved with unit 'pixels'. Bio-Formats then ignores the pixel calibration! Using the standard ImageJ opener below is a workaround, but it doesn't work for multiseries files.
	//			if (microscope == "Confocal TCSPC / TauSeparation") run("Bio-Formats Importer", "open=["+file_list[f]+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+1);
				if (microscope == "Confocal TCSPC / TauSeparation") {
					if(substring(file_list[f], lastIndexOf(file_list[f], ".")+1) == "tif") open(file_list[f]);
					else run("Bio-Formats Importer", "open=["+file_list[f]+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+1);
				}
				else if (microscope == "TauContrast") run("Bio-Formats Importer", "open=["+file_list[f]+"] color_mode=Default rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT series_"+s+1);
				else if (microscope == "Fast FLIM") open(file_list[f]);
				else if (microscope == "Frequency Domain FLIM") openfli(file_list[f], saveName);
				else if (microscope == "Ratio Imaging") open(file_list[f]);
				else if (microscope == "Intensity only") open(file_list[f]);
			}
			setBatchMode("show");
			//Run this only once if segmentFirstAnalyzeLater == true
	
			if(segmentFirstAnalyzeLater == true && alreadySegmented == false && microscope != "Frequency Domain FLIM") {
				for(file=0; file<file_list.length; file++) {
					open(file_list[file]);
					input_image = getTitle();
					getDimensions(width, height, channels, slices, frames);
					print("Opening "+File.getName(file_list[file]));
					if(removeLastFrame_boolean == true && frames>1) {
						Stack.setFrame(frames);
						run("Delete Slice", "delete=frame");
						getDimensions(width, height, channels, slices, frames);
					}
					if(correctDrift_boolean == true && frames > 1) correct_drift(input_image);
				}
				if(file_list.length>1) run("Concatenate...", "all_open title=hyperstack open");
				getDimensions(width, height, channels, slices, frames);
				run("Stack to Hyperstack...", "order=xyctz channels="+channels+" slices="+slices*file_list.length+" frames="+frames/file_list.length+" display=Grayscale");
				if(frames>1 && slices>1) run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
				Stack.setDisplayMode("grayscale");
				rename("Stack_all");
				input_image = getTitle();
			}
			else if(segmentFirstAnalyzeLater == false || alreadySegmented == true) {
				rename(saveName);
				input_image = getTitle();
				getDimensions(width, height, channels, slices, frames);
				print("\nAnalyzing "+saveName+"\n");
	
				if(removeLastFrame_boolean == true && frames>1) {
					Stack.setFrame(frames);
					run("Delete Slice", "delete=frame");
					getDimensions(width, height, channels, slices, frames);
				}
		
				if(correctDrift_boolean == true && frames > 1) correct_drift(input_image);
			}
			//Nuclei image and pixel size (from .lif file)
			if(nucleiChannel != -1 && segmentFirstAnalyzeLater == false) {
				run("Bio-Formats Importer", "open=["+lif_file+"] autoscale color_mode=Default rois_import=[ROI manager] specify_range view=Hyperstack stack_order=XYCZT series_list="+s+1+" c_begin_"+s+1+"="+nucleiChannel+" c_end_"+s+1+"="+nucleiChannel+" c_step_"+s+1+"=1");
				nuclei_stack = getTitle();
				getPixelSize(unit, pixelWidth, pixelHeight);
			}
			else if(nucleiChannel != -1 && segmentFirstAnalyzeLater == true) exit("Sorry, segmenting all images first before analyzing is not supported in combination with a nuclei channel!");
			else if(nucleiChannel == -1) {	//Get pixel size from .tif file
				getPixelSize(unit, pixelWidth, pixelHeight);
			}

//In case the nuclei are an extra channel in the exported .tif file instead of in the .lif file: remove the comments from the following 'if statement'
		/*	if(nucleiChannel != -1) {
				run("Duplicate...", "duplicate title="+saveName+"_nuclei channels="+nucleiChannel);
				nuclei_stack = getTitle();
				selectWindow(input_image);
				Stack.setChannel(nucleiChannel);
				run("Delete Slice", "delete=channel");
			}
		*/
			selectWindow(input_image);
			getDimensions(width, height, channels, slices, frames);
			if(frames == 1) singleFrame == true;
			
			frameInterval = Stack.getFrameInterval();
			if(frameInterval != 0) print("Frame Interval detected: "+frameInterval+" s");
			else {
				print("WARNING: Frame Interval not found! Using manual value of "+default_frameInterval+" s.");
				frameInterval = default_frameInterval; 
			}

		//TO DO: get frame interval from metadata for FDFLIM:
		//  FLIMIMAGE: TIMESTAMPS - t0 = 30905379 4265506959
		//  FLIMIMAGE: TIMESTAMPS - t1 = 30905380 526102
		
			for(c=channels; c>=1; c--) {
				Stack.setChannel(c);
				run("Enhance Contrast", "saturated=0.35");
			}
			
			//Retrieve the intensity and lifetime stacks using the correct modality
			if(microscope == "Confocal TCSPC / TauSeparation") lifetime_and_intensity_stacks = calculate_lifetime_and_intensity_TCSPC(input_image);
			else if(microscope == "TauContrast") lifetime_and_intensity_stacks = calculate_lifetime_and_intensity_TauContrast(input_image);
			else if(microscope == "FAST FLIM") lifetime_and_intensity_stacks = calculate_lifetime_and_intensity_FASTFLIM(input_image);
			else if(microscope == "Frequency Domain FLIM") lifetime_and_intensity_stacks = calculate_lifetime_and_intensity_FDFLIM(input_image);
			else if(microscope == "Ratio Imaging") lifetime_and_intensity_stacks = calculate_ratio_and_intensity(input_image);
			else if(microscope == "Intensity only") lifetime_and_intensity_stacks = calculate_intensity(input_image);
			intensity_stack = lifetime_and_intensity_stacks[0];
			
			if(segmentationChannel != -1) {
				selectImage(input_image);
				run("Duplicate...", "title=["+input_image+"_for_segmentation] duplicate channels="+segmentationChannel);
				intensity_stack_for_segmentation = getTitle();
				run("Grays");
				setBatchMode("show");
			}
			else intensity_stack_for_segmentation = intensity_stack;
			
			selectWindow(intensity_stack);
			run("Grays");
			setBatchMode("show");
			lifetime_stack = lifetime_and_intensity_stacks[1];
			selectWindow(lifetime_stack);
			setBatchMode("show");

			//Correct bi-directional scanning phase mismatch
			if(correctBidirect_boolean == true) {
				intensity_stack_corrected = correct_bidirectional_phase(intensity_stack);
				close(intensity_stack);
				selectImage(intensity_stack_corrected);
				rename(intensity_stack);
				setBatchMode("show");
			}

			//Segment the cells
			if(alreadySegmented == false) {
				if(CellposeModel != "manual segmentation") {
					if(nucleiChannel != -1) labelmaps = segment_cells(intensity_stack_for_segmentation, nuclei_stack);
					else labelmaps = segment_cells_no_nuclei(intensity_stack_for_segmentation);
					labelmap_cells = labelmaps[0];
					nr_cells = getValue("Max");
					if(nucleiChannel != -1) labelmap_nuclei = labelmaps[1];
				}
				else {
					labelmaps = segment_cells_manually(intensity_stack_for_segmentation);
					labelmap_cells = labelmaps[0];
					nr_cells = getValue("Max");
					if(nr_cells > 255) resetMinAndMax();
					else setMinAndMax(0, 255);
				}
			}

			//Save all labelmaps if segmentFirstAnalyzeLater == true and break out of the loop
			if(segmentFirstAnalyzeLater == true && alreadySegmented == false) {
				File.makeDirectory(output+"labelmaps");
				selectWindow(labelmap_cells);
				rename("labelmap_cells_");
				run("Image Sequence... ", "dir=["+output+"labelmaps] format=TIFF digits=4");
				File.makeDirectory(output+"intensities");
				selectWindow("intensity_image_for_Cellpose");
				rename("intensity_");
				run("Image Sequence... ", "dir=["+output+"intensities] format=TIFF digits=4");
//				File.makeDirectory(output+"lifetimes");
//				selectWindow(lifetime_stack);
//				run("Image Sequence... ", "dir="+output+"lifetimes format=TIFF digits=4");				
				alreadySegmented = true;
				f--;
				continue;
			}
			//Open a single labelmap and intensity image from disk
			if(segmentFirstAnalyzeLater == true && alreadySegmented == true) {
				print("Analyzing "+File.getName(file_list[f]));
				open(output+"labelmaps"+File.separator+"labelmap_cells_"+IJ.pad(f,4)+".tif");
				labelmap_cells = getTitle();
				nr_cells = getValue("Max");
				setBatchMode("show");
				open(output+"intensities"+File.separator+"intensity_"+IJ.pad(f,4)+".tif");
				rename("intensity_image_for_Cellpose");
				setBatchMode("show");
				run("Duplicate...", "title=intensity");	//Necessary later to create the RGB overlay
				setBatchMode("show");
//				open(output+"lifetimes"+File.separator+"intensity_stack"+IJ.pad(f,4)+".tif");
//				lifetime_stack = getTitle();
//				setBatchMode("show");
			}
			
			if(nr_cells == 0) {
				print("No cells found in this image!");
				s++;
				break;
			}
			//Overlay labelmap with intensity image
			selectWindow("intensity_image_for_Cellpose");
			run("Add Image...", "image="+labelmap_cells+" x=0 y=0 opacity=33 zero");

			//Measure the lifetime traces
			lifetimeInfo = measure_lifetime_traces(intensity_stack, lifetime_stack, labelmap_cells, saveName, nr_cells);
			lifetimeTable = lifetimeInfo[0];	//Table name, but the table does not exist if generateTables is false
			kymograph = lifetimeInfo[1];
			kymograph_smoothed = lifetimeInfo[2];
			kymographSorted = lifetimeInfo[3];
			rank_responseVector = lifetimeInfo[4];

			if(additionalChannel > 0) {
				additionalChannelInfo = measure_intensity_in_additional_channel(input_image, additionalChannel);
				kymograph_additionalChannel = additionalChannelInfo[0];
				additionalChannelMaxIntensity = additionalChannelInfo[1];
				if(createScatterPlots == true) scatterPlotAdditional = makeScatterPlot2D(kymograph, kymograph_additionalChannel, additionalChannelMaxIntensity, saveName, nr_cells);
			}
			
			if(nucleiChannel != -1) labelmap_reassigned_nuclei_and_empty_cells = assign_nuclei_to_cells(labelmap_cells, labelmap_nuclei);

			if(microscope != "Intensity only" && createRGBOverlay == true) {
				if(overlayMovie == "intensity movie") RGB_overlay = overlay_intensity(intensity_stack, lifetime_stack, saveName, smoothRadiusOverlayXY, smoothRadiusOverlayTime);
				else if(overlayMovie == "intensity image") RGB_overlay = overlay_intensity("intensity", lifetime_stack, saveName, smoothRadiusOverlayXY, smoothRadiusOverlayTime);
			}

			if(createHistogram == true) {
				histograms = makeHistogram(kymograph, saveName);
				lifetimesHistogram = histograms[0];
				densityPlot = histograms[1];
			}
			
			if(createScatterPlots == true) {
				scatterPlotIntensity = makeScatterPlot(kymograph, "Cell_statistics", "MEAN_INTENSITY", saveName, nr_cells);
				scatterPlotArea = makeScatterPlot(kymograph, "Cell_statistics", "PIXEL_COUNT", saveName, nr_cells);
			}
		}
		else if(reanalyze_boolean == true) {
			if(seriesName != File.getName(file_list[f])) {
				saveName = File.getNameWithoutExtension(file_list[f]) + " - " + seriesName;
			}
			else saveName = File.getNameWithoutExtension(file_list[f]);
			inputDir = output;

			//Open saved tables and images
			lifetimeTable = saveName + "_lifetime.tsv";
			if(File.exists(inputDir + lifetimeTable)) open(inputDir + lifetimeTable);
			intensityTable = saveName + "_intensity.tsv";
			if(File.exists(inputDir + intensityTable)) open(inputDir + intensityTable);
			additionalChannelTable = saveName + "_intensity_ch_"+additionalChannel+".tsv";
			if(File.exists(inputDir + additionalChannelTable)) {
				open(inputDir + additionalChannelTable);
			}
			
			labelmap_cells = saveName + " (labelmap_cells).tif";
			open(inputDir + labelmap_cells);	//N.B. Currently 'labelmap_reassigned_nuclei_and_empty_cells' is not saved (yet). Fix this in the future.
			getDimensions(width, height, channels, slices, frames);
			getPixelSize(unit, pixelWidth, pixelHeight);
			nr_cells = getValue("Max");
			setBatchMode("show");
			open(inputDir + saveName + " (ROIs).zip");
			
			kymograph = saveName + " (kymograph).tif";
			open(inputDir + kymograph);
			setBatchMode("show");
			
			if(smoothRadiusTraces > 0) {
				kymograph_smoothed = saveName + " (kymograph smoothed).tif";
				open(inputDir + kymograph_smoothed);
				setBatchMode("show");
			}

			kymograph_sorted = saveName + " (kymograph sorted).tif";
			if(File.exists(inputDir + kymograph_sorted)) open(inputDir + kymograph_sorted);
			setBatchMode("show");

			rank_responseVector = saveName + " (rank vector).tif";
			if(File.exists(inputDir + rank_responseVector)) open(inputDir + rank_responseVector);
			setBatchMode("show");
			
			if(additionalChannel > 0) {
				kymograph_additionalChannel = saveName + " (kymograph additional channel).tif";
				if(File.exists(inputDir + kymograph_additionalChannel)) {
					open(inputDir + kymograph_additionalChannel);
					setBatchMode("show");
				}
			}
			
			RGB_overlay = saveName + " (lifetime & intensity RGB overlay).tif";
			if(createRGBOverlay) open(inputDir + RGB_overlay);
//			getDimensions(width, height, channels, slices, frames);
//			getPixelSize(unit, pixelWidth, pixelHeight);
			tileSizeX = (width / (nrTilesX - (nrTilesX-1)*tileOverlap) * pixelWidth)/1E6;
			tileSizeY = (height / (nrTilesY - (nrTilesY-1)*tileOverlap) * pixelHeight)/1E6;
			setBatchMode("show");

			lifetimesHistogram = saveName + " (lifetime histogram).tif";
			if(File.exists(inputDir + lifetimesHistogram)) {
				open(inputDir + lifetimesHistogram);
				setBatchMode("show");
			}
			densityPlot = saveName + " (density plot).tif";
			if(File.exists(inputDir + densityPlot)) {
				open(inputDir + densityPlot);
				setBatchMode("show");
			}
			
			scatterPlotIntensity = saveName + " (scatterplot MEAN_INTENSITY).tif";
			if(File.exists(inputDir + scatterPlotIntensity)) {
				open(inputDir + scatterPlotIntensity);
				setBatchMode("show");
			}
			
			scatterPlotArea = saveName + " (scatterplot PIXEL_COUNT).tif";
			if(File.exists(inputDir + scatterPlotArea)) {
				open(inputDir + scatterPlotArea);
				setBatchMode("show");
			}	
				
			plot = saveName + " (lifetime traces plot).tif";
			if(File.exists(inputDir + plot)) {
				open(inputDir + plot);
				setBatchMode("show");
			}
			else plot = saveName + " (lifetime histogram plot).tif";
			if(File.exists(inputDir + plot)) {
				open(inputDir + plot);
				setBatchMode("show");
			}
			Plot.getValues(timeArray, lifetimeArray);
			frameInterval = timeArray[1] - timeArray[0];
			
			string = File.openAsString(inputDir + saveName + "_Stim_&_Cal_frames.txt");
			splitString = split(string, "\t");
			stimulationFrame = parseInt(splitString[0]);
			calibrationFrame = parseInt(splitString[1]);
			if(manualStimCalFrames != "0") print("Stimulation and calibration frames (from file): "+stimulationFrame+" and "+calibrationFrame+" ("+d2s(stimulationFrame*frameInterval,1)+" and "+d2s(calibrationFrame*frameInterval,1)+" sec)");
			else {
				print("Baseline-only mode active");
				baselineOnly = true;
			}
		}

		//Find hits in a screen
		if(runScreen == true) {
			if(reanalyze_boolean == false) selectWindow(input_image);
			else selectWindow(labelmap_cells);
			metadata = split(getImageInfo(), "\n");
			nrTilesX = 1;
			nrTilesY = 1;
			if(metadata[0] == "Stitched image with tile layout:") {	//Get stage coordinates from metadata
				tiles_ = split(metadata[1], ',');
				nrTilesX = (tiles_[0]);
				nrTilesY = (tiles_[1]);
				tiles = (nrTilesX)*(nrTilesY);
				print("Stitched image detected, with tile layout: ["+nrTilesX+","+nrTilesY+"]");
				stagePositions = split(metadata[3], ",");
			}
			else if(file_list.length == 1) {	//Get stage coordinates for single (stitched) image from .lif file
				//Get stage coordinates of all the tiles. N.B. Counting always starts at the first tile in the .lif file.
				oldLogWindow = getInfo("log");
				print("\\Clear");
				tiles = 1;
				run("NKI get stage coordinates to log window", "file=["+lif_file+"], nrtiles="+tiles);	//Run a Jython script to print the stage coordinates to the log window
				logWindow = getInfo("log");
				print("\\Clear");
				print(oldLogWindow);
				if(logWindow == "") {
					print("[WARNING] Stage positions not found! Did you enter a correct path to a .lif file? Using coordinates (0,0)");
					logWindow = "0\n0";
				}
				stagePositions = split(logWindow, "\n");
				print("Found stage coordinates for "+nrTilesX*nrTilesY+" tiles");
				showStatus("Found stage coordinates for "+nrTilesX*nrTilesY+" tiles");
			}	
			else {
				layout = find_tile_layout(stagePositions);
				nrTilesX = layout[0];
				nrTilesY = layout[1];
			}
			
			//Get tile size in meters; N.B. allows only square tiles!
			if(metadata[4] == "Tile size (um):") {
				tileSizeArray = split(metadata[5], ",");
				tileSizeX = (tileSizeArray[0])/1E6;
				tileSizeY = (tileSizeArray[1])/1E6;
			}
			else {
				tileSizeX = (width / (nrTilesX - (nrTilesX-1)*tileOverlap) * pixelWidth)/1E6;
				tileSizeY = (height / (nrTilesY - (nrTilesY-1)*tileOverlap) * pixelHeight)/1E6;
			}
			if(metadata[6] == "overlap (%):") {
				tileOverlap = (metadata[7]);
				print("Found overlap in metadata: "+tileOverlap+"%");
			}
			if(debugMode) {
				print("tileSizeX = "+tileSizeX);
				print("tileSizeY = "+tileSizeY);
				print("nrTilesX = "+nrTilesX);
				print("nrTilesX = "+nrTilesY);
			}
			if(additionalChannel == -1) kymograph_additionalChannel = "";
			if(tileOverlap >= 0) {
				if(nucleiChannel != -1) hitList = find_hits(lifetimeTable, labelmap_reassigned_nuclei_and_empty_cells, kymograph, f);
				else hitList = find_hits(labelmap_cells, kymograph, kymograph_additionalChannel, f);
			}
			else if(tileOverlap < 0 && file_list.length == 1) {
				if(nucleiChannel != -1) hitList = find_hits(lifetimeTable, labelmap_reassigned_nuclei_and_empty_cells, kymograph, nrTilesX-1);	//Send the upper leftmost tile. Scanning with overlap < 0 goes from right to left.
				else hitList = find_hits(labelmap_cells, kymograph, kymograph_additionalChannel, nrTilesX-1);		
			}

			if(createRGBOverlay == true && singleFrame == false) color_hits_on_RGB_overlay(RGB_overlay, hitList, "white", f);
			else if(singleFrame == true) color_hits_on_RGB_overlay(input_image, hitList, "white", f);
			if(singleFrame == false && smoothRadiusTraces == 0) plot_hit_timetraces(kymograph, hitList, f, saveName, nr_cells);
			else if(singleFrame == false && smoothRadiusTraces > 0) plot_hit_timetraces(kymograph_smoothed, hitList, f, saveName, nr_cells);
			if(singleFrame == false && createNormalizedKymograph == true) plot_normalized_timetraces(kymograph+"_baseline_normalized", saveName, nr_cells);		
			if(createHistogram == true && nrHits>0) {
				hit_cell_IDs = Table.getColumn("Cell", hitList);
				selectWindow(kymograph);
				getDimensions(kwidth, kheight, kchannels, kslices, kframes);
				kymographID = getImageID();
				//Create hit-only kymograph
				kymograph_hits = "Kymograph of Hits";
				newImage(kymograph_hits, "32-bit black", hit_cell_IDs.length, kheight, 1);
				kymograph_nonhits = "Kymograph of non-hits";
				newImage(kymograph_nonhits, "32-bit black", nr_cells - hit_cell_IDs.length, kheight, 1);
				hitnr = 0;
				nonhitnr = 0;
				
				for(i=0; i<hit_cell_IDs.length; i++) {
					selectImage(kymographID);
					makeRectangle(hit_cell_IDs[i]-1, 0, 1, kheight);
					run("Copy");
					selectImage(kymograph_hits);
					makeRectangle(i, 0, 1, kheight);
					run("Paste");
				}

				for(i=0; i<nr_cells; i++) {
					selectImage(kymographID);
					makeRectangle(i, 0, 1, kheight);
					if(!occursInArray(hit_cell_IDs, i+1)) { //a hit
						run("Copy");
						selectImage(kymograph_nonhits);		//not a hit
						makeRectangle(nonhitnr, 0, 1, kheight);
						run("Paste");
						nonhitnr++;
					}
				}
				selectWindow(kymograph_hits);
				run(lut);
				setMinAndMax(minLifetime, maxLifetime);
				setBatchMode("show");
				histograms = makeHistogram(kymograph_hits, saveName+"_hits");
				lifetimesHistogramHits = histograms[0];
				densityPlotHits = histograms[1];

				selectWindow(kymograph_nonhits);
				run(lut);
				setMinAndMax(minLifetime, maxLifetime);
				setBatchMode("show");
				histograms = makeHistogram(kymograph_nonhits, saveName+"_nonhits");
				lifetimesHistogramNonHits = histograms[0];
				densityPlotNonHits = histograms[1];
			}
		}

		//Save images and data tables
		run("Set Measurements...", "mean redirect=None decimal=9");	//Make sure that enough decimals are saved!
		if(reanalyze_boolean == false) {	//Don't save these when re-analyzing
			selectWindow(lifetime_stack);
			saveAs("tiff", output + saveName + " (weighted lifetime)");

			selectWindow(labelmap_cells);
			if(runScreen == true) setMetadata("Info", metadata[0]+"\n"+metadata[1]+"\n"+metadata[2]+"\n"+metadata[3]+"\n"+metadata[4]+"\n"+metadata[5]+"\n"+metadata[6]+"\n"+metadata[7]+"\n");	//Write the tile layout and stage positions to metadata of the labelmap
			Stack.setXUnit("um");
			run("Properties...", "pixel_width="+pixelWidth+" pixel_height="+pixelHeight+" voxel_depth=1.0000 frame="+frameInterval);
			saveAs("tiff", output + saveName + " (labelmap_cells)");

			selectWindow("intensity_image_for_Cellpose");
			Stack.setXUnit("um");
			run("Properties...", "pixel_width="+pixelWidth+" pixel_height="+pixelHeight+" voxel_depth=1.0000 frame="+frameInterval);
			saveAs("tiff", output + saveName + " (intensity & labelmap)");

			if(frames>1) {
				Ext.CLIJ2_pull(rank_responseVector);
				run(lut);
				setBatchMode("show");
				saveAs("tiff", output + saveName + " (rank vector)");
				selectWindow(kymograph);
//				setMinAndMax(minLifetime, maxLifetime);
				saveAs("tiff", output + saveName + " (kymograph)");
				selectWindow(kymographSorted);
//				setMinAndMax(minLifetime, maxLifetime);
				saveAs("tiff", output + saveName + " (kymograph sorted)");
				if(smoothRadiusTraces > 0) {
					selectWindow(kymograph_smoothed);
//					setMinAndMax(minLifetime, maxLifetime);
					saveAs("tiff", output + saveName + " (kymograph smoothed)");
				}
				if(additionalChannel > 0) {
					selectWindow(kymograph_additionalChannel);
//					setMinAndMax(minLifetime, maxLifetime);
					saveAs("tiff", output + saveName + " (kymograph additional channel)");
				}
			}

			if(createRGBOverlay == true && microscope != "Intensity only") {
				selectWindow(RGB_overlay);
				if(frames>1) Stack.setFrame(1);
				updateDisplay();
				saveAs("tiff", output + saveName + " (lifetime & intensity RGB overlay)");
			}
			if(createHistogram == true) {	
				selectWindow(lifetimesHistogram);
				saveAs("tiff", output + lifetimesHistogram);
				selectWindow(densityPlot);
				saveAs("tiff", output + densityPlot);
				if(runScreen == true && nrHits > 0) {
					selectWindow(lifetimesHistogramHits);
					saveAs("tiff", output + saveName + lifetimesHistogramHits);
					selectWindow(lifetimesHistogramNonHits);
					saveAs("tiff", output + saveName + lifetimesHistogramNonHits);
					selectWindow(densityPlotHits);
					saveAs("tiff", output + saveName + densityPlotHits);
					selectWindow(densityPlotNonHits);
					saveAs("tiff", output + saveName + densityPlotNonHits);
				}
			}
			if(createScatterPlots == true) {	
				selectWindow(scatterPlotIntensity);
				saveAs("tiff", output + scatterPlotIntensity);
				selectWindow(scatterPlotArea);
				saveAs("tiff", output + scatterPlotArea);
				if(additionalChannel > 0) {
					selectWindow(scatterPlotAdditional);
					saveAs("tiff", output + scatterPlotAdditional);
				}
			}
			roiManager("deselect");
			roiManager("Remove Frame Info");
			roiManager("save", output + saveName + " (ROIs).zip");

			//Save the log window as text file
			txtfile = File.open(output + saveName + "_log.txt");
			print(txtfile, getInfo("log"));
			File.close(txtfile);

			if(runScreen == true) {
				stagePositionsFile = File.open(output + "Stage_positions.txt");
				print(stagePositionsFile, arrayToString(stagePositions, "\n"));
				File.close(stagePositionsFile);
	
				selectWindow("Cell_coordinates");
				Table.save(output + "Cell_coordinates.tsv");
			}
		}
		Ext.CLIJ2_clear();
	}
}

//Create .rgn files in chucks of 'maxNrHits' events.
if(runScreen == true) {
//	outputName = File.getNameWithoutExtension(lif_file);
	outputName = saveName;
	selectWindow(hitList);
//	if(sort_hits_column != "do not sort" && Table.size > 0) Table.sort("measurement");	//To do: create reverseTable(tableName) function to be able to sort top-down
	Table.save(output + outputName + " (Hit list).tsv");
	Table.update;

	generate_rgn_file(outputName, hitList);
	if(Table.size > 0) {
		nrHits = Table.size;
		n = 0;
		while(n < nrHits) {
			showStatus("Creating Hit tables...");
			showProgress(n, nrHits);
			open(output + outputName + " (Hit list).tsv");	//re-open file to generate the top N hits table (couldn't find a command to duplicate the table)
			currentHitListName = outputName + " (Hit list "+n+"-"+n+minOf(maxNrHits, nrHits)-1+").tsv";
			Table.deleteRows(n+maxNrHits, nrHits);
			Table.update;
			Table.deleteRows(0, n-1);
			Table.update;
			Table.rename(outputName + " (Hit list).tsv", currentHitListName);
			Table.save(output + currentHitListName);
			generate_rgn_file(outputName + "_" +n+"-"+n+minOf(maxNrHits, nrHits)-1, currentHitListName);
			close(currentHitListName);
			n += maxNrHits;
		}
	}
}

//Create top-N hits .rgn file
if(runScreen == true) {
//	outputName = File.getNameWithoutExtension(lif_file);
	selectWindow(hitList);
	if(Table.size > 0) {
		open(output + outputName + " (Hit list).tsv");	//re-open file to generate the top N hits table (couldn't find a command to duplicate the table)
		topHitList = outputName + " (Top "+topNHits+" Hit list).tsv";
		Table.deleteRows(topNHits, Table.size);
		Table.update;
		Table.rename(outputName + " (Hit list).tsv", topHitList);
		Table.save(output + outputName + " (Top "+topNHits+" Hit list).tsv");
		generate_rgn_file(outputName + "_TOP_"+topNHits+"_", topHitList);
		plot_hit_locations(hitList, topHitList, tiles);
		saveAs("tiff", output + saveName + " (hit positions)");
	}
}

restoreSettings();
print("Done without errors");




/////////////////////// END /////////////////////////





function makeHistogram(kymographImage, saveName) {
	selectWindow(kymographImage);
	getDimensions(width, height, channels, slices, frames);
	nrTimePoints = height;
	plotString = "";
	maxBinHeight = 0;
	
	densityPlot = saveName + " (density plot)";
	newImage(densityPlot, "16-bit black", nrTimePoints, histogramBins, 1);
//	Plot.create(densityPlot, "Lifetime (ns)", "time");
//	Plot.show();
//	timeArray = newArray(nrTimePoints);
	
	selectWindow(kymographImage);
	kymographImageID = getImageID();
	getDimensions(width, height, channels, slices, frames);
	for(t=0; t<nrTimePoints; t++) {
		showStatus("Creating histogram plot...");
		showProgress(t, nrTimePoints);
		Plot.create(t, "Lifetime (ns)", "counts");
		Plot.setFrameSize(500, 500);
		Plot.setAxisLabelSize(axisFontSize);
		Plot.setFontSize(axisFontSize);
		selectImage(kymographImageID);
		makeRectangle(0, t, width, 1);
		getStatistics(area, mean, min, max, std, histogram);
		getHistogram(values, counts, histogramBins, minLifetime, maxLifetime);
		Plot.add("bar", values, counts);
		Plot.setStyle(0, "#0000a0, #00a0ff, 1.0, Separated Bars");
		Plot.setJustification("right");
		Plot.setFontSize(axisFontSize);
		Plot.addText(d2s(t*frameInterval,0)+" sec", 0.17, 0.08);
		Plot.addText("Mean ± sd: "+d2s(mean,2)+" ± "+d2s(std,2), 0.97, 0.08);
		Array.getStatistics(counts, minCounts, maxCounts);
		maxBinHeight = maxOf(maxBinHeight, maxCounts);
//		if(displayGrid == true) Plot.setFormatFlags("11000000111111");
//		else Plot.setFormatFlags("11000000001111");
		Plot.setFormatFlags("11000000001111");
		plotString += " image"+t+1+"=["+t+"]";
		Plot.show();
		
//		timeArray = Array.fill(timeArray, t);
//		selectWindow(densityPlot);
//		Plot.add("box", timeArray, values);
		selectWindow(densityPlot);
		for (y = 0; y < histogramBins; y++) setPixel(t, y, counts[y]);
	}
//	run("Calibration Bar...", "location=[Upper Left] fill=None label=White number=4 decimal=0 font=14 zoom=1 bold overlay");
	for(t=0; t<nrTimePoints; t++) {
		selectWindow(t);
		Plot.setLimits(minLifetime, maxLifetime, 0, maxBinHeight);
	}
	plotName = saveName + " (lifetime histogram)";
	if(nrTimePoints>1) run("Concatenate...", "  title=["+plotName+"] open "+plotString);
	else rename(plotName);
	setBatchMode("show");
	
	selectWindow(densityPlot);
	setBatchMode("show");
	run("Flip Vertically");
	resetMinAndMax();
	run("Turbo");
	run("In [+]");
	run("In [+]");
	return newArray(plotName, densityPlot);
}


//Plot lifetime vs additional channel intensity
function makeScatterPlot2D(kymograph, kymograph_additionalChannel, maxIntensity, saveName, nr_cells) {
	selectWindow(labelmap_cells);	//Any image with the 'glasbey on dark' LUT
	getLut(reds, greens, blues);
	
	plotString = "";
	
	lifetimes = newArray(nr_cells);
	intensities = newArray(nr_cells);
	selectWindow(kymograph);
	getDimensions(kwidth, kheight, kchannels, kslices, kframes);
	nrTimePoints = kheight;
	for(t=0; t<nrTimePoints; t++) {
		showStatus("Creating additional channel scatter plot...");
		showProgress(t, nrTimePoints);
		Plot.create(t, "Intensity_ch"+additionalChannel, "Lifetime (ns)");
		Plot.setFrameSize(500, 500);
		Plot.setAxisLabelSize(axisFontSize);
		Plot.setFontSize(axisFontSize);
		for(i=0; i<nr_cells; i++) {
			selectImage(kymograph);
			lifetimes[i] = getPixel(i, t);
			selectImage(kymograph_additionalChannel);
			intensities[i] = getPixel(i, t);
			color = getLabelColor(i, nr_cells);
			Plot.setColor(color);
			x = newArray(1);
			x[0] = intensities[i];
			maxIntensity = maxOf(maxIntensity, intensities[i]);
			y = newArray(1);
			y[0] = lifetimes[i];
			Plot.add("box", x, y);
			Plot.setStyle(i, ""+color+","+color+", 1.0, Box");
		}
		Plot.setLimits(0, maxIntensity, minLifetime, maxLifetime);
		Plot.addText(d2s(t*frameInterval,0)+" sec", 0.17, 0.08);
		if(displayGrid == true) Plot.setFormatFlags("11000000111111");
		else Plot.setFormatFlags("11000000001111");
		plotString += " image"+t+1+"=["+t+"]";
		Plot.show();
	}
	plotName = saveName + " (scatterplot ch"+additionalChannel+")";
	if(nrTimePoints>1) run("Concatenate...", "  title=["+plotName+"] open "+plotString);
	else rename(plotName);
	setBatchMode("show");
	return plotName;
}


function makeScatterPlot(kymograph, cellStatsTable, parameter, saveName, nr_cells) {
	selectWindow(labelmap_cells);	//Any image with the 'glasbey on dark' LUT
	getLut(reds, greens, blues);

	selectWindow(cellStatsTable);
	xValues = Table.getColumn(parameter);
	Array.getStatistics(xValues, minValue, maxValue, meanValue, stdDev);
	
	plotString = "";

	lifetimes = newArray(nr_cells);
	selectWindow(kymograph);
	getDimensions(kwidth, kheight, kchannels, kslices, kframes);
	nrTimePoints = kheight;
	for(t=0; t<nrTimePoints; t++) {
		showStatus("Creating "+parameter+" scatter plot...");
		showProgress(t, nrTimePoints);
		Plot.create(t, parameter, "Lifetime (ns)");
		Plot.setFrameSize(500, 500);
		Plot.setAxisLabelSize(axisFontSize);
		Plot.setFontSize(axisFontSize);
		for(i=0; i<nr_cells; i++) {
			lifetimes[i] = getPixel(i, t);
			color = getLabelColor(i, nr_cells);
			Plot.setColor(color);
			x = newArray(1);
			x[0] = xValues[i];
			y = newArray(1);
			y[0] = lifetimes[i];
			Plot.add("box", x, y);
			Plot.setStyle(i, ""+color+","+color+", 1.0, Box");
		}
		Plot.setLimits(0, maxValue, minLifetime, maxLifetime);
		Plot.addText(d2s(t*frameInterval,0)+" sec", 0.17, 0.08);
		if(displayGrid == true) Plot.setFormatFlags("11000000111111");
		else Plot.setFormatFlags("11000000001111");
		plotString += " image"+t+1+"=["+t+"]";
		Plot.show();
	}
	plotName = saveName + " (scatterplot "+parameter+")";
	if(nrTimePoints>1) run("Concatenate...", "  title=["+plotName+"] open "+plotString);
	else rename(plotName);
	setBatchMode("show");
	return plotName;
}


function plot_normalized_timetraces(kymograph_baseline_normalized, saveName, nr_cells) {
	selectWindow(kymograph_baseline_normalized);
	getDimensions(kwidth, kheight, kchannels, kslices, kframes);
	timeArray = Array.getSequence(kheight);
	timeArray = multiplyArraywithScalar(timeArray, frameInterval);

	plotName = saveName + " (lifetime traces plot NORMALIZED)";
	Plot.create(plotName, "time (s)", "Lifetime (ns)");
	Plot.setFrameSize(900, 600);
	Plot.setAxisLabelSize(axisFontSize);
	Plot.setFontSize(axisFontSize);
	if(microscope != "Intensity only") Plot.setLimits(0, kheight*frameInterval, minLifetime, maxLifetime);
	
	selectWindow(labelmap_cells);	//Any image with the 'glasbey on dark' LUT
	getLut(reds, greens, blues);
	
	selectWindow(kymograph_baseline_normalized);
	for(i=0;i<nr_cells;i++) {
		makeRectangle(i, 0, 1, kheight);
		setKeyDown("alt");	//ensure vertical profile
		lifetimeData = getProfile();
		setKeyDown("none");
		//Generate lines with 'glasbey on dark' colors
		color = getLabelColor(i, nr_cells);
		Plot.setColor(color);
		Plot.add("line", timeArray, lifetimeData);
	}
	Plot.addFromPlot(saveName + " (lifetime traces plot).tif", nr_cells);
	Plot.setStyle(nr_cells, "black,none,3.0,Line");
	Plot.show();
	setBatchMode("show");
	saveAs("tiff", output + saveName + " (lifetime traces plot baseline-normalized)");
}


function plot_hit_timetraces(kymograph, hitList, tile, saveName, nr_cells) {
	selectWindow(hitList);
	nrHits = Table.size;
	selectWindow(kymograph);
	getDimensions(kwidth, kheight, kchannels, kslices, kframes);
	timeArray = Array.getSequence(kheight);
	timeArray = multiplyArraywithScalar(timeArray, frameInterval);

	//Add detection window to the original plot - Doesn't work, because ImageJ thinks that it is not a plot.
//	selectWindow(saveName + " (lifetime traces plot).tif");
//	Plot.setColor("#003399");
//	Plot.getLimits(xMin, xMax, yMin, yMax);
//	Plot.setLineWidth(3);
//	Plot.drawShapes("rectangles", responseTimeStart*frameInterval, yMax, (responseTimeStart + hit_avg_response_time_window - 1)*frameInterval, yMin);
	
	plotName = saveName + " (lifetime traces plot HITS only)";
	Plot.create(plotName, "time (s)", "Lifetime (ns)");
	if(microscope != "Intensity only") Plot.setLimits(0, kheight*frameInterval, minLifetime, maxLifetime);
	Plot.setFrameSize(900, 600);
	Plot.setAxisLabelSize(axisFontSize);
	Plot.setFontSize(axisFontSize);
		
	selectWindow(labelmap_cells);	//Any image with the 'glasbey on dark' LUT
	getLut(reds, greens, blues);
	
	if(nrHits>0) hit_cell_nr = Table.getColumn("Cell", hitList);	//Get the cell number of all the hits
	selectWindow(kymograph);
	getDimensions(kwidth, kheight, kchannels, kslices, kframes);
	for(i=0;i<nrHits;i++) {
		current_tile = Table.get("Tile", i, hitList);
		if(current_tile == tile) {
			//lifetimeData = Table.getColumn("cell_"+IJ.pad(hit_cell_nr[i],5), lifetimeTable);
			makeRectangle(hit_cell_nr[i]-1, 0, 1, kheight);
			setKeyDown("alt");	//ensure vertical profile
			lifetimeData = getProfile();
			setKeyDown("none");
			//Generate traces with 'glasbey on dark' colors
			color = getLabelColor(hit_cell_nr[i]-1, nr_cells);
			Plot.setColor(color);
			Plot.add("line", timeArray, lifetimeData);
			roiManager("select", hit_cell_nr[i]-1);
			if(useColoredHitOverlay) roiManager("Set Color", color);	//Works, but looks very confusing
			else roiManager("Set Color", "white");
			roiManager("Set Line Width", 3);
		}
	}

	if(displayGrid == true) Plot.setFormatFlags("11000000111111");
	else Plot.setFormatFlags("11000000001111");
	Plot.setColor("#003399");
	Plot.setFontSize(axisFontSize*0.8);
	Plot.setJustification("left");
	Plot.addText(nrHits+" traces", 0.02, 0.05);
	Plot.addText("smoothing: "+smoothRadiusTraces+ " frames", 0.02, 0.08);
	Plot.setFontSize(axisFontSize);

	//Draw lines at the hit finding windows
//	selectWindow(saveName + " (lifetime traces plot).tif");
	Plot.setColor("#003399");
	Plot.getLimits(xMin, xMax, yMin, yMax);
	Plot.setLineWidth(2);

	if(hit_baseline_boolean == true) {
		Plot.drawLine(0, -1/0, 0, 1/0);
		Plot.drawLine(baselineLength*frameInterval, -1/0, baselineLength*frameInterval, 1/0);
	}
	if(hit_avg_response_boolean || hit_max_abs_response_boolean || hit_max_response_boolean || hit_rrr_boolean) {
		Plot.drawLine(responseTimeStart*frameInterval, -1/0, responseTimeStart*frameInterval, 1/0);
		Plot.drawLine((responseTimeStart + hit_avg_response_time_window - 1)*frameInterval, -1/0, (responseTimeStart + hit_avg_response_time_window - 1)*frameInterval, 1/0);
	}
	Plot.setLineWidth(1);
	Plot.show();
	
	//Add semi-transparent rectangle to the hit finding windows
	Plot.getFrameBounds(plotX, plotY, plotWidth, plotHeight);
	plotPixelWidth = kheight*frameInterval/plotWidth;
	if(hit_baseline_boolean == true) {
		makeRectangle(plotX, plotY, baselineLength*frameInterval/plotPixelWidth, plotHeight);
		Overlay.addSelection("", 0, "#22003399");
		run("Select None");
	}
	if(hit_avg_response_boolean || hit_max_abs_response_boolean || hit_max_response_boolean || hit_rrr_boolean) {
		//Add semi-transparent rectangle to the respinse hit finding window
		makeRectangle(plotX + responseTimeStart*frameInterval/plotPixelWidth, plotY, (hit_avg_response_time_window - 1)*frameInterval/plotPixelWidth, plotHeight);
		Overlay.addSelection("", 0, "#33003399");
		run("Select None");
	}
	setBatchMode("show");
	saveAs("tiff", output + saveName + " (lifetime traces plot HITS only)");
}


function getLabelColor(label, nrOfLabels) {
	if(nrOfLabels >= 255) {
		color1 = IJ.pad(toHex(reds[label/nrOfLabels*255]),2);
		color2 = IJ.pad(toHex(greens[label/nrOfLabels*255]),2);
		color3 = IJ.pad(toHex(blues[label/nrOfLabels*255]),2);
	}
	else {
		color1 = IJ.pad(toHex(reds[label]),2);
		color2 = IJ.pad(toHex(greens[label]),2);
		color3 = IJ.pad(toHex(blues[label]),2);		
	}
	labelColor = "#"+color1+color2+color3;
	return labelColor;
}



function color_hits_on_RGB_overlay(RGB_overlay, hitList, color, tile) {
	selectWindow(hitList);
	nrHits = Table.size;
//	hit_cell_numbers = Table.getColumn("Cell");
//	hit_cell_numbers_corrected = addScalarToArray(hit_cell_numbers, -1);
	selectWindow(RGB_overlay);
	roiManager("Set Color", "#555555");	//Set all selections to gray
	roiManager("Set Line Width", 1); //Set all linewitdths to 1
	roiManager("Deselect");
	RoiManager.setGroup(1);
	for (i=0; i<nrHits; i++) {
		current_tile = Table.get("Tile", i);
		if(current_tile == tile) {
			cell = Table.get("Cell", i, hitList);
			roiManager("select", cell-1);
			roiManager("Set Color", color);
			roiManager("Set Line Width", 2);
			RoiManager.setGroup(2);
		}
	}
	Stack.setFrameInterval(frameInterval);
	roiManager("show all without labels");
}


function plot_hit_locations(hitList, topHitList, tiles) {
	Plot.create("Hit positions", "X (mm)", "Y (mm)");
	Plot.setFrameSize(768, 768);

	selectWindow(hitList);
	Plot.setLineWidth(1);
	Plot.setColor("red");
	X = Table.getColumn("AbsPosX (m)");
	Y = Table.getColumn("AbsPosY (m)");
	X = multiplyArraywithScalar(X, 1000);
	Y = multiplyArraywithScalar(Y, -1000);
	Plot.add("circle", X, Y);
	
	selectWindow(topHitList);
	Plot.setLineWidth(1);
	Plot.setColor("blue");
	X = Table.getColumn("AbsPosX (m)");
	Y = Table.getColumn("AbsPosY (m)");
	X = multiplyArraywithScalar(X, 1000);
	Y = multiplyArraywithScalar(Y, -1000);
	Plot.add("circle", X, Y);
	
	//Draw tiles as gray squares
	lefts = newArray(tiles);
	tops = newArray(tiles);
	rights = newArray(tiles);
	bottoms = newArray(tiles);
	for (i = 0; i < tiles; i++) {
		lefts[i] = (d2s((-parseFloat(stagePositions[2*i+1]) - (tileSizeX)/2)*1000, 10));	//extra brackets, because ImageJ treats these as numbers
		tops[i] = (d2s((-parseFloat(stagePositions[2*i]) - (tileSizeY)/2)*-1000, 10));
		rights[i] = (d2s((-parseFloat(stagePositions[2*i+1]) - (tileSizeX)/2 + tileSizeX)*1000, 10));
		bottoms[i] = (d2s((-parseFloat(stagePositions[2*i]) - (tileSizeY)/2 + tileSizeX)*-1000, 10));
	}
	Plot.setColor("gray");
	Plot.setLineWidth(2);
	Plot.drawShapes("rectangles", lefts, tops, rights, bottoms);
//	if(displayGrid == true) Plot.setFormatFlags("11000000111111");
//	else Plot.setFormatFlags("11000000001111");
	Plot.setFormatFlags("11000000001111");
	Plot.setAxisLabelSize(axisFontSize);
	Plot.setLegend("all hits\ttop "+topNHits+" hits", "top-right");
	Plot.setFontSize(axisFontSize);	
	Plot.show();
	Plot.getLimits(xMin, xMax, yMin, yMax);
	xSpan = xMax - xMin;
	ySpan = yMax - yMin;
	span = maxOf(xSpan, ySpan);
	Plot.setLimits(xMin - tileSizeX*100, xMin + span + tileSizeX*100, yMin - tileSizeY*100, yMin + span + tileSizeY*100);
	setBatchMode("show");
}


function generate_rgn_file(name, hitList) {
	preamble = "<StageOverviewRegions>\n<Regions>\n<ShapeList>\n<Items>\n<Item0>\n<Name>hits</Name>\n<Identifier>"+generate_UUID()+"</Identifier>\n<Type>CompoundShape</Type>\n<Font />\n<Verticies>\n<Items />\n</Verticies>\n<DecoratorColors>\n<Items />\n</DecoratorColors>\n<ExtendedProperties>\n<Items />\n</ExtendedProperties>\n<Children>\n<Items>\n";
	postamble = "</Items>\n</Children>\n</Item0>\n</Items>\n<FillMaskMode>None</FillMaskMode>\n<VertexUnitMode>Pixels</VertexUnitMode>\n</ShapeList>\n</Regions>\n<StackList>\n";
	stacklist = "";
	if(File.exists(output + "HITS_"+name+".rgn")) {
		temp = File.delete(output + "HITS_"+name+".rgn");
	}
	file = File.open(output + "HITS_"+name+".rgn");
	print(file, preamble);
	selectWindow(hitList);
	nrHitsTemp = Table.size;
	for (i = 0; i < nrHitsTemp; i++) {
		UUID = generate_UUID();
		start = "<Item"+i+">\n<Number>"+i+1+"</Number>\n<Tag>Cell_"+i+1+"</Tag>\n<Identifier>"+UUID+"</Identifier>\n<Type>Point</Type>\n<Fill>R:1,G:0,B:0,A:0</Fill>\n<Font />\n<Verticies>\n<Items>\n<Item0>\n";
		middle = "\<X\>" + d2s(Table.get("AbsPosX (m)", i),10) + "\</X\>\n\<Y\>"+ d2s(Table.get("AbsPosY (m)", i),10) + "\</Y\>\n";
		end = "</Item0>\n</Items>\n</Verticies>\n<DecoratorColors>\n<Items />\n</DecoratorColors>\n<ExtendedProperties>\n<Items />\n</ExtendedProperties>\n</Item"+i+">\n";
		print(file, start+middle+end);
		stacklist += "<Entry Identifier=\""+UUID+"\" Begin=\"0.0000000000\" End=\"0.0000000000\" SectionCount=\"0\" ReferenceX=\"0.0000000000\" ReferenceY=\"0.0000000000\" FocusStabilizerOffset=\"0.0000000000\" FocusStabilizerOffsetFixed=\"false\" StackValid=\"false\" Marked=\"false\" />\n";
	}
	print(file, postamble);
	print(file, stacklist);
	print(file, "</StackList>\n</StageOverviewRegions>\n");
	File.close(file);
}


function generate_UUID() {
	string = "";
	for (i = 0; i < 36; i++) {
		if(i==8 || i==13 || i==18 || i==23) string += "-";
		else string += toHex(random*16);
	}
	return string;
}


//Find Tile layout from stage positions. Returns an array (tilesX, tilesY).
function find_tile_layout(stagePositions) {
	tilesX = 1;
	tilesY = 1;
	yPos = newArray(stagePositions.length/2);
	xPos = newArray(stagePositions.length/2);

	for(i=0; i<xPos.length; i++) {
		if(i>0) { if(!occursInArray(yPos, stagePositions[2*i])) tilesY++; }
		if(i>0) { if(!occursInArray(xPos, stagePositions[2*i+1])) tilesX++; }
		yPos[i] = stagePositions[2*i];
		xPos[i] = stagePositions[2*i+1];
	}
	if(tilesX > 1 || tilesY >1) print("Found tile layout: "+tilesX+" x "+tilesY);
	else print("No tile layout found");
	return newArray(tilesX, tilesY);
}


function find_hits(labelmap, kymograph, kymograph_additionalChannel, tile) {
	hitList = "Hit List";
	if(!isOpen(hitList)) {
		Table.create(hitList);
		Table.showRowNumbers(true);
		Table.setLocationAndSize(50, 50, 672, 1000);
		nrHits = 0;
	}
	else {
		selectWindow(hitList);
		nrHits = Table.size;
	}
	run("Set Measurements...", "mean redirect=None decimal=9");

	selectWindow(kymograph);
	getDimensions(kwidth, kheight, kchannels, kslices, kframes);
	nrTimePoints = kheight;
	if(nrTimePoints == 1) singleFrame = true;	//Then it is not a timelapse

	if(calibrationFrame == -1 && singleFrame == false) calibrationFrame = kheight - 1;	//Set calibration frame to the last frame
	else if (calibrationFrame == -1 && singleFrame == true) {
		stimulationFrame = 0;
		calibrationFrame = 1; 
	}

	if(singleFrame == false) {
		//Get average baseline lifetime
		if(calibrationOnly == false) baselineLength = stimulationFrame;
		else if(calibrationOnly == true) baselineLength = calibrationFrame - 2;	//Leave 1 frame for response
		if(baselineLength == 0) baselineLength = 1;
		selectWindow(kymograph);
		makeRectangle(0, 0, kwidth, baselineLength);
		run("Plots...", "vertical");	//Make sure the profile is vertical
		baselineAverageTrace = getProfile();
		Array.getStatistics(baselineAverageTrace, min, max, baselineAverage, baselineStdDev);
		baselineAverage = averageArray(baselineAverageTrace);	//Get average without NaNs
		//print("N.B. For calculating the averages a 1-frame margin from stimulation and calibration is always used!");
		print("Average baseline lifetime, based on frames 1 - "+baselineLength+": "+ baselineAverage + " ± " + baselineStdDev + " ns");
		
		if(baselineOnly == false && calibrationOnly == false) {
			//Get mean stimulated lifetime
			selectWindow(kymograph);
			makeRectangle(0, stimulationFrame + 1, kwidth, calibrationFrame - stimulationFrame - 1);	//N.B. For the average: start stimulation one frame later and skip the last frame
			run("Plots...", "vertical");	//Make sure the profile is vertical
			stimulationMeanTrace = getProfile();
			Array.getStatistics(stimulationMeanTrace, min, max, stimulationMean, stimulationStdDev);
			stimulationMean = averageArray(stimulationMeanTrace);	//Get average without NaNs
			print("Average response lifetime, based on frames "+stimulationFrame + 2 +" - "+calibrationFrame+": " + stimulationMean + " ± " + stimulationStdDev + " ns");	//N.B. This is with a 1 frame margin!

			//Get mean calibration lifetime
			selectWindow(kymograph);
			makeRectangle(0, calibrationFrame + 1, kwidth, kheight - calibrationFrame);
			run("Plots...", "vertical");	//Make sure the profile is vertical
			calibrationMeanTrace = getProfile();
			Array.getStatistics(calibrationMeanTrace, min, max, calibrationMean, calibrationStdDev);
			calibrationMean = averageArray(calibrationMeanTrace);	//Get average without NaNs
			print("Average calibration lifetime, based on frames "+calibrationFrame + 1 +" - "+nrTimePoints
			+": " + calibrationMean + " ± " + calibrationStdDev + " ns");
		}
		if(baselineOnly == false && calibrationOnly == true) {
			//Get mean stimulated lifetime
			selectWindow(kymograph);
			makeRectangle(0, calibrationFrame - 2, width, 1);
//print(0, calibrationFrame - 2, width, 1);
			run("Plots...", "vertical");	//Make sure the profile is vertical
			stimulationMeanTrace = getProfile();
			Array.getStatistics(stimulationMeanTrace, min, max, stimulationMean, stimulationStdDev);
			stimulationMean = averageArray(stimulationMeanTrace);	//Get average without NaNs
			print("Average response lifetime, based on frames "+calibrationFrame - 1 +" - "+calibrationFrame - 1+": " + stimulationMean + " ± " + stimulationStdDev + " ns");

			//Get mean calibration lifetime
			selectWindow(kymograph);
			makeRectangle(0, calibrationFrame - 1, width, nrTimePoints - calibrationFrame + 1);
//print(0, calibrationFrame - 1, width, nrTimePoints - calibrationFrame + 1);
			run("Plots...", "vertical");	//Make sure the profile is vertical
			calibrationMeanTrace = getProfile();
			Array.getStatistics(calibrationMeanTrace, min, max, calibrationMean, calibrationStdDev);
			calibrationMean = averageArray(calibrationMeanTrace);	//Get average without NaNs
			print("Average calibration lifetime, based on frames "+calibrationFrame +" - "+nrTimePoints+": " + calibrationMean + " ± " + calibrationStdDev + " ns");
		}
		else if(baselineOnly == true) stimulationFrame = 0; //Set back to 1 to prevent crashes later on
	}
	else if (singleFrame == true) {
		selectWindow(kymograph);
		//makeRectangle(0, 0, width, height);
		run("Select All");
		run("Plots...", "vertical");	//Make sure the profile is vertical
		baselineLength=0;
		baselineAverageTrace = getProfile();
		Array.getStatistics(baselineAverageTrace, min, max, baselineAverage, baselineStdDev);
		stimulationMeanTrace = getProfile();
		Array.getStatistics(stimulationMeanTrace, min, max, stimulationMean, stimulationStdDev);
		calibrationMeanTrace = getProfile();
		Array.getStatistics(calibrationMeanTrace, min, max, calibrationMean, calibrationStdDev);
	}

	//Measure pixel coordinates of all cells within the tile
	selectWindow(labelmap);
	Ext.CLIJ2_push(labelmap);
	run("Clear Results");
	Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(labelmap, labelmap);	//Possibility to measure the intensity as well - currently not implemented
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap, nr_cells);
	Ext.CLIJ2_release(labelmap);
	x_coords = Table.getColumn("CENTROID_X", "Results");	//N.B. First entry [0] is the background
	y_coords = Table.getColumn("CENTROID_Y", "Results");

	//Fill table with absolute coordinates of all cells of all tiles
	cellCoordinatesTable = "Cell_coordinates";
	if(!isOpen(cellCoordinatesTable)) {
		Table.create(cellCoordinatesTable);
		Table.setLocationAndSize(100, 100, 672, 1000);
		total_cells_all_tiles = 0;
	}
	else {
		selectWindow(cellCoordinatesTable);
		total_cells_all_tiles = Table.size;
	}
	cellNr_array = Array.getSequence(nr_cells+1);	//Create cell number array
	cellNr_array_total = addScalarToArray(cellNr_array, total_cells_all_tiles);
	for (i = 1; i < nr_cells+1; i++) {
		Table.set("Tile", total_cells_all_tiles + i-1, tile, cellCoordinatesTable);
		Table.set("Cell", total_cells_all_tiles + i-1, cellNr_array[i], cellCoordinatesTable);
		Table.set("Cell total", total_cells_all_tiles + i-1, cellNr_array_total[i], cellCoordinatesTable);
		Table.set("PosX (px)", total_cells_all_tiles + i-1, x_coords[i], cellCoordinatesTable);
		Table.set("PosY (px)", total_cells_all_tiles + i-1, y_coords[i], cellCoordinatesTable);
//		Table.set("AbsPosX (m)", total_cells_all_tiles + i-1, d2s(-parseFloat(stagePositions[2*tile+1]) - tileSizeX/2 + parseFloat((x_coords[i]*pixelWidth)/1E6),10), cellCoordinatesTable);
//		Table.set("AbsPosY (m)", total_cells_all_tiles + i-1, d2s(-parseFloat(stagePositions[2*tile]) - tileSizeY/2 + parseFloat((y_coords[i]*pixelHeight)/1E6),10), cellCoordinatesTable);
		Table.set("AbsPosX (m)", total_cells_all_tiles + i-1, d2s(-parseFloat(stagePositions[2*tile+1]) - tileSizeX/2 + parseFloat((x_coords[i]*pixelWidth)/1E6),10), cellCoordinatesTable);
		Table.set("AbsPosY (m)", total_cells_all_tiles + i-1, d2s(-parseFloat(stagePositions[2*tile]) - tileSizeY/2 + parseFloat((y_coords[i]*pixelHeight)/1E6),10), cellCoordinatesTable);
	}
	Table.update;

/*
	selectWindow(lifetimeTable);
	nrTimePoints = Table.size;
	headings = Table.headings(lifetimeTable);
	headers = split(headings, "\t");
*/

	//Correct response time window in case it falls outside the data boundaries. The margin is regarded as more important than the window.
	hit_avg_response_time_window = minOf(hit_avg_response_time_window, calibrationFrame - stimulationFrame);
	if(hit_avg_response_time_window == -1) {
		hit_avg_response_time_window = calibrationFrame - stimulationFrame;
		hit_avg_response_time_margin = 0;		
	}
	if(meanResponseTimingOption == "after stimulation" && (stimulationFrame + hit_avg_response_time_window + hit_avg_response_time_margin) > calibrationFrame) hit_avg_response_time_window = minOf(calibrationFrame - stimulationFrame - hit_avg_response_time_margin, calibrationFrame - stimulationFrame);
	else if(meanResponseTimingOption == "before calibration" && (calibrationFrame - hit_avg_response_time_margin - hit_avg_response_time_window) < stimulationFrame) hit_avg_response_time_window = maxOf(calibrationFrame - stimulationFrame - hit_avg_response_time_margin, 1);
	if(meanResponseTimingOption == "after stimulation") responseTimeStart = minOf(stimulationFrame + hit_avg_response_time_margin + 1, calibrationFrame);
	else if(meanResponseTimingOption == "before calibration") responseTimeStart = maxOf(calibrationFrame - hit_avg_response_time_margin - hit_avg_response_time_window + 1, stimulationFrame + 1);

	//Overrule when single frame
	if(singleFrame == true) {
		responseTimeStart = 0;
		hit_avg_response_time_window = 1;
		hit_avg_response_time_margin = 0;
	}
	if(baselineOnly == true) hit_avg_response_time_window = 0;
	if(hit_avg_response_boolean || hit_max_abs_response_boolean || hit_max_response_boolean || hit_rrr_boolean) print("Detecting response hits in a time window (frames): "+responseTimeStart+" - "+responseTimeStart + hit_avg_response_time_window - 1+".");
	if(hit_baseline_boolean) print("Detecting baseline hits in a time window (frames): 1 - "+baselineLength+".");
	if(hit_additional_channel_boolean && additionalChannel > 0) {
		if(hit_additional_channel_choice == "lower" || hit_additional_channel_choice == "higher") print("Detecting hits in additional channel ("+additionalChannel+"): "+hit_additional_channel_timing+" "+hit_additional_channel_choice+" than "+hit_additional_channel_number1+".");
		else if(hit_additional_channel_choice == "between") print("Detecting hits in additional channel ("+additionalChannel+"): "+hit_additional_channel_timing+" "+hit_additional_channel_choice+" "+hit_additional_channel_number1+" and "+hit_additional_channel_number2+".");
	}
	else if(hit_additional_channel_boolean) {
		print("WARNING: Find hits in additional channel is checked, but this parameter is not correctly set! (current value: "+additionalChannel+"). No hits will be registered.");
		hit_additional_channel_boolean = false;
	}

	//Preparations for fitting
	timeArray = Array.getSequence(nrTimePoints);
	timeArray = multiplyArraywithScalar(timeArray, frameInterval);

	timeArray_fit = Array.getSequence(hit_avg_response_time_window);
	timeArray_fit = addScalarToArray(timeArray_fit, responseTimeStart);
	timeArray_fit = multiplyArraywithScalar(timeArray_fit, frameInterval);

//	initialGuesses = newArray(2.3 0.7 200*frameInterval -0.05);	//Already defined in the dialog
	fitP0 = newArray(nr_cells);
	fitP1 = newArray(nr_cells);
	fitP2 = newArray(nr_cells);
	fitP3 = newArray(nr_cells);
	fitP4 = newArray(nr_cells);
	fitRSquare = newArray(nr_cells);
	fitPlotString = "";
	fitparametersTable = "Fit parameters - "+fit_equation;
	if(isOpen(fitparametersTable)) close(fitparametersTable);
	else if(fit_traces == true) Table.create(fitparametersTable);

	if(createNormalizedKymograph == true) {
		//Duplicate kymograph for normalization on baseline
		kymograph_baseline_normalized = kymograph+"_baseline_normalized";
		selectImage(kymograph);
		run("Select None");
		run("Duplicate...", "title=["+kymograph_baseline_normalized+"] ignore");
		setBatchMode("show");
	}

	//Prefill table with relevant headers
	Table.set("Cell", 0, "", hitList);
	if(hit_baseline_boolean == true) Table.set("mean baseline", nrHits, 0, hitList);
	if(hit_avg_response_boolean == true) Table.set("mean response", nrHits, 0, hitList);
	if(hit_max_abs_response_boolean == true) Table.set("max response", nrHits, 0, hitList);
	if(hit_max_response_boolean == true) Table.set("max response difference", nrHits, 0, hitList);
	if(hit_risetime_boolean == true) Table.set("rise time (frames)", nrHits, 0, hitList);
	if(hit_rrr_boolean == true) Table.set("rapid response ratio", nrHits, 0, hitList);
	if(hit_additional_channel_boolean == true) Table.set("mean additional channel intensity", nrHits, 0, hitList);
	Table.deleteRows(0, 0, hitList);
	Table.update(hitList);
	nrHits=0;

	if(generateRandomHits_boolean == false) {
		for(i=0; i<nr_cells; i++) {
	//		if(!reanalyze_boolean) {
				selectImage(kymograph);
				makeRectangle(i, 0, 1, kheight);
	//			lifetime_cell = Table.getColumn(headers[i+1], lifetimeTable);
				setKeyDown("alt");	//Make the profile vertical (stupid, but it works).
				lifetime_cell = getProfile();
				setKeyDown("none");
	//		}
	//		else lifetime_cell = Table.getColumn(headers[i], lifetimeTable);	//The first column (number) is not saved to disk 
			if(singleFrame == false) {
				lifetime_cell_baseline = Array.slice(lifetime_cell, 0, baselineLength);							//Array.slice does [...)
				if(calibrationOnly == false) lifetime_cell_response = Array.slice(lifetime_cell, baselineLength + 1, calibrationFrame);		//Array.slice does [...)
				else if (calibrationOnly == true) lifetime_cell_response = Array.slice(lifetime_cell, baselineLength, calibrationFrame-1);
			}
			else if(singleFrame == true) {
				lifetime_cell_baseline = lifetime_cell;
				lifetime_cell_response = lifetime_cell;
			}
			if(hit_avg_response_time_window != -1 && singleFrame == false) {	//User defined timing
	//			lifetime_cell_response_timed = Array.slice(lifetime_cell, baselineLength + 1 + responseTimeStart, baselineLength + 1 + responseTimeEnd);		//Array.slice does [...) 
				lifetime_cell_response_timed = Array.slice(lifetime_cell, responseTimeStart - 1, responseTimeStart + hit_avg_response_time_window - 1);		//Array.slice does [...) 
			}
			else lifetime_cell_response_timed = lifetime_cell_response;
	
			if(calibrationOnly == false) lifetime_cell_calibration = Array.slice(lifetime_cell, calibrationFrame, lifetime_cell.length);		//Array.slice does [...)
			else if(calibrationOnly == true) lifetime_cell_calibration = Array.slice(lifetime_cell, calibrationFrame-1, lifetime_cell.length);
			Array.getStatistics(lifetime_cell_baseline, min_baseline, max_baseline, mean_baseline, stdDev_baseline);
			Array.getStatistics(lifetime_cell_response, min_response, max_response, mean_response, stdDev_response);
			Array.getStatistics(lifetime_cell_response_timed, min_response_timed, max_response_timed, mean_response_timed, stdDev_response_timed);
			Array.getStatistics(lifetime_cell_calibration, min_calibration, max_calibration, mean_calibration, stdDev_calibration);
	//Array.show(lifetime_cell_baseline, lifetime_cell_response, lifetime_cell_calibration);		
			lifetime_cell_response_diff = addScalarToArray(lifetime_cell_response, -mean_baseline); //Subtract the cell's baseline from the response
			lifetime_cell_response_diff_timed = addScalarToArray(lifetime_cell_response_timed, -mean_baseline); //Subtract the cell's baseline from the response
			lifetime_cell_response_diff_to_avg_baseline = addScalarToArray(lifetime_cell_response, -baselineAverage); //Subtract the average baseline of all cells from the response
			Array.getStatistics(lifetime_cell_response_diff_timed, min_response_diff_timed, max_response_diff_timed, mean_response_diff_timed, stdDev_response_diff_timed);
			Array.getStatistics(lifetime_cell_response_diff, min_response_diff, max_response_diff, mean_response_diff, stdDev_response_diff);
			Array.getStatistics(lifetime_cell_response_diff_to_avg_baseline, min_response_diff_to_avg_baseline, max_response_diff_to_avg_baseline, mean_response_diff_to_avg_baseline, stdDev_response_diff_to_avg_baseline);
	
			//Get the first and second half of the response (for the rapid response ratio)
			lifetime_cell_response_diff_1 = Array.slice(lifetime_cell_response_diff, 0, lifetime_cell_response_diff.length/2);
			lifetime_cell_response_diff_2 = Array.slice(lifetime_cell_response_diff, lifetime_cell_response_diff.length/2, lifetime_cell_response_diff.length);
			Array.getStatistics(lifetime_cell_response_diff_1, min_response_diff_1, max_response_diff_1, mean_response_diff_1, stdDev_response_diff_1);
			Array.getStatistics(lifetime_cell_response_diff_2, min_response_diff_2, max_response_diff_2, mean_response_diff_2, stdDev_response_diff_2);
	
			//Calculate the slope of the response up to the maximum value (assuming linear rise)
	//max_response_lifetime_index = maxIndexOfArray(lifetime_cell_response_diff_timed);
	//slope = max_response_diff / (max_response_lifetime_index+1);
	//if(max_response_timed < 2.8) print("cell "+i+": "+max_response_timed, max_response_diff, max_response_lifetime_index, slope*100);
	//else print(slope*100);
	
			if(createNormalizedKymograph == true) {
				//Fill normalized kymograph
				tempArray = addScalarToArray(lifetime_cell, -mean_baseline);
				normalized_lifetime_cell = addScalarToArray(tempArray, baselineAverage);
				selectImage(kymograph_baseline_normalized);
				for (f=0; f<nrTimePoints; f++) setPixel(i, f, normalized_lifetime_cell[f]);
			}
			
			//Test if the cell is a hit
	
			//Check for valid traces
			validTrace = true;
			if(baseline_calibration_difference_boolean == true && ( (mean_calibration - mean_baseline) < baseline_calibration_difference_number_low || (mean_calibration - mean_baseline) > baseline_calibration_difference_number_high) ) validTrace = false;
			if(mean_calibration - mean_baseline == NaN) validTrace = false;
			if(baseline_avg_baseline_difference_boolean == true && mean_baseline - baselineAverage > baseline_avg_baseline_difference_number) validTrace = false;
			if(mean_calibration - mean_baseline == NaN) print("WARNING: NaN found in the data. Average values are NaN as well.");

			
			potential_hit = false;
			if(validTrace == true) {		
				//Mean baseline
				if(hit_baseline_boolean == true) {
					if(hit_baseline_choice == "lower") {
						if(mean_baseline < hit_baseline_number) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with baseline lifetime "+mean_baseline+ " ns.");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("mean baseline", nrHits, mean_baseline, hitList);
							potential_hit = true;
						}

					}
					else if(hit_baseline_choice == "higher") {
						if(mean_baseline > hit_baseline_number) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with baseline lifetime "+mean_baseline+ " ns.");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("mean baseline", nrHits, mean_baseline, hitList);
							potential_hit = true;
						}
					}
				}

				//Mean response
				if(absoluteOrRelativeOption == "mean response lifetime difference with baseline") {
					if(hit_avg_response_boolean == true) {
						if(hit_avg_response_choice == "lower") {
							if(mean_response_diff < hit_avg_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean response lifetime difference "+mean_response_diff+ " ns ("+mean_response+" - "+mean_baseline+").");
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("mean response", nrHits, mean_response_diff, hitList);
								potential_hit = true;
							}

						}
						else if(hit_avg_response_choice == "higher") {
							if(mean_response_diff > hit_avg_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean response lifetime difference "+mean_response_diff+ " ns ("+mean_response+" - "+mean_baseline+").");
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("mean response", nrHits, mean_response_diff, hitList);
								potential_hit = true;
							}
						}
					}
				}
				else if(absoluteOrRelativeOption == "mean absolute response lifetime") {
					if(hit_avg_response_boolean == true) {
						if(hit_avg_response_choice == "lower") {
							if(mean_response_timed < hit_avg_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean response lifetime "+mean_response_timed+ " ns.");
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("mean response", nrHits, mean_response, hitList);								
								potential_hit = true;
							}
	
						}
						else if(hit_avg_response_choice == "higher") {
							if(mean_response_timed > hit_avg_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean response lifetime "+mean_response_timed+ " ns.");
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("mean response", nrHits, mean_response, hitList);	
								potential_hit = true;
							}
						}
					}
				}
				else if(absoluteOrRelativeOption == "fraction of max response minus baseline") {
					if(hit_avg_response_boolean == true) {
						if(hit_avg_response_choice == "lower") {
							if((mean_response_diff_timed / max_response_diff) < hit_avg_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean fraction of max response difference with baseline "+mean_response_diff_timed / max_response_diff+" ("+mean_response_timed+" / "+max_response_diff+")"+".");
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("mean response", nrHits, mean_response, hitList);								
								potential_hit = true;
							}

						}
						else if(hit_avg_response_choice == "higher") {
							if((mean_response_diff_timed / max_response_diff) > hit_avg_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean fraction of max response difference with baseline "+mean_response_diff_timed / max_response_diff+" ("+mean_response_timed+" / "+max_response_diff+")"+".");
		print(mean_response_diff_timed, max_response_diff, hit_avg_response_number * max_response_diff);
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("mean response", nrHits, mean_response, hitList);	
								potential_hit = true;
							}
						}
					}
				}

	//
	//			//absolute mean response lifetime
	//			if(hit_avg_abs_response_boolean == true) {
	//				if(hit_avg_abs_response_choice == "lower") {
	//						if(mean_response_timed < hit_avg_abs_response_number) {
	//							print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean absolute response lifetime "+mean_response_timed+".");
	//							register_hit(hitList, nrHits, i, tile, x_coords[i+1], y_coords[i+1], mean_response_timed);
	//						nrHits++;
	//					}
	//				}
	//				else if(hit_avg_abs_response_choice == "higher") {
	//						if(mean_response_timed > hit_avg_abs_response_number) {
	//							print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with max absolute response lifetime difference "+mean_response_timed+".");
	//							register_hit(hitList, nrHits, i, tile, x_coords[i+1], y_coords[i+1], mean_response_timed);
	//						nrHits++;
	//					}
	//				}
	//			}
	
	
				//Maximum absolute response lifetime
				if(hit_max_abs_response_boolean == true) {
					if(hit_max_abs_response_choice == "lower") {
						if(max_response < hit_max_abs_response_number) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with max absolute response lifetime difference "+max_response+".");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("max response", nrHits, max_response, hitList);	
							potential_hit = true;
						}
					}
					else if(hit_max_abs_response_choice == "higher") {
						if(max_response > hit_max_abs_response_number) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with max absolute response lifetime difference "+max_response+".");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("max response", nrHits, max_response, hitList);	
							potential_hit = true;
						}
					}
				}

				//Max difference to baseline
				if(hit_max_response_boolean == true) {
					if(hit_max_response_to_avg_baseline_boolean == false) {
						if(hit_max_response_choice == "lower") {
							if(max_response_diff < hit_max_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with max response lifetime difference "+max_response_diff+ " ns ("+max_response+" - "+mean_baseline+").");
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("max response difference", nrHits, max_response_diff, hitList);
								potential_hit = true;
							}
						}
						else if(hit_max_response_choice == "higher") {
							if(max_response_diff > hit_max_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with max response lifetime difference "+max_response_diff+ " ns ("+max_response+" - "+mean_baseline+").");
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("max response difference", nrHits, max_response_diff, hitList);
								potential_hit = true;
							}
						}
					}
					else if(hit_max_response_to_avg_baseline_boolean == true){	//Max difference to *average* baseline (of all cells)
						if(hit_max_response_choice == "lower") {
							if(max_response_diff_to_avg_baseline < hit_max_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with max response lifetime difference to average baseline "+max_response_diff_to_avg_baseline+ " ns ("+max_response+" - "+baselineAverage+").");
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("max response diff to avg baseline", nrHits, max_response_diff_to_avg_baseline, hitList);	
								potential_hit = true;
							}
						}
						else if(hit_max_response_to_avg_baseline_choice == "higher") {
							if(max_response_diff_to_avg_baseline > hit_max_response_number) {
								if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with max response lifetime difference to average baseline "+max_response_diff_to_avg_baseline+ " ns ("+max_response+" - "+baselineAverage+").");
								Table.set("Cell", nrHits, i+1, hitList);
								Table.set("max response diff to avg baseline", nrHits, max_response_diff_to_avg_baseline, hitList);	
								potential_hit = true;
							}
						}
					}
				}

				//Rise time after stimulation
				if(hit_risetime_boolean == true) {
					z=0;
					while(lifetime_cell_response_diff[z] < max_response_diff * risetime_fraction) {	//currently set to 75% of max reponse value
						z++;
						if(z == lifetime_cell_response_diff.length-1) break;
					}
					if(hit_risetime_choice == "slower") {
						if(z > hit_risetime_number-1) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with rise time "+z+1+ " frames (span of response (max-baseline): "+max_response+" - "+mean_baseline+" ns).");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("rise time (frames)", nrHits, z+1, hitList);	
							potential_hit = true;
						}
					}
					if(hit_risetime_choice == "faster") {
						if(z < hit_risetime_number-1) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with rise time "+z+1+ " frames (span of response (max-baseline): "+max_response+" - "+mean_baseline+" ns).");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("rise time (frames)", nrHits, z+1, hitList);	
							potential_hit = true;
						}
					}
				}
	
				//Rapid response ratio (mean second half / mean first half of the response)
				if(hit_rrr_boolean == true) {
					rrr = mean_response_diff_2 / mean_response_diff_1;
					if(hit_rrr_choice == "lower") {
						if(rrr < hit_rrr_number) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with rapid response ratio "+rrr+".");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("rapid response ratio", nrHits, rrr, hitList);	
							potential_hit = true;
						}
					}
					else if(hit_rrr_choice == "higher") {
						if(rrr > hit_rrr_number) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with rapid response ratio "+rrr+".");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("rapid response ratio", nrHits, rrr, hitList);	
							potential_hit = true;
						}
					}
				}

				//Mean intensity in the additional channel
				if(hit_additional_channel_boolean == true && additionalChannel != -1) {
					selectImage(kymograph_additionalChannel);
					if(hit_additional_channel_timing == "Intensity in the first frame") {
						mean_additional_channel_intensity = getPixel(i, 0);
					}
					else if (hit_additional_channel_timing == "Average intensity") {
						makeRectangle(i, 0, 1, kheight);
						mean_additional_channel_intensity = getValue("Mean");
					}
					Table.set("intensity ch"+additionalChannel, i, mean_additional_channel_intensity, cellCoordinatesTable);
					if(hit_additional_channel_choice == "lower") {
						if(mean_additional_channel_intensity < hit_additional_channel_number1) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean intensity "+mean_additional_channel_intensity+".");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("mean additional channel intensity", nrHits, mean_additional_channel_intensity, hitList);	
							potential_hit = true;
						}
					}
					else if(hit_additional_channel_choice == "higher") {
						if(mean_additional_channel_intensity > hit_additional_channel_number1) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean intensity "+mean_additional_channel_intensity+".");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("mean additional channel intensity", nrHits, mean_additional_channel_intensity, hitList);	
							potential_hit = true;
						}
					}
					else if(hit_additional_channel_choice == "between") {
						if(mean_additional_channel_intensity > hit_additional_channel_number1 && mean_additional_channel_intensity < hit_additional_channel_number2) {
							if(debugMode) print("Potential hit found: Tile "+tile+", cell "+i+1+ ", with mean intensity "+mean_additional_channel_intensity+".");
							Table.set("Cell", nrHits, i+1, hitList);
							Table.set("mean additional channel intensity", nrHits, mean_additional_channel_intensity, hitList);	
							potential_hit = true;
						}
					}
					run("Select None");
				}

//				if(hit_measurements_.length > 0) {
//					if(sort_hits_column != "do not sort") measurement = hit_measurements_[(sort_hits_column)-1];	//Fill the hit_measurement_ entry with the measurement to be sorted on
//					register_hit(hitList, nrHits, i, tile, x_coords[i+1], y_coords[i+1], measurement);
//				}
			}	//End of hit conditions loop
			
			if(fit_traces == true) {
				showStatus("Fitting cell "+i+1+"/"+nr_cells);
				showProgress(i, nr_cells);
				if(initialGuesses.length != 0) Fit.doFit(fit_equation, timeArray_fit, lifetime_cell_response_timed, initialGuesses);
				else Fit.doFit(fit_equation, timeArray_fit, lifetime_cell_response_timed);
				fitP0[i] = Fit.p(0);
				if(nrFitParameters > 1) fitP1[i] = Fit.p(1);
				if(nrFitParameters > 2) fitP2[i] = Fit.p(2);
				if(nrFitParameters > 3) fitP3[i] = Fit.p(3);
				if(nrFitParameters > 4) fitP4[i] = Fit.p(4);
				fitRSquare[i] = Fit.rSquared;
				
				//Plotting the fits takes a long time, so is currently disabled
	/*			Fit.plot;
				rename("fit");
				Plot.getValues(xfit, yfit);	//stupid but fast way to create fit data
				close("fit");
				Plot.create("Fit cell "+i+1, "time (s)", "counts");
				Plot.setFrameSize(900, 600);
				Plot.setLineWidth(1);
				Plot.setColor("blue");
				Plot.add("circle", timeArray, lifetime_cell);
				Plot.setColor("red");
				Plot.add("line", xfit, yfit);
				Plot.setLineWidth(1);		
				Plot.setLimits(NaN, NaN, minLifetime, maxLifetime);
				Plot.setColor("black");
				textposX = 0.21;
				textposY = 0.79;
				lineDistance = 0.03;
				Plot.addText(fit_equation, textposX, textposY);
				Plot.addText("a =  "+fitP0[i], textposX, textposY + lineDistance);
				Plot.addText("b =  "+fitP1[i], textposX, textposY + lineDistance*2);
	//			Plot.addText("a+b (max value)=  "+fitP0[i] + fitP1[i], textposX, textposY + lineDistance*3);
				Plot.addText("c =  "+fitP2[i], textposX, textposY + lineDistance*3);
				Plot.addText("d =  "+fitP3[i], textposX, textposY + lineDistance*4);
				Plot.addText("e =  "+fitP4[i], textposX, textposY + lineDistance*5);
				Plot.addText("R2 = "+fitRSquare[i], textposX, textposY + lineDistance*6);
				Plot.show();
				rename("Fit cell "+i+1);
				fitPlotString += " image"+i+1+"=[Fit cell "+i+1+"]";	*/
			}

			if(potential_hit == true) {
				register_hit(hitList, nrHits, i+1, tile,  x_coords[i+1], y_coords[i+1]);
				nrHits++;
			}
		}	//End of loop over all cells

		if(nrHits > 0) {
			Table.update(hitList);
			headings = split(Table.headings(hitList), "\t");
			
			//Remove hits if 'AND' is not satisfied
			if(hit_find_logic == "AND") {
				nrCriteria = firstIndexOfArray(headings, "Tile") - 2;	//The first column is empty, the second one is 'Cell'
				total = nrHits;
				print(total+" potential hits discovered.");
				for(k=total-1; k>=0; k--) {
					showStatus("Purging false hits...");
					showProgress(total-1-k, total-1);
					OK = 0;
					for(n=2; n<nrCriteria+2; n++) {
						if(Table.get(headings[n], k, hitList) != 0) OK++;
					}
					if(OK != nrCriteria) {
						Table.deleteRows(k, k);
						nrHits--;
					}
				}
				print(total - nrHits+" were removed because they don't satisfy all criteria.");
			}

			//Sort the hits
			if(occursInArray(headings, sort_hits_column)) {
				Table.sort(sort_hits_column, hitList);
				if(sort_hits_direction == "highest first") reverse_table(hitList);
			}
			else if(sort_hits_column != "do not sort") print("[WARNING] Column '"+sort_hits_column+"' not found in Potential hits table. Cannot sort the hits!");

		}
	}
	else if(generateRandomHits_boolean == true) {
		sequentialArray = Array.getSequence(nr_cells);
		randomArray = shuffle_array(sequentialArray);
		randomHits = Array.trim(randomArray, nrOfRandomHits);
		print("Generating "+nrOfRandomHits+" random hits");
		for (i = 0; i < nrOfRandomHits; i++) {
			Table.set("Cell", nrHits, randomHits[i], hitList);
			register_hit(hitList, nrHits, randomHits[i]-1, tile, x_coords[randomHits[i]], y_coords[randomHits[i]]);
			nrHits++;
		}
	}

	if(createNormalizedKymograph == true) {
		selectImage(kymograph_baseline_normalized);
		updateDisplay();
	}
	
	//Display the additional channel average intensity
	if(hit_additional_channel_boolean == true && additionalChannel > 0) {
		if(hit_additional_channel_timing == "Average intensity" && singleFrame == false) {
			selectWindow(input_image);
			run("Duplicate...", "title=Channel"+additionalChannel+" duplicate channels="+additionalChannel);
			run("Z Project...", "projection=[Average Intensity]");
			rename("Channel"+additionalChannel+" average intensity");
			setBatchMode("show");
			roiManager("Show all without labels");
			saveAs("tiff", output + saveName + " (Average intensity channel "+additionalChannel+")");
			close("Channel"+additionalChannel);
		}
	}
	Table.update(cellCoordinatesTable);
	
/*	if(fit_traces == true) {
		run("Concatenate...", "  title=fit_stack open "+fitPlotString);
		for(i=0; i<nr_cells; i++) {
			Stack.setFrame(i+1);
			run("Set Label...", "label=[cell "+i+1+"]");
		}
		setBatchMode("show");
*/		
	if(fit_traces == true) {
		selectWindow(fitparametersTable);
		sequenceArray = addScalarToArray(Array.getSequence(nr_cells), 1);
		Table.setColumn("cell nr", sequenceArray);
		Table.setColumn("a", fitP0);
		if(nrFitParameters > 1) Table.setColumn("b", fitP1);
		if(nrFitParameters > 2) Table.setColumn("c", fitP2);
		if(nrFitParameters > 3) Table.setColumn("d", fitP3);
		if(nrFitParameters > 4) Table.setColumn("e", fitP4);
		Table.setColumn("R2", fitRSquare);
		Table.update;
	}
	Table.update(hitList);
	nrHits = Table.size(hitList);
	print("Screening result: "+nrHits+" hits, in a total of "+nr_cells+" cells. ("+d2s(nrHits/nr_cells*100,1)+" %)");
/*	
	Plot.create("Midpoint", "time (s)", "counts");
	Plot.addHistogram(fitP2, 50);
	Plot.show();
	setBatchMode("show");
*/
	return hitList;
}


function register_hit(hitList, nrHits, cell, tile, x, y) {
//	Table.set("nr", nrHits, nrHits, hitList);
	Table.set("Tile", nrHits, tile, hitList);
//	Table.set("Cell", nrHits, cell+1, hitList);
	Table.set("X (um)", nrHits, x*pixelWidth, hitList);
	Table.set("Y (um)", nrHits, y*pixelWidth, hitList);
	Table.set("TilePosX", nrHits, d2s(-parseFloat(stagePositions[2*tile+1]) - tileSizeX/2, 10), hitList);	//X and Y seem to be reversed!
	Table.set("TilePosY", nrHits, d2s(-parseFloat(stagePositions[2*tile]) - tileSizeY/2, 10), hitList);		//X and Y seem to be reversed!
	Table.set("AbsPosX (m)", nrHits, d2s(-parseFloat(stagePositions[2*tile+1]) - tileSizeX/2 + parseFloat((x*pixelWidth)/1E6),10), hitList);
	Table.set("AbsPosY (m)", nrHits, d2s(-parseFloat(stagePositions[2*tile]) - tileSizeY/2 + parseFloat((y*pixelHeight)/1E6),10), hitList);
}


function assign_nuclei_to_cells(labelmap_cells, labelmap_nuclei) {
	Ext.CLIJ2_push(labelmap_cells);
	Ext.CLIJ2_push(labelmap_nuclei);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_cells, nr_cells);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap_nuclei, nr_nuclei);
	print("nr of cells: "+nr_cells);
	print("nr of nuclei: "+nr_nuclei);

	//Compute the overlap index of each cell with each nucleus. Pair every cell to the nucleus with the largest overlap.
	Ext.CLIJ2_generateJaccardIndexMatrix(labelmap_cells, labelmap_nuclei, jaccard_index_matrix);
	Ext.CLIJ2_transposeXZ(jaccard_index_matrix, jaccard_index_matrix_transposed);
	Ext.CLIJ2_argMaximumZProjection(jaccard_index_matrix_transposed, max_overlap, index_max_overlap);	//The index of the maximum is the nucleus with the largest overlap
	Ext.CLIJ2_transposeXY(index_max_overlap, index_max_overlap_transposed);
	labelmap_reassigned_nuclei = "labelmap_reassigned_nuclei";
	Ext.CLIJ2_replaceIntensities(labelmap_nuclei, index_max_overlap_transposed, labelmap_reassigned_nuclei);
	Ext.CLIJ2_pull(labelmap_reassigned_nuclei);
	setBatchMode("show");
	Ext.CLIJ2_release(jaccard_index_matrix_transposed);
	Ext.CLIJ2_release(jaccard_index_matrix);
	Ext.CLIJ2_release(index_max_overlap);
	Ext.CLIJ2_release(max_overlap);
	Ext.CLIJ2_release(index_max_overlap_transposed);

	//Find missing values (cells without nuclei) and add them to the labelmap
	//Create binary array of cells without a nucleus (empty cells)
	selectWindow(labelmap_reassigned_nuclei);
	run("Glasbey on dark");
	getRawStatistics(nPixels, mean, min, max, std, histogram);	//histogram size is equal to the nucleus with largest label that overlaps with a cell.
	empty_cells = Array.getSequence(nr_cells);
	for (i = 0; i < nr_cells; i++) {
		if(i<=max) {
			if(histogram[i] != 0) empty_cells[i] = 0;
		}
	}
	//Array.print(empty_cells);
	//print(empty_cells.length);
	
	Ext.CLIJ2_pushArray(empty_cells_image, empty_cells, nr_cells, 1, 1);
	labelmap_empty_cells = "labelmap_empty_cells";
	Ext.CLIJ2_replaceIntensities(labelmap_cells, empty_cells_image, labelmap_empty_cells);
	Ext.CLIJ2_pull(labelmap_empty_cells);
	run("Glasbey on dark");
	setMinAndMax(0, nr_cells);
	setBatchMode("show");
	
	//Mask part of empty cells that overlap with nuclei
	Ext.CLIJ2_mask(labelmap_reassigned_nuclei, labelmap_empty_cells, masked);
	Ext.CLIJ2_binaryIntersection(labelmap_reassigned_nuclei, labelmap_empty_cells, mask);
	Ext.CLIJ2_binaryNot(mask, mask_inverted);
	Ext.CLIJ2_mask(labelmap_empty_cells, mask_inverted, labelmap_empty_cells_masked);
	
	//Add reassigned nuclei and empty cells labelmaps
	labelmap_reassigned_nuclei_and_empty_cells = "labelmap_reassigned_nuclei_and_empty_cells";
	Ext.CLIJ2_addImages(labelmap_reassigned_nuclei, labelmap_empty_cells_masked, labelmap_reassigned_nuclei_and_empty_cells);
	Ext.CLIJ2_pull(labelmap_reassigned_nuclei_and_empty_cells);
	run("Glasbey on dark");
	setMinAndMax(0, nr_cells);
	setBatchMode("show");
	
	Ext.CLIJ2_clear();
	return labelmap_reassigned_nuclei_and_empty_cells;
}


function labels_to_ROI_Manager(labelmap) {
//	selectWindow(labelmap);
	Ext.CLIJ2_push(labelmap);
//	run("Clear Results");
//	Ext.CLIJ2_statisticsOfLabelledPixels(labelmap, labelmap);	//This is in fact redundant, because the statistics are already present.
	boundingBox_X = Table.getColumn("BOUNDING_BOX_X", "Cell_statistics");
	boundingBox_Y = Table.getColumn("BOUNDING_BOX_Y", "Cell_statistics");
	boundingBox_width = Table.getColumn("BOUNDING_BOX_WIDTH", "Cell_statistics");
	boundingBox_height = Table.getColumn("BOUNDING_BOX_HEIGHT", "Cell_statistics");
	Array.getStatistics(boundingBox_width, min, boundingBoxMax_X, mean, stdDev);
	Array.getStatistics(boundingBox_height, min, boundingBoxMax_Y, mean, stdDev);
	Ext.CLIJ2_getMaximumOfAllPixels(labelmap, nr_cells);

	if(isOpen("ROI Manager")) { selectWindow("ROI Manager"); }//run("Close"); }	//This step goes faster when the ROI manager is not visible.
	for (i = 0; i < nr_cells; i++) {
		showStatus("Converting "+nr_cells+" labels to ROIs...");
		Ext.CLIJ2_crop2D(labelmap, label_cropped, boundingBox_X[i], boundingBox_Y[i], boundingBoxMax_X, boundingBoxMax_Y);
		Ext.CLIJ2_labelToMask(label_cropped, mask_label, i+1);
		Ext.CLIJ2_pullToROIManager(mask_label);
		roiManager("Select",i);
		Roi.move(boundingBox_X[i], boundingBox_Y[i]);
		roiManager("update");
//		else {	//Else the label is square and doesn't contain a zero, causing a crash because it is not added to the ROI manager.
//			//TO DO: ReplaceIntensities on the labelmap in GPU memory, closeIndexGapsInLabelMap and, ultimately, pull and overwrite the final labelmap
//		}
	}
	run("Select None");
	Ext.CLIJ2_release(mask_label);
	Ext.CLIJ2_release(label_cropped);
	roiManager("Deselect");
	roiManager("Remove Frame Info");
}


function measure_intensity_in_additional_channel(image, additionalChannel) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
	run("Duplicate...", "title=["+saveName+"_ch"+additionalChannel+"] duplicate channels="+additionalChannel);

	kymograph_additionalChannel = measure_intensities(saveName+"_ch"+additionalChannel, nr_cells, generateTables);
	maxIntensity = getValue("Max");
	if(generateTables) {
//		for (i = 0; i < nr_cells; i++) Table.renameColumn("Mean(cell_"+i+1+")", "cell_"+IJ.pad(i+1,5));
		if(frames > 1) {
			timeArray = Array.getSequence(frames);
			timeArray = multiplyArraywithScalar(timeArray, frameInterval);
			Table.setColumn("time (s)", timeArray);
		}
		Table.update;
		Table.save(output + saveName + "_intensity_ch"+additionalChannel+".tsv");
	}
	
//	if(isOpen("Results")) { selectWindow("Results"); run("Close"); }	//If the Results window is open measurements take longer!
//	run("Set Measurements...", "mean redirect=None decimal=3");
//	roiManager("deselect");
//	showStatus("Measuring intensities...");
//	roiManager("Multi Measure");
//	selectWindow("Results");
//
//	//Get maximum intensity for Plot limits
//	intensities_additional_channel = "intensities_additional_channel";
//	Ext.CLIJ2_pushResultsTable(intensities_additional_channel);
//	Ext.CLIJ2_getMaximumOfAllPixels(intensities_additional_channel, maxIntensity);
//	Ext.CLIJ2_release(intensities_additional_channel);
//
//	Table.rename("Results", "Intensity_table_ch"+additionalChannel);

	close(saveName+"_ch"+additionalChannel);
	return newArray(kymograph_additionalChannel, maxIntensity);
}


function measure_lifetime_traces(intensity_stack, lifetime_stack, labelmap_cells, saveName, nr_cells) {

	//Create masked intensity stack in case not all pixels in the lifetime image have a value
	Ext.CLIJ2_push(intensity_stack);
	Ext.CLIJ2_push(lifetime_stack);
	intensity_stack_masked = "intensity_masked";
	Ext.CLIJ2_mask(intensity_stack, lifetime_stack, intensity_stack_masked);
	Ext.CLIJ2_release(intensity_stack);
	Ext.CLIJ2_pull(intensity_stack_masked);
	run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	setBatchMode("show");
	
	if(microscope == "Frequency Domain FLIM") frames = frames/phases;
	
	run("Enhance Contrast", "saturated=0.35");

//N.B. Unfortunately we cannot have NaNs any more, since CLIJ2 doesn't support measuring only non-NaN pixels. 
/*
	for(i=1;i<=frames;i++) {
		if(frames>1) Stack.setFrame(i);
		showStatus("Converting zeros to NaNs in intensity image...");
		changeValues(0,0,NaN);	//Set zeroes to NaN in intensity stack
		showProgress(i, frames);
	}
*/
	run("Set Measurements...", "area mean redirect=None decimal=3");

	//Convert labelmap to ROIs
	selectImage(labelmap_cells);
	getLut(reds, greens, blues);
setBatchMode(false);	//Somehow it doesn't want to select this image in Batch Mode
selectImage(labelmap_cells);
setBatchMode(true);
	if(!speedup) labels_to_ROI_Manager(labelmap_cells);
	//run("NKI Labelmap to ROI Manager");
	else if (speedup) {
		run("Label image to ROIs", "rm=[RoiManager[visible=true]]");	//fast BIOP plugin, but composite ROIs from Cellpose mess it up
		run("Clear Results");
		roiManager("Measure");
		ROIsize = Table.getColumn("Area", "Results");
		for (i = 0; i < ROIsize.length; i++) {
			if(ROIsize[i] < minCellSize) {
				roiManager("Select", i);
				run("Clear", "slice");
			}
		}
		run("Select None");
		run("Label image to ROIs", "rm=[RoiManager[visible=true]]");	//fast BIOP plugin, but single pixel ROIs from Cellpose mess it up
		max = getValue("Max");
		if(max == roiManager("count")) print("Nr. of labels matches nr. of ROIs.");
		else {
			waitForUser("WARNING!  Nr. of labels ("+max+") does not match nr. of ROIs. ("+roiManager("count")+"). Press OK to continue.");
			nr_cells = max;	//TO DO: CHECK if this actually works out!
		}
	}
	setBatchMode("show");

//	nr_cells = roiManager("count");
//setBatchMode(false);
selectImage(labelmap_cells);
Ext.CLIJ2_push(labelmap_cells);
setBatchMode(true);
//run("Label image to ROIs", "rm=[RoiManager[visible=true]]");
//Or:
//run("Label image to composite ROIs", "rm=[RoiManager[visible=true]]");
	//measure and save intensities
	selectWindow(intensity_stack_masked);
//	if(isOpen("Results")) { selectWindow("Results"); run("Close"); }	//If the Results window is open measurements take longer!
	run("Set Measurements...", "mean redirect=None decimal=3");
	selectWindow(intensity_stack_masked);
	roiManager("deselect");

	showStatus("Measuring intensities...");
	intensity_kymograph = measure_intensities(intensity_stack_masked, nr_cells, generateTables);
	if(generateTables) Table.save(output + saveName + "_intensity.tsv");

	//multiply lifetime with intensity (for normalization)
	imageCalculator("Multiply create 32-bit stack", lifetime_stack, intensity_stack_masked);
	rename("Lifetime_times_intensity");
	
	Ext.CLIJ2_release(intensity_stack_masked);
	close(intensity_stack_masked);

	//measure lifetimes times intensity
	if(isOpen("Results")) { selectWindow("Results"); run("Close"); }
	run("Set Measurements...", "mean redirect=None decimal=3");
	selectWindow("Lifetime_times_intensity");
	roiManager("deselect");
	showStatus("Measuring lifetime*intensity...");
	lifetimeTimesIntensity_kymograph = measure_intensities("Lifetime_times_intensity", nr_cells, false);

	//Normalize lifetime with cell intensity and create tables and plot
	kymograph = "kymograph";
	showStatus("Computing lifetimes...");
	imageCalculator("Divide create 32-bit stack", lifetimeTimesIntensity_kymograph, intensity_kymograph);
	rename(kymograph);
//	getDimensions(kwidth, kheight, channels, slices, kframes);
//	nrTimePoints = kheight;
	
	Ext.CLIJ2_push(kymograph);
	run("Clear Results");
	kymograph_smoothed = "kymograph smoothed";
	if(smoothRadiusTraces >0) {
		Ext.CLIJ2_mean2DBox(kymograph, kymograph_smoothed, 0, smoothRadiusTraces);
		Ext.CLIJ2_pull(kymograph_smoothed);
		setMinAndMax(minLifetime, maxLifetime);
		setBatchMode("show");
		if(generateTables) Ext.CLIJ2_pullToResultsTable(kymograph_smoothed);
	}
	else {
		if(generateTables) {
			Ext.CLIJ2_pullToResultsTable(kymograph);
			
		}
//		Ext.CLIJ2_release(kymograph_transposed);
	}

	//Construct plot of all traces, or a histogram if there is only one time point.
	if(microscope == "Confocal TCSPC / TauSeparation") 	y_axis = "Lifetime (ns)";
	else if(microscope == "TauContrast") 				y_axis = "Lifetime (ns)";
	else if(microscope == "Fast FLIM") 					y_axis = "Lifetime (ns)";
	else if(microscope == "Frequency Domain FLIM")		y_axis = "Lifetime (ns)";
	else if(microscope == "Ratio Imaging")				y_axis = "Ratio (ch1 / ch2)";
	else if(microscope == "Intensity only")				y_axis = "Intensity";

	if(frames > 1) {
		timeArray = Array.getSequence(frames);
		timeArray = multiplyArraywithScalar(timeArray, frameInterval);
		plotName = saveName + " (lifetime traces plot)";
		Plot.create(plotName, "time (s)", y_axis);
		if(microscope != "Intensity only") Plot.setLimits(0, frames*frameInterval, minLifetime, maxLifetime);
		Plot.setFrameSize(900, 600);
		Plot.setAxisLabelSize(axisFontSize);
		Plot.setFontSize(axisFontSize);
		//Plot.setXYLabels("Distance (microns)", "Gray Value");
	}
	else {
		plotName = saveName + " (lifetime histogram plot)";
		Plot.create(plotName, y_axis, "count");
		Plot.setFrameSize(640, 480);
		//Plot.setLimits(0, frames*frameInterval, minLifetime, maxLifetime);
	}

	selectWindow(kymograph);
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	run("Rotate 90 Degrees Left");
	run("Flip Vertically");
	setBatchMode("hide");
	if(smoothRadiusTraces > 0) {
		selectWindow(kymograph_smoothed);
		run(lut);
		run("Rotate 90 Degrees Left");
		run("Flip Vertically");
		setBatchMode("hide");
	}
	for(i=0;i<nr_cells;i++) {
		showStatus("Computing lifetime... "+i+"/"+nr_cells);
		showProgress(i/nr_cells);
		makeRectangle(0, i, frames, 1);
		lifetimeData = getProfile();

		if(frames > 1) {
			//Generate traces with 'glasbey on dark' colors
			color = getLabelColor(i, nr_cells);
			Plot.setColor(color);
			Plot.add("line", timeArray, lifetimeData);
		}
		roiManager("select",i);
		roiManager("rename", "cell_"+i+1);
//		roiManager("Set Color", "#"+color1+color2+color3);	//Works, but looks very confusing
	}
	run("Select None");

	//Calculate average trace
	run("Select All");
	run("Plots...", "minimum=0");	//Make sure the profile is horizontal
	average_lifetime_trace = getProfile();
	run("Select None");

	//Calculate lifetime standard deviation for all cells at all timepoints
	Ext.CLIJ2_transposeXZ(kymograph, kymograph_XZ);
	Ext.CLIJ2_standardDeviationZProjection(kymograph_XZ, kymograph_stddev_XZ);
	Ext.CLIJ2_transposeXY(kymograph_stddev_XZ, kymograph_stddev);
	Ext.CLIJ2_pull(kymograph_stddev);
	Ext.CLIJ2_release(kymograph_XZ);
	Ext.CLIJ2_release(kymograph_stddev_XZ);
	Ext.CLIJ2_release(kymograph_stddev);
	selectWindow(kymograph_stddev);
	run("Select All");
	stddev_avg_lifetime_trace = getProfile();
	confidence_interval = multiplyArraywithScalar(stddev_avg_lifetime_trace, confidence_interval_sigma);

	//Add average trace to the plot and save
	if(frames > 1) {
		Plot.setColor("black");
		Plot.setLineWidth(3);
		Plot.add("line", timeArray, average_lifetime_trace);
		Plot.setLineWidth(1);
		average_trace_confidence_max = addArrays(average_lifetime_trace, confidence_interval);
		average_trace_confidence_min = subtractArrays(average_lifetime_trace, confidence_interval);
		if(displayConfidenceInterval == true) {
			Plot.setColor("#dddddd");
			Plot.setLineWidth(3);
			Plot.add("line", timeArray, average_trace_confidence_max);
			Plot.add("line", timeArray, average_trace_confidence_min);
		}
		if(microscope != "Intensity only") Plot.setLimits(0, frames*frameInterval, minLifetime, maxLifetime);
		else Plot.setLimitsToFit();
		if(displayGrid == true) Plot.setFormatFlags("11000000111111");
		else Plot.setFormatFlags("11000000001111");
		Plot.setColor("#003399");
		Plot.setFontSize(axisFontSize*0.8);
		Plot.setJustification("left");
		Plot.addText(nr_cells+" traces", 0.02, 0.05);
		Plot.addText("smoothing: "+smoothRadiusTraces+ " frames", 0.02, 0.08);
		Plot.setFontSize(axisFontSize);
	}
	else {	//Create histogram instead of time trace plot
		selectWindow(kymograph);
		getStatistics(area, mean, min, max, std);
		getHistogram(values, counts, histogramBins, minLifetime, maxLifetime);
		Plot.add("bar", values, counts);
		Plot.setStyle(0, "#0000a0, #00a0ff, 1.0, Separated Bars");
		Plot.setJustification("right");
		Plot.setFontSize(axisFontSize);
		Plot.addText("Mean ± sd: "+d2s(mean,2)+" ± "+d2s(std,2), 0.97, 0.08);
	}
	updateDisplay();
	Plot.show();
	setBatchMode("show");
	saveAs("Tiff", output + plotName);

	if(frames > 1) {
		//Detect stimulation and calibration (if present) and sort kymograph
		if(manualStimCalFrames != "") {
			stimCalFrames = split(manualStimCalFrames, ",");
			if(stimCalFrames.length == 2) {
				stimulationFrame = (stimCalFrames[0]);
				calibrationFrame = (stimCalFrames[1]);
				if(stimulationFrame == 0) {
					calibrationOnly = true;
					print("Baseline and calibration only, at frame "+calibrationFrame);
				}
				print("Using manual stimulation and calibration frames: "+stimulationFrame+" and "+calibrationFrame);
			}
			else if (stimCalFrames.length == 1 && stimCalFrames[0] == 0) {
				baselineOnly = true;
				stimulationFrame = frames;
				calibrationFrame = -1;
				print("Baseline-only mode active");
			}
			else if (stimCalFrames.length == 1  && stimCalFrames[0] != 0) {
				print("WARNING!  Stimulation and calibration frames are set to manual, but the parameter field does not have two values!\nProceeding with automatic detection.");
				manualStimCalFrames = "";
			}
		}
		if(manualStimCalFrames == "") {
			diff_average_lifetime_trace = differentiateArray(average_lifetime_trace);
			diff2_average_lifetime_trace = differentiateArray(diff_average_lifetime_trace);
	
			if(debugMode) {
				Plot.create("2nd derivative trace", "time (s)", "dtau/dt | d2tau/dt2", timeArray, diff2_average_lifetime_trace);
				Plot.add("line", timeArray, diff_average_lifetime_trace);
				Plot.setStyle(0, "blue,none,2.0,Line");
				Plot.setStyle(1, "red,none,2.0,Line");
				Plot.setLegend("2nd derivative\t1st derivative", "bottom-to-top");
				Plot.setLimitsToFit();
				Plot.setLineWidth(1);
				Plot.show();
				setBatchMode("show");
			}
			Array.getStatistics(diff2_average_lifetime_trace, min_, max_, mean_, stdDev);
			maxima = Array.findMaxima(diff2_average_lifetime_trace, stimulationSensitivity * stdDev);	//Detect peaks in the second derivative with a certain minimum amount of change in the slope.
			if(debugMode) {
				print("\n2nd derivative stddev: "+stdDev);
				print("tolerance: "+stimulationSensitivity * stdDev);
				print(maxima.length+" peaks in the 2nd derivative found, at positions:");
				Array.print(maxima);
			}
			if(maxima.length == 0) {	//No stimulation and no calibration found
				stimulationFrame = 0;
				calibrationFrame = -1;
				print("WARNING: no stimulation frame and no calibration frame detected! Using the full trace for sorting the traces."); 
			}
			else if(maxima.length == 1) {	//No calibration found or no stimulation found
				stimulationFrame = maxima[0]-1;
				calibrationFrame = -1;
				print("WARNING: no calibration frame detected! Stimulation detected at "+stimulationFrame*frameInterval+" seconds (frame "+stimulationFrame+").");
			}
			else if(maxima.length >= 2 && maxima[0] < maxima[1]) { //Then stimulation causes a faster change than calibration: reverse them if this is not the case.
				stimulationFrame = maxima[0]-1;
				calibrationFrame = maxima[1]-1;
				print("Stimulation and calibration detected at frames "+stimulationFrame+1+" and "+calibrationFrame+1+" ("+d2s(stimulationFrame*frameInterval,1)+" and "+d2s(calibrationFrame*frameInterval,1)+" seconds).");
			}
			else if(maxima.length >= 2 && maxima[1] < maxima[0]) {
				stimulationFrame = maxima[1]-1;
				calibrationFrame = maxima[0]-1;
				print("Stimulation and calibration detected at frames "+stimulationFrame+1+" and "+calibrationFrame+1+" ("+d2s(stimulationFrame*frameInterval,1)+" and "+d2s(calibrationFrame*frameInterval,1)+" seconds).");
			}
		}

		selectWindow(kymograph);
		response = newArray(nr_cells);

		for (i = 0; i < nr_cells; i++) {
			showStatus("Measuring responses... "+i+"/"+nr_cells);
			//Sort on response
			if(calibrationFrame != -1) makeRectangle(stimulationFrame, i, calibrationFrame-1 -stimulationFrame, 1);
			//Alternative: sort on baseline:
	//		if(calibrationFrame != -1) makeRectangle(0, i, stimulationFrame, 1);
			else if (calibrationFrame == -1 && baselineOnly == false) makeRectangle(stimulationFrame, i, frames-stimulationFrame, 1);
			else if (baselineOnly == true) makeRectangle(0, i, frames, 1);
			profile = getProfile();
			Array.getStatistics(profile, min, max, mean, stdDev);
			response[i] = mean;
			run("Select None");
		}

		rank_response = Array.rankPositions(response);
	//	Array.show("positions", response, rank_response);
		kymographSorted = "kymograph sorted";
		selectWindow(kymograph);
		run("Duplicate...", "title=["+kymographSorted+"]");
		for (i = 0; i < nr_cells; i++) {
			selectWindow(kymograph);
			showStatus("Sorting responses... "+i+"/"+nr_cells);
			makeRectangle(0, rank_response[i], frames, 1);
			run("Copy");
			selectWindow(kymographSorted);
			makeRectangle(0, i, frames, 1);
			run("Paste");
		}
		run("Select None");

		rank_responseVector = "rank_responseVector";
		Ext.CLIJ2_pushArray(rank_responseVector, rank_response, nr_cells, 1, 1);
		//Save and rename back - TO DO: move all the saving to the end in the main function (?)

		//Rotate kymographs to horizontal position
		selectWindow(kymograph);
		run("Rotate 90 Degrees Right");
		saveAs("tiff", output + saveName + " (kymograph)");
		rename(kymograph);
		if(smoothRadiusTraces > 0) {
			selectWindow(kymograph_smoothed);
			run("Rotate 90 Degrees Right");
		}
		selectWindow(kymographSorted);
		run("Rotate 90 Degrees Right");
		run("Flip Horizontally");
		setBatchMode("show");
		saveAs("tiff", output + saveName + " (kymograph sorted)");
		rename(kymographSorted);
		
		selectWindow(kymograph);
		run("Flip Horizontally");
		setBatchMode("show");
		if(smoothRadiusTraces > 0) {
			selectWindow(kymograph_smoothed);
			run("Flip Horizontally");
			setBatchMode("show");
		}	
	}
	else {
		selectWindow(kymograph);
		run("Rotate 90 Degrees Right");
		run("Flip Horizontally");
		
		kymographSorted = "";
		rank_responseVector = "";
	}
	lifetimeTable = "Lifetime_Data";
	if(generateTables) {
		Table.rename("Results", lifetimeTable);
		//if(frames > 1) Table.setColumn("time (s)", timeArray, lifetimeTable);
		selectWindow(lifetimeTable);
		for (i = 0; i < nr_cells; i++) {
			Table.renameColumn("X"+i, "cell_"+IJ.pad(i+1,5));
		}
		if(frames > 1) {
			timeArray = Array.getSequence(frames);
			timeArray = multiplyArraywithScalar(timeArray, frameInterval);
			Table.setColumn("time (s)", timeArray);
		}
		Table.update;
		Table.save(output + saveName + "_lifetime.tsv");
	}

	//Save a text file with stimulationFrame and calibrationFrame
	txtfile = File.open(output + saveName + "_Stim_&_Cal_frames.txt");
	print(txtfile, stimulationFrame+"\t"+calibrationFrame);
	File.close(txtfile);

	return newArray(lifetimeTable, kymograph, kymograph_smoothed, kymographSorted, rank_responseVector);
}


//Measure the intensities of all the ROIs in the input timelapse image 
function measure_intensities(input_timelapse, nr_cells, generateTables) {

	if(generateTables) { dataTable = input_timelapse+"_data"; Table.create(dataTable); }
	run("Clear Results");
	selectImage(input_timelapse);
	getDimensions(width, height, channels, slices, frames);
	newImage(input_timelapse+"_kymograph", "32-bit black", nr_cells, frames, 1);
	
	selectImage(input_timelapse);
	roiManager("multi-measure measure_all append");
	dataAll = Table.getColumn("Mean", "Results");

	selectImage(input_timelapse+"_kymograph");
	for (f=0; f<frames; f++) {
		showProgress(f, frames);
		for(x=0; x<nr_cells; x++) setPixel(x, f, dataAll[f*nr_cells+x]);
		if(generateTables) {
			dataSingleFrame = Array.slice(dataAll, f*nr_cells, (f+1)*nr_cells);
			Table.setColumn("frame "+f+1, dataSingleFrame, dataTable);
		}
	}
	return input_timelapse+"_kymograph";
}


//Measure the intensities of all the labels in the input timelapse image on the GPU. Assumes that the labelmap is already present in the GPU memory.  [Currently not used]
function measure_intensities_on_GPU(labelmap, input_timelapse) {
	dataTable = input_timelapse+"_data";
	Table.create(dataTable);
	run("Clear Results");
	selectImage(input_timelapse);
	getDimensions(width, height, channels, slices, frames);
	for(i=1; i<=frames; i++) {
		showStatus("Measuring "+input_timelapse+", frame "+i+"/"+frames);
		Stack.setFrame(i);
		Ext.CLIJ2_pushCurrentSlice(input_timelapse);
		Ext.CLIJ2_statisticsOfLabelledPixels(input_timelapse, labelmap);
		mean = Table.getColumn("MEAN_INTENSITY", "Results");
		Table.setColumn("Mean_"+i, mean, dataTable);
		run("Clear Results");
		showProgress(i/frames);
	}
	Table.rename(dataTable, "Results");
}


function calculate_lifetime_and_intensity_TCSPC(image) {
	selectWindow(image);
	//Get rid of the extension in the name (if any)
	//name = substring(image, 0, lastIndexOf(image, "."));
	run("Duplicate...", "title=IntensityStack duplicate");
	run("32-bit");
	run("Split Channels");
	//Add component channels (here called 'intensity' and 'lifetime' channels - not really correct names)
	imageCalculator("Add create 32-bit stack", "C"+intensityChannel+"-IntensityStack","C"+lifetimeChannel+"-IntensityStack");
	rename("Added_Amplitudes");

	//Multiply amplitudes with lifetime components to get total intensities of the components
	selectWindow("C"+intensityChannel+"-IntensityStack");
	run("Multiply...", "value="+tau1+" stack");
	selectWindow("C"+lifetimeChannel+"-IntensityStack");
	run("Multiply...", "value="+tau2+" stack");
	imageCalculator("Add create 32-bit stack", "C"+intensityChannel+"-IntensityStack","C"+lifetimeChannel+"-IntensityStack");
	rename("Intensity");
	imageCalculator("Divide create 32-bit stack", "Intensity","Added_Amplitudes");
	rename(image+"_weighted_lifetime");
	close("C"+intensityChannel+"-IntensityStack");
	close("C"+lifetimeChannel+"-IntensityStack");

	//Set NaNs to 0 - Otherwise drift correction / movement (?) can cause artifacts.
	getDimensions(width, height, channels, slices, frames);
	if(frames>1) {
		for (f = 1; f <= frames; f++) {
	    	Stack.setFrame(f);
			changeValues(NaN, NaN, 0);
		}
	}
	else changeValues(NaN, NaN, 0);
	
	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_weighted_lifetime");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_weighted_lifetime");
}


function calculate_lifetime_and_intensity_TauContrast(image) {
	selectWindow(image);
	run("Duplicate...", "title=IntensityStack duplicate channels="+intensityChannel);
	run("32-bit");
	rename("Intensity");

	selectWindow(image);
	run("Duplicate...", "title=["+image+"_weighted_lifetime] duplicate channels="+lifetimeChannel);
	run("32-bit");
	//Scale TauContrast numbers to lifetime (-1 to 24 ns range equals 0-255 range)
	run("Multiply...", "value=0.097 stack");
	run("Subtract...", "value=1 stack");
	
//	Set pixels with value -1 (created after drift correction) to NaN - Creates artifacts!
//	setThreshold(-0.999, 24);
//	run("NaN Background", "stack");
//	resetThreshold();

	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_weighted_lifetime");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_weighted_lifetime");
}


function calculate_lifetime_and_intensity_FASTFLIM(image) {
	selectWindow(image);
	run("Duplicate...", "title=IntensityStack duplicate channels="+intensityChannel);
	run("32-bit");
	rename("Intensity");

	selectWindow(image);
	run("Duplicate...", "title=["+image+"_weighted_lifetime] duplicate channels="+lifetimeChannel);
	run("32-bit");
	//Scale TauContrast numbers to lifetime (-1 to 24 ns range equals 0-255 range)
	run("Multiply...", "value=0.001 stack");
	
//	Set pixels with value -1 (created after drift correction) to NaN - Creates artifacts!
//	setThreshold(-0.999, 24);
//	run("NaN Background", "stack");
//	resetThreshold();

	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_weighted_lifetime");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_weighted_lifetime");
}


function calculate_lifetime_and_intensity_FDFLIM(image) {
	openfli(reference,"reference");
	selectWindow(image);
//	name = substring(image, 0, lastIndexOf(image, "."));
//	getDimensions(w, h, channels, slices, frames);

	run("Duplicate...", "title=temp duplicate");
	frames = frames/phases;
	// loop over all time frames
	for (i = 0; i < (frames); i++) {
		selectWindow("temp");
		if (i<((frames)-1)){ //last frame
			run("Make Substack...", "delete slices=1-"+phases);
		}
		rename("temp_img");
		showProgress(i, frames);
		run("fdFLIM", "image1=[temp_img] boolphimod=false image2=reference tau_ref="+tau_ref+" freq="+freq);	
		close("temp_img");
		if (i==0) rename("Lifetimes_final");
		else run("Concatenate...", "  title=Lifetimes_final open image1=Lifetimes_final image2=Lifetimes");
	}
	close("reference");
//	close(image);

	//Separate intensity and lifetime
	selectWindow("Lifetimes_final");
	getDimensions(width, height, channels, slices, frames);
	if(tauPhiOrTauMod == "phase") lifetimeSlice = 1;
	else if(tauPhiOrTauMod == "modulation") lifetimeSlice = 2;
	if(frames>1) run("Duplicate...", "duplicate slices=3");
	else run("Duplicate...", "duplicate range=3-3");
	rename("Intensity");
	selectWindow("Lifetimes_final");
	if(frames>1) run("Duplicate...", "duplicate slices="+lifetimeSlice);
	else run("Duplicate...", "duplicate range="+lifetimeSlice+"-"+lifetimeSlice);
	rename(image+"_phase_lifetime");
	close("Lifetimes_final");

	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_phase_lifetime");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");

	return newArray(image+"_intensity", image+"_phase_lifetime");
}

function calculate_ratio_and_intensity(image) {
	run("Duplicate...", "title=IntensityStack duplicate");
	run("32-bit");
	run("Split Channels");
	imageCalculator("Add create 32-bit stack", "C1-IntensityStack","C2-IntensityStack");
	rename("Intensity");
	selectWindow("C1-IntensityStack");
	setThreshold(1,65536);
	run("NaN Background", "stack");
	selectWindow("C2-IntensityStack");
	setThreshold(1,65536);
	run("NaN Background", "stack");

	//Divide intensities to get the ratio
	imageCalculator("Divide create 32-bit stack", "C1-IntensityStack","C2-IntensityStack");
	rename(image+"_ratio");
	close("C1-IntensityStack");
	close("C2-IntensityStack");

	//rename and apply display settings
	selectWindow("Intensity");
	rename(image+"_intensity");
	run("Enhance Contrast", "saturated=0.35");	
	selectWindow(image+"_ratio");
	run(lut);
	setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_ratio");
}

function calculate_intensity(image) {
	selectWindow(image);
	run("Duplicate...", "title=["+image+"_intensity] duplicate");
	run("Enhance Contrast", "saturated=0.35");
	run("Duplicate...", "title=["+image+"_intensity_nobackground] duplicate");
	run("32-bit");
	setThreshold(1,65536);
	run("NaN Background", "stack");

	//rename and apply display settings
	//run(lut);
	//setMinAndMax(minLifetime, maxLifetime);
	setBatchMode("show");
	return newArray(image+"_intensity", image+"_intensity_nobackground");
}

// Open both sample and background from a .fli file and subtract background
function openfli(input,name) {
	if(FDFLIM_SPAD_camera == false) {
		run("Bio-Formats", "open=["+input+"] view=Hyperstack stack_order=XYCZT series_1");
		rename(name);
		run("Bio-Formats", "open=["+input+"] view=Hyperstack stack_order=XYCZT series_2");
		rename("BG");
		imageCalculator("Subtract 32-bit stack", name, "BG");
		close("BG");
		close(name);
		rename(name);
	}
	else if(FDFLIM_SPAD_camera == true) {
		if(endsWith(input, ".fli")) run("Bio-Formats", "open=["+input+"] view=Hyperstack stack_order=XYCZT");
		else open(input);
		run("Median...", "radius=1 stack");	//some immediate noise reduction
		rename(name);
	}
}


function segment_cells_manually(intensity_stack) {
	selectWindow(intensity_stack);
	getDimensions(width, height, channels, slices, frames);
	if(frames > 1 && slices == 1) {			//timelapse, single image
		if(CellposeStartFrame == -1 && CellposeEndFrame == -1) run("Z Project...", "projection=[Sum Slices] all");
		else if(CellposeStartFrame == -1) run("Z Project...", "stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		else if(CellposeEndFrame == -1) run("Z Project...", "start="+CellposeStartFrame+" projection=[Sum Slices] all");
		else run("Z Project...", "start="+CellposeStartFrame+" stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		rename("intensity");
	}
	else if(slices > 1 && frames > 1) {		//timelapse, multiple images
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		if(CellposeStartFrame == -1 && CellposeEndFrame == -1) run("Z Project...", "projection=[Sum Slices] all");
		else if(CellposeStartFrame == -1) run("Z Project...", "stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		else if(CellposeEndFrame == -1) run("Z Project...", "start="+CellposeStartFrame+" projection=[Sum Slices] all");
		else run("Z Project...", "start="+CellposeStartFrame+" stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		rename("intensity");
	}
	else if(slices > 1 && frames == 1) {	//no timelapse, multiple images
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		run("Duplicate...", "title=intensity duplicate");
	}										//single image
	else run("Duplicate...", "title=intensity duplicate");
	
	intensity_image = "intensity";
	Ext.CLIJ2_push(intensity_image);	//Push to GPU for intensity measurements
	
	if(equalize_contrast_cp == true) {
		//run("Square Root");
		run("Gamma...", "value=0.50");
		//run("Enhance Local Contrast (CLAHE)", "blocksize=64 histogram=256 maximum=3 mask=*None*");	//The fancy way, but not really better
		run("Enhance Contrast", "saturated=0.35");
	}
	selectWindow(intensity_image);
	setBatchMode("show");

	run("Duplicate...", "title=intensity_image_for_Cellpose duplicate");	//Keep "intensity" for RGB overlay later
	setBatchMode("show");

	//Manual segmentation
	run("Colors...", "foreground=white background=black selection=cyan");
	tool = IJ.getToolName();
	setTool("Freehand");
	roiManager("reset");
	roiManager("Show All with labels");
	waitForUser("Draw manual regions and add to the ROI Manager (shortcut: 't'), or click OK to use the whole image.");
	setTool(tool);
	if(roiManager("count") == 0) { run("Select All"); roiManager("Add"); }
	run("Select None");
	roiManager("Show None");
	run("ROI Manager to LabelMap(2D)");
	rename("labelmap_cells");
	roiManager("reset");
	run("glasbey on dark");

	labelmap_cells = "labelmap_cells";
	Ext.CLIJ2_push(labelmap_cells);
	Ext.CLIJ2_statisticsOfLabelledPixels(intensity_image, labelmap_cells);
	Ext.CLIJ2_release(labelmap_cells);

	Table.rename("Results", "Cell_statistics");
	Table.save(output + saveName + "_Cell_statistics.tsv");
	
	return newArray("labelmap_cells");
}


function segment_cells_no_nuclei(intensity_stack) {
	//Segment cells using Cellpose
	selectWindow(intensity_stack);
	getDimensions(width, height, channels, slices, frames);
	if(frames > 1 && slices == 1) {			//timelapse, single image
		if(CellposeStartFrame == -1 && CellposeEndFrame == -1) run("Z Project...", "projection=[Sum Slices] all");
		else if(CellposeStartFrame == -1) run("Z Project...", "stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		else if(CellposeEndFrame == -1) run("Z Project...", "start="+CellposeStartFrame+" projection=[Sum Slices] all");
		else run("Z Project...", "start="+CellposeStartFrame+" stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		rename("intensity");
	}
	else if(slices > 1 && frames > 1) {		//timelapse, multiple images
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		if(CellposeStartFrame == -1 && CellposeEndFrame == -1) run("Z Project...", "projection=[Sum Slices] all");
		else if(CellposeStartFrame == -1) run("Z Project...", "stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		else if(CellposeEndFrame == -1) run("Z Project...", "start="+CellposeStartFrame+" projection=[Sum Slices] all");
		else run("Z Project...", "start="+CellposeStartFrame+" stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		rename("intensity");
	}
	else if(slices > 1 && frames == 1) {	//no timelapse, multiple images
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		run("Duplicate...", "title=intensity duplicate");
	}										//single image
	else run("Duplicate...", "title=intensity duplicate");
	
	intensity_image = "intensity";
	Ext.CLIJ2_push(intensity_image);	//Push to GPU for intensity measurements
	
	if(equalize_contrast_cp == true) {
		//run("Square Root");
		run("Gamma...", "value=0.50");
		//run("Enhance Local Contrast (CLAHE)", "blocksize=64 histogram=256 maximum=3 mask=*None*");	//The fancy way, but not really better
		run("Enhance Contrast", "saturated=0.35");
	}
	selectWindow(intensity_image);
	setBatchMode("show");

	run("Duplicate...", "title=intensity_image_for_Cellpose duplicate");	//Keep "intensity" for RGB overlay later

	//Scale up for Cellpose and StarDist
	if(upSampleFactor > 1) {
		run("Scale...", "x="+upSampleFactor+" y="+upSampleFactor+" interpolation=Bicubic process create");
	}

//If 32-bit, first convert to 16-bit, with some contrast enhancement - better segmentation (sometimes)
	run("Set Measurements...", "mean redirect=None decimal=3");
	percentile_threshold(0.05);	//Threshold on the 5% lowest values
	List.setMeasurements("limit");
	mean = List.getValue("Mean");
	changeValues(NaN, NaN, mean);
	run("Enhance Contrast", "saturated="+saturatedPixels);	//TO DO: another way to reliably set the B&C settings
	run("Conversions...", "scale");
	run("16-bit");

	//Run Cellpose without nuclei channel
	selectWindow("intensity_image_for_Cellpose");
	setBatchMode("show");
	
	if(load_labelmap_boolean == false) {
		setBatchMode(false);	//Cellpose doesn't return an image in batch mode
		cellpose_success = false;
		selectWindow("intensity_image_for_Cellpose");
		while(cellpose_success == false) {
			if(newCellposeWrapper == false) {
				if(CellposeModel != "custom") run("Cellpose Advanced", "diameter="+CellposeDiameter+" cellproba_threshold="+CellposeProbability+" flow_threshold="+CellposeFlowThreshold+" anisotropy=1.0 diam_threshold=12.0 model="+CellposeModel+" nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
				else run("Cellpose Advanced (custom model)", "diameter="+CellposeDiameter+" cellproba_threshold="+CellposeProbability+" flow_threshold=0.7 anisotropy=1.0 diam_threshold=12.0 model_path=["+cellposeModelPath+"] model=["+cellposeModelPath+"] nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");
			}
			else if(newCellposeWrapper) {
				if(CellposeModel != "custom") run("Cellpose ...", "env_path="+env_path+" env_type="+env_type+" model="+CellposeModel+" model_path=path\\to\\own_cellpose_model diameter="+CellposeDiameter+" ch1=1 ch2=0 additional_flags=[--use_gpu, --flow_threshold="+CellposeFlowThreshold+", --cellprob_threshold="+CellposeProbability+"]");
				else run("Cellpose ...", "env_path="+env_path+" env_type="+env_type+" model=["+cellposeModelPath+"] model_path=["+cellposeModelPath+"] diameter="+CellposeDiameter+" ch1=1 ch2=0 additional_flags=[--use_gpu, --flow_threshold="+CellposeFlowThreshold+", --cellprob_threshold="+CellposeProbability+"]");
			}

			if(getTitle() != "intensity_image_for_Cellpose-cellpose") {
				print("Cellpose failed!\nCheck the Fiji console for more information. Waiting 60 seconds and then trying again...");
				wait(60000);
			}
			else cellpose_success = true;
		}
		setBatchMode("hide");

		rename("labelmap_cells_upsampled");

		//Downscale labelmap
		selectWindow("labelmap_cells_upsampled");
		if(upSampleFactor > 1) run("Scale...", "x="+1/upSampleFactor+" y="+1/upSampleFactor+" interpolation=None process create");
		rename("labelmap_cells");
		close("labelmap_cells_upsampled");
	}
	else {
		if(File.exists(labelmapPath)) open(labelmapPath);
		else exit("Error loading labelmap from disk: the file '"+labelmapPath+"' does not exist!");
		rename("labelmap_cells");
		setBatchMode("show");
	}

	//Measure circularity
	run("Analyze Regions", "circularity");
	if(Table.size("labelmap_cells-Morphometry") == 0) break;
	circularity_ = Table.getColumn("Circularity", "labelmap_cells-Morphometry");
	circularity_mask = createBinaryArrayWhereSmallerThanCutoff(circularity_, maxCircularity);
	circularity_mask_expanded = insertElementIntoArrayAtPosition(0, circularity_mask, 0);
	Ext.CLIJ2_pushArray(circularity_mask_GPU, circularity_mask_expanded, circularity_mask_expanded.length, 1, 1);

	labelmap_cells = "labelmap_cells";
	Ext.CLIJ2_push(labelmap_cells);
	close("labelmap_cells");
	
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_cells, labelmap_cells_sizeFiltered, minCellSize, 1e8);

	run("Clear Results");
	Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(intensity_image, labelmap_cells_sizeFiltered);
	Ext.CLIJ2_pushResultsTableColumn(mean_Intensity_vector, "MEAN_INTENSITY");
	Ext.CLIJ2_mask(mean_Intensity_vector, circularity_mask_GPU, mean_Intensity_vector_masked);
	Ext.CLIJ2_excludeLabelsWithValuesOutOfRange(mean_Intensity_vector_masked, labelmap_cells_sizeFiltered, labelmap_cells_sizeAndIntensityFiltered, minCellBrightness, 1e30);
	Ext.CLIJ2_release(labelmap_cells);
	Ext.CLIJ2_release(labelmap_cells_sizeFiltered);
	run("Clear Results");
	Ext.CLIJ2_statisticsOfLabelledPixels(intensity_image, labelmap_cells_sizeAndIntensityFiltered);
	Table.rename("Results", "Cell_statistics");
	Table.save(output + saveName + "_Cell_statistics.tsv");
	
	Ext.CLIJ2_release(circularity_mask_GPU);
	Ext.CLIJ2_release(mean_Intensity_vector);
	Ext.CLIJ2_release(mean_Intensity_vector_masked);
	Ext.CLIJ2_release(intensity_image);
	Ext.CLIJ2_pull(labelmap_cells_sizeAndIntensityFiltered);
	Ext.CLIJ2_release(labelmap_cells_sizeAndIntensityFiltered);
	
	rename("labelmap_cells");
	run("glasbey on dark");
	if(getValue("Max") <= 255) setMinAndMax(0, 255);
	else resetMinAndMax();
	setBatchMode("show");
	
	return newArray("labelmap_cells");
}


function segment_cells(intensity_stack, nuclei_stack) {
	//Segment cells using Cellpose and segment nuclei using Stardist
	selectWindow(intensity_stack);
	if(frames > 1 && slices == 1) {			//timelapse, single image
		if(CellposeStartFrame == -1 && CellposeEndFrame == -1) run("Z Project...", "projection=[Sum Slices] all");
		else if(CellposeStartFrame == -1) run("Z Project...", "stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		else if(CellposeEndFrame == -1) run("Z Project...", "start="+CellposeStartFrame+" projection=[Sum Slices] all");
		else run("Z Project...", "start="+CellposeStartFrame+" stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		rename("intensity");
	}
	else if(slices > 1 && frames > 1) {		//timelapse, multiple images
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		if(CellposeStartFrame == -1 && CellposeEndFrame == -1) run("Z Project...", "projection=[Sum Slices] all");
		else if(CellposeStartFrame == -1) run("Z Project...", "stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		else if(CellposeEndFrame == -1) run("Z Project...", "start="+CellposeStartFrame+" projection=[Sum Slices] all");
		else run("Z Project...", "start="+CellposeStartFrame+" stop="+CellposeEndFrame+" projection=[Sum Slices] all");
		rename("intensity");
	}
	else if(slices > 1 && frames == 1) {	//no timelapse, multiple images
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		run("Duplicate...", "title=intensity duplicate");
	}										//single image
	else run("Duplicate...", "title=intensity duplicate");

	intensity_image= "intensity";
	selectWindow(intensity_image);
	run("Grays");
	run("Enhance Contrast", "saturated=0.35");
	Ext.CLIJ2_push(intensity_image);	//Push to GPIU for intensity measurements
	
	if(equalize_contrast_cp == true) {
		//run("Square Root");
		run("Gamma...", "value=0.50");
		//run("Enhance Local Contrast (CLAHE)", "blocksize=64 histogram=256 maximum=3 mask=*None*");	//The fancy way, but not really better
		run("Enhance Contrast", "saturated=0.35");
	}
	setBatchMode("show");
	run("Duplicate...", "title=intensity_image_for_Cellpose duplicate");	//Keep "intensity" for RGB overlay later

	selectWindow(nuclei_stack);
	if(frames > 1) {
		if(CellposeEndFrame == -1) run("Z Project...", "projection=[Sum Slices]");
		else run("Z Project...", "stop="+CellposeEndFrame+" projection=[Sum Slices]");
		rename("nuclei");
	}
	else run("Duplicate...", "title=nuclei");

	//Correct for bidirectional phase mismatch - currently only on the projection, for nuclei segmentation
	if(correctBidirect_boolean == true) nuclei_corrected = correct_bidirectional_phase("nuclei");
	else nuclei_corrected = "nuclei";

	run("Merge Channels...", "c1="+intensity_image+" c2="+nuclei_corrected+" create keep");
	rename("merged");
	Stack.setChannel(2);
	run("biop-Azure");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(1);
	run("Grays");
	run("Enhance Contrast", "saturated=0.35");
	setBatchMode("show");

	//Scale up for Cellpose and StarDist
	if(upSampleFactor > 1) {
		run("Scale...", "x="+upSampleFactor+" y="+upSampleFactor+" interpolation=Bicubic process create");
	}
	rename("intensity_image_for_Cellpose");

	//Run Cellpose with nuclei channel
	selectWindow("intensity_image_for_Cellpose");
	setBatchMode("show");
	setBatchMode("exit and display");	//Cellpose doesn't return an image in batch mode
	selectWindow("intensity_image_for_Cellpose");


//TO DO: If 32-bit, first convert to 16-bit, with some contrast enhancement - better segmentation!
	run("Enhance Contrast", "saturated="+saturatedPixels);	//TO DO: another way to reliably set the B&C settings
	run("Conversions...", "scale");
	run("16-bit");
	Stack.setChannel(2);
	run("biop-Azure");
	run("Enhance Contrast", "saturated=0.35");
	Stack.setChannel(1);
	run("Grays");
	run("Enhance Contrast", "saturated=0.35");

//!!!! TO DO: CHANGE THIS FOR NEW CELLPOSE VERSIONS AND INCLUDE NUCLEI
	run("Cellpose Advanced", "diameter="+CellposeDiameter+" cellproba_threshold="+CellposeProbability+" flow_threshold="+CellposeFlowThreshold+" anisotropy=1.0 diam_threshold=12.0 model="+CellposeModel+" nuclei_channel=0 cyto_channel=1 dimensionmode=2D stitch_threshold=-1.0 omni=false cluster=false additional_flags=");

	setBatchMode(true);
	if(getTitle() != "intensity_image_for_Cellpose-cellpose") exit("Cellpose failed!\nExiting macro.");

	rename("labelmap_cells_upsampled");
	run("glasbey on dark");
	resetMinAndMax();

	//Run StarDist
	selectWindow("intensity_image_for_Cellpose");
	Stack.setChannel(2);
	run("Duplicate...", "title=nuclei_upsampled");
	run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'nuclei_upsampled', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'2', 'percentileTop':'98', 'probThresh':'"+StarDistProbThreshold+"', 'nmsThresh':'0.3', 'outputType':'ROI Manager', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Stack', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
	run("ROI Manager to LabelMap(2D)");
	roiManager("reset");
	setBatchMode("show");
	run("glasbey on dark");
	resetMinAndMax();
	rename("labelmap_nuclei_upsampled");

	//Downscale labelmaps
	selectWindow("labelmap_cells_upsampled");
	if(upSampleFactor > 1) run("Scale...", "x="+1/upSampleFactor+" y="+1/upSampleFactor+" interpolation=None process create");
	rename("labelmap_cells");
	setBatchMode("show");
	close("labelmap_cells_upsampled");
	selectWindow("labelmap_nuclei_upsampled");
	if(upSampleFactor > 1) run("Scale...", "x="+1/upSampleFactor+" y="+1/upSampleFactor+" interpolation=None process create");
	rename("labelmap_nuclei");
	setBatchMode("show");
	close("labelmap_nuclei_upsampled");
	close("merged");

	labelmap_cells = "labelmap_cells";
	Ext.CLIJ2_push(labelmap_cells);
	close("labelmap_cells");
	
	Ext.CLIJ2_excludeLabelsOutsideSizeRange(labelmap_cells, labelmap_cells_sizeFiltered, minCellSize, 1e8);

	run("Clear Results");
	Ext.CLIJ2_statisticsOfBackgroundAndLabelledPixels(intensity_image, labelmap_cells_sizeFiltered);
	Ext.CLIJ2_pushResultsTableColumn(mean_Intensity_vector, "MEAN_INTENSITY");
	Ext.CLIJ2_excludeLabelsWithValuesOutOfRange(mean_Intensity_vector, labelmap_cells_sizeFiltered, labelmap_cells_sizeAndIntensityFiltered, minCellBrightness, 1e30);
	Ext.CLIJ2_release(labelmap_cells);
	Ext.CLIJ2_release(labelmap_cells_sizeFiltered);
	run("Clear Results");
	Ext.CLIJ2_statisticsOfLabelledPixels(intensity_image, labelmap_cells_sizeAndIntensityFiltered);
	Table.rename("Results", "Cell_statistics");
	Table.save(output + saveName + "_Cell_statistics.tsv");
	
	Ext.CLIJ2_release(intensity_image);
	Ext.CLIJ2_pull(labelmap_cells_sizeAndIntensityFiltered);
	Ext.CLIJ2_release(labelmap_cells_sizeAndIntensityFiltered);
	rename("labelmap_cells");
	run("glasbey on dark");
	resetMinAndMax();
	setBatchMode("show");

	return newArray("labelmap_cells", "labelmap_nuclei");
}


function correct_drift(image) {
	selectWindow(image);
	getDimensions(width, height, channels, slices, frames);
/*	
	//Change z and t dimensions, if applicable
	Stack.getPosition(channel, slice, frame);
	getDimensions(width, height, channels, slices, frames);
	if(frames == 1 && slices == 1) exit("No time- or z-dimension found.");
	else if(frames == 1 && slices >= 1) {
		print("Swapping z and t dimensions.");
		run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
		getDimensions(width, height, channels, slices, frames);
	}
*/
	//Prepare intensity image used for registration
	if (microscope == "Confocal TCSPC / TauSeparation" || microscope == "Ratio Imaging") {
		selectWindow(image);
		run("Duplicate...", "title=C1_Intensity duplicate channels="+parseInt(intensityChannel));
		selectWindow(image);
		run("Duplicate...", "title=C2_Intensity duplicate channels="+parseInt(intensityChannel+1));
		imageCalculator("Add create 32-bit stack", "C1_Intensity","C2_Intensity");
		rename("Intensity_for_registration");
		close("C1_Intensity");
		close("C2_Intensity");
	}
	else run("Duplicate...", "title=Intensity_for_registration duplicate channels="+parseInt(intensityChannel));

	//Apply variance filter when edge detection is enabled
	if(edge_detect) {
		run("Duplicate...", "title=[For_FFT_"+image+"] duplicate");
		run("Variance...", "radius=2 stack");
	}
	
	//Change image size if required (has to be a power of two for the FFT)
	if(!edge_detect) run("Duplicate...", "title=[For_FFT_"+image+"] duplicate");
	dim = pad_image_edges("For_FFT_"+image, "Top-Left");

	reg_image = getTitle();	//image used for registration
	
	Xpos_max = newArray(frames);
	Ypos_max = newArray(frames);
	
	concat_string = "";
	if(registration_against == "First frame") run("Duplicate...", "title=frameA duplicate range=1-1");
	if(registration_against == "Last frame") run("Duplicate...", "title=frameA duplicate range="+frames+"-"+frames);
	for(f=1;f<frames;f++) {
		showProgress(f, frames);
		if(registration_against == "Previous frame") {
			selectWindow(reg_image);
			run("Duplicate...", "title=frameA duplicate range="+f+"-"+f);
		}
		selectWindow(reg_image);
		run("Duplicate...", "title=frameB duplicate range="+f+1+"-"+f+1);
		//Calculate cross-correlation
		run("FD Math...", "image1=frameA operation=Correlate image2=frameB result=CC"+f+" do");
		concat_string+="image"+f+"=CC"+f+" ";
		if(registration_against == "Previous frame") close("FrameA");
		close("FrameB");
	}
	close("FrameA");
	if(registration_against != "Previous frame") close("frameB");
	close("For_FFT_"+image);
	
	//concatenate cross correlation images
	crossCorrelationImage = "Cross Correlation";
	run("Concatenate...", "  title=["+crossCorrelationImage+"]" + concat_string);
	run("Gaussian Blur...", "sigma="+sigma_CC+" stack");
	
	run("Clear Results");
	run("CLIJ2 Macro Extensions", "cl_device=");
	Ext.CLIJ2_clear();
	Ext.CLIJ2_push(crossCorrelationImage);
	
	//determine max X positions for all frames
	Ext.CLIJ2_maximumYProjection(crossCorrelationImage, max_Y_projection);
	Ext.CLIJ2_transposeXZ(max_Y_projection, max_Y_TransPoseXZ);
	Ext.CLIJ2_argMaximumZProjection(max_Y_TransPoseXZ, max_XY_projection, arg_max_X);
	//determine max Y positions for all frames
	Ext.CLIJ2_maximumXProjection(crossCorrelationImage, max_X_projection);
	Ext.CLIJ2_transposeYZ(max_X_projection, max_X_TransPoseYZ);
	Ext.CLIJ2_argMaximumZProjection(max_X_TransPoseYZ, max_XY_projection, arg_max_Y);
	Ext.CLIJ2_transposeXY(arg_max_Y, arg_max_Y_transposed);
	//Put X and Y coordinates in the Results table and then in arrays
	Ext.CLIJ2_combineHorizontally(arg_max_X, arg_max_Y_transposed, max_XY_positions);
	Ext.CLIJ2_pullToResultsTable(max_XY_positions);
	Xpos_max = Table.getColumn("X0", "Results");
	Ypos_max = Table.getColumn("X1", "Results");
	
	//Subpixel (center of mass) determination
	if(subpixel) {
		run("Clear Results");
		if(show_CC) {
			crossCorrelationCropImage = "Normalized Cross Correlation Crop";
			Ext.CLIJ2_create3D(crossCorrelationCropImage, cropSizeCC, cropSizeCC, frames-1, 32);
		}
		for(f=1; f<frames; f++) {
			Ext.CLIJ2_copySlice(crossCorrelationImage, crossCorrelationFrame, f-1);
			Ext.CLIJ2_crop2D(crossCorrelationFrame, crossCorrelationFrameCrop, Xpos_max[f-1]-floor(cropSizeCC/2), Ypos_max[f-1]-floor(cropSizeCC/2), cropSizeCC, cropSizeCC);
			Ext.CLIJx_normalize(crossCorrelationFrameCrop, crossCorrelationFrameCropNormalized);
			Ext.CLIJ2_centerOfMass(crossCorrelationFrameCropNormalized);
			Xpos_max[f-1] = Xpos_max[f-1] + getResult("MassX") - floor(cropSizeCC/2);
			Ypos_max[f-1] = Ypos_max[f-1] + getResult("MassY") - floor(cropSizeCC/2);
			if(show_CC) Ext.CLIJ2_copySlice(crossCorrelationFrameCropNormalized, crossCorrelationCropImage, f-1);
		}
		if(show_CC) Ext.CLIJ2_pull(crossCorrelationCropImage);
		setBatchMode("show");
	}
	Ext.CLIJ2_clear();
	
	//Calculate translations
	translate_x = newArray(frames);
	translate_y = newArray(frames);
	for(f=1; f<frames; f++) {
		selectWindow(crossCorrelationImage);
		if(registration_against == "Previous frame") {
			if(f==1) {
				translate_x[f-1] = Xpos_max[f-1]-dim/2;
				translate_y[f-1] = Ypos_max[f-1]-dim/2;
			}
			if(f>=2) {
				translate_x[f-1] = translate_x[f-2] + Xpos_max[f-1]-dim/2;
				translate_y[f-1] = translate_y[f-2] + Ypos_max[f-1]-dim/2;
			}
		}
		else {
			translate_x[f-1] = Xpos_max[f-1]-dim/2;
			translate_y[f-1] = Ypos_max[f-1]-dim/2;
		}
		if(show_CC) {
			setSlice(f);
			makePoint(Xpos_max[f-1], Ypos_max[f-1], "medium red dot add");
		}
	}
	
	if(show_CC == false) close(crossCorrelationImage);
	
	if(edge_detect || channels>1) close(reg_image);
	close("Intensity_for_registration");
	
	selectWindow(image);
	run("Select None");
	run("Duplicate...", "title=["+image+"_driftcorr] duplicate");
	//print("Translation (pixels):");
	for(c=1; c<=channels; c++) {
		Stack.setChannel(c);
		for(f=1; f<frames; f++) {
			showProgress(f + (c-1)*frames, frames*channels);
			showStatus("Translating channel "+c+", frame "+f+"/"+frames);
			selectWindow(image+"_driftcorr");
			if(registration_against != "Last frame") Stack.setFrame(f+1);
			else Stack.setFrame(f);
			//print("channel "+c+", frame "+f+": "+translate_x[f-1]+", "+translate_y[f-1]);
			if(subpixel) run("Translate...", "x="+translate_x[f-1]+" y="+translate_y[f-1]+" interpolation=Bicubic slice");
			else run("Translate...", "x="+translate_x[f-1]+" y="+translate_y[f-1]+" interpolation=None slice");
		}
	}
	run("Clear Results");
	selectWindow(image);
	rename(image+"_NOT_drift_corrected");
	selectWindow(image+"_driftcorr");
	rename(image);
	setBatchMode("show");
	return image;
}


//expand the canvas to the nearest power of 2 and return the dimension
function pad_image_edges(image, location) {
	selectWindow(image);
	width = getWidth();
	height = getHeight();
	w=1;
	h=1;
	while(width>pow(2,w)) w++;
	while(height>pow(2,h)) h++;
	dim = pow(2,maxOf(w,h));
	//print("dimension of FFT: "+dim);
	run("Canvas Size...", "width="+dim+" height="+dim+" position="+location+" zero");
	return dim;
}


function correct_bidirectional_phase(image) {
//Adjust the phase of bidirectional confocal images

	image = getTitle();
	getDimensions(orgWidth, orgHeight, channels, slices, frames);
	getVoxelSize(pixelWidth, pixelHeight, depth, unit);
	run("Properties...", "pixel_width=1 pixel_height=1 voxel_depth=1.0000");	//Remove pixel calibration (if any)

	//Change image size if required (has to be a power of two for the FFT)

	if(orgWidth > 4096 || orgHeight > 4096) {
		print("The image is very large ("+orgWidth+" x "+orgHeight+"). Bidirectional correction may e slow.");
	}
	original = getTitle();
	run("Duplicate...", "title=[For_FFT] duplicate");
	dim = pad_image_edges("For_FFT", "Top-Left");
	getDimensions(width, height, channels, slices, frames);
	image = getTitle();

	run("Reslice [/]...", "output=1 start=Left avoid");
	run("Deinterleave", "how=2 keep");

	//Determine pixel shift with cross-correlation
	selectWindow("Reslice of "+image+" #1");
	run("Reslice [/]...", "output=1 start=Top rotate avoid");
	run("Canvas Size...", "width="+width+" height="+width+" position=Center zero");
	rename("Odd");
	getDimensions(_width, _height, _channels, _slices, _frames);
	if(_slices > 1) {
		run("Z Project...", "projection=[Sum Slices]");
		close("Odd");
	}
	rename("Odd");
	
	selectWindow("Reslice of "+image+" #2");
	run("Reslice [/]...", "output=1 start=Top rotate avoid");
	run("Canvas Size...", "width="+width+" height="+width+" position=Center zero");
	rename("Even");
	getDimensions(_width, _height, _channels, _slices, _frames);
	if(_slices > 1) {
		run("Z Project...", "projection=[Sum Slices]");
		close("Even");
	}
	rename("Even");
	
	image1 = "Odd";
	image2 = "Even";
	run("FD Math...", "image1=["+image1+"] operation=Correlate image2=["+image2+"] result=CC do");
	getStatistics(area, mean, min, max, std, histogram);
	run("Subtract...", "value="+min);
	run("Divide...", "value="+max-min);
	resetMinAndMax();
	
	//fit profile to a rectangle
	makeRectangle(width/2-16, width/2-16, 32, 32);
	run("Crop");
	run("Rotate 90 Degrees Right");
	run("Select All");
	profile = getProfile();
	x_points = Array.getSequence(32);
	Fit.doFit("Gaussian", x_points, profile);
	//Fit.plot;
	y0 = Fit.p(2);
	shift = (y0 - 15);
	print("Correcting bidirectional pixel shift: "+shift+", Rsq = "+Fit.rSquared);
	
	//Correct and recombine
	selectWindow("Reslice of "+image+" #1");
	//run("Translate...", "x="+shift+" y=0 interpolation=Bicubic stack");	//Doesn't allow subpixel translation on images 1 pixels high
	run("TransformJ Translate", "x-distance="+shift/2+" y-distance=0.0 z-distance=0.0 interpolation=[Cubic Convolution] background=0.0");
	selectWindow("Reslice of "+image+" #2");
	run("TransformJ Translate", "x-distance="+-shift/2+" y-distance=0.0 z-distance=0.0 interpolation=[Cubic Convolution] background=0.0");
	run("Interleave", "stack_1=[Reslice of "+image+" #1 translated] stack_2=[Reslice of "+image+" #2 translated]");
	run("Reslice [/]...", "output=1.000 start=Top rotate avoid");
	
	close("For_FFT");
	close("Reslice of "+image);
	close("Reslice of "+image+" #1");
	close("Reslice of "+image+" #1 translated");
	close("Reslice of "+image+" #2");
	close("Reslice of "+image+" #2 translated");
	close("Combined Stacks");
	close("Even");
	close("Odd");
	close("CC");

	corrected_image = image + "_phase_corrected";
	rename(corrected_image);
	setVoxelSize(pixelWidth, pixelHeight, depth, unit);
	makeRectangle(0, 0, orgWidth, orgHeight);
	run("Crop");
	run("Select None");
	
	getDimensions(width, height, channels, slices, frames);
	if(slices>1 && frames==1) run("Re-order Hyperstack ...", "channels=[Channels (c)] slices=[Frames (t)] frames=[Slices (z)]");
	run("Grays");
	return corrected_image;
}


//Overlay intensity image with the smoothed lifetime image, and add a calibration bar 
function overlay_intensity(intensity_image, lifetime_stack, saveName, smoothRadiusOverlayXY, smoothRadiusOverlayTime) {
	selectWindow(lifetime_stack);
	getDimensions(width, height, channels, slices, frames);
	getPixelSize(unit, pixelWidth, pixelHeight);
	if(smoothRadiusOverlayXY > 0 || smoothRadiusOverlayTime > 0) {
		lifetime_stack_for_smoothing = lifetime_stack+"_noNaNs";
		run("Duplicate...", "title=["+lifetime_stack_for_smoothing+"] duplicate");
		//Create NaNs from zeros and then 'pad' them before smoothing on the GPU
		getDimensions(width, height, channels, slices, frames);
		if(frames>1) {
			for (f = 1; f <= frames; f++) {
		    	Stack.setFrame(f);
				changeValues(0, 0, NaN);
			}
		}
		else changeValues(0, 0, NaN);
		run("Remove NaNs...", "radius="+smoothRadiusOverlayXY+" stack");
		Ext.CLIJ2_push(lifetime_stack_for_smoothing);
		close(lifetime_stack_for_smoothing);
		if(frames>1) Ext.CLIJ2_mean3DSphere(lifetime_stack_for_smoothing, lifetime_stack_filtered, smoothRadiusOverlayXY, smoothRadiusOverlayXY, smoothRadiusOverlayTime);
		else Ext.CLIJ2_mean2DSphere(lifetime_stack_for_smoothing, lifetime_stack_filtered, smoothRadiusOverlayXY, smoothRadiusOverlayXY);
		Ext.CLIJ2_pull(lifetime_stack_filtered);
		run(lut);
		setMinAndMax(minLifetime, maxLifetime);
		Ext.CLIJ2_release(lifetime_stack_for_smoothing);
		Ext.CLIJ2_release(lifetime_stack_filtered);
	}
	else run("Duplicate...", "title=splits duplicate");

	run("Calibration Bar...", "location=["+calibrationBarPosition+"] fill=Black label=White number=5 decimal=1 font=12 zoom=1 overlay");
	Overlay.copy;
	rename("splits");
//	run("Duplicate...", "title=splits duplicate");
//setBatchMode("show");
//setBatchMode(false);
	run("RGB Color");
	run("Split Channels");
	selectWindow(intensity_image);
	if(nucleiChannel != -1) {
		run("Duplicate...", "duplicate channels=1");	//Only the cells, not the nuclei
		intensity_image = getTitle();
	}

	run("Enhance Contrast", "saturated="+RGB_brightness);	//TO DO: another way to reliably set the B&C settings
	run("Conversions...", "scale");
	run("16-bit");
//	run("Apply LUT", "stack");
	showStatus("Generating overlay...");
	imageCalculator("Multiply 32-bit stack", "splits (red)", intensity_image);
	rename("Red");
	setMinAndMax(0, 65536*256);
	imageCalculator("Multiply 32-bit stack", "splits (green)", intensity_image);
	rename("Green");
	setMinAndMax(0, 65536*256);
	imageCalculator("Multiply 32-bit stack", "splits (blue)", intensity_image);
	rename("Blue");
	setMinAndMax(0, 65536*256);

	run("Merge Channels...", "c1=Red c2=Green c3=Blue");
	rename(lifetime_stack + " (RGB overlay)");
	Overlay.paste;
	rename("RGB_overlay");
	if(frames>1) Stack.setSlice(1);
	setBatchMode("show");
	Stack.setXUnit("um");
	run("Properties...", "channels=1 slices=1 frames="+frames+" pixel_width="+pixelWidth+" pixel_height="+pixelHeight+" voxel_depth=1.0000 frame="+frameInterval);

	close("splits (red)");
	close("splits (green)");
	close("splits (blue)");

	return "RGB_overlay";
}


function close_windows(name) {
	windowList = getList("window.titles");
	for(i=0 ; i<windowList.length ; i++) {
		if(matches(windowList[i],".*"+name+".*")) {
			selectWindow(windowList[i]);
			run("Close");
		}
	}
}

//Lower threshold an image at a certain percentile
function percentile_threshold(percentile) {
	getRawStatistics(nPixels, mean, min, max, std, histogram);
	total = 0;
	bin=0;
	while (total < nPixels*percentile) {
		total += histogram[bin];
		bin++;
	} 
	setThreshold(bin-1, max);
}


//Rename table headers (remove "Mean([cell_xxx])")
function renameTableHeaders(table) {
	headings = Table.headings(table);
	headers = split(headings, "\t");
	for(i=0;i<headers.length;i++) {
		newHeader = substring(headers[i+1],5,lengthOf(headers[i+1])-1);
		Table.renameColumn(headers[i+1], newHeader);
	}
	Table.update;
}


//Adds two arrays of equal length element-wise
function addArrays(array1, array2) {
	added_array=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		added_array[a]=array1[a] + array2[a];
	}
	return added_array;
}


//Subtracts two arrays of equal length element-wise
function subtractArrays(array1, array2) {
	subtracted_array=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		subtracted_array[a]=array1[a] - array2[a];
	}
	return subtracted_array;
}


//Returns the differentiated array. The first element is set to zero
function differentiateArray(array) {
	diffArray = newArray(array.length-1);
	for (i = 1; i < array.length-1; i++) {
		diffArray[i] = array[i] - array[i-1];
	}
	diffArray[0] = 0;
	return diffArray;
}


//Adds a scalar to all elements of an array
function addScalarToArray(array, scalar) {
	added_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		added_array[a] = (array[a]) + (scalar);
	}
	return added_array;
}


//Divides the elements of two arrays and returns the new array
function divideArrays(array1, array2) {
	divArray=newArray(lengthOf(array1));
	for (a=0; a<lengthOf(array1); a++) {
		divArray[a]=array1[a]/array2[a];
	}
	return divArray;
}


//Multiplies all elements of an array with a scalar
function multiplyArraywithScalar(array, scalar) {
	multiplied_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		multiplied_array[a]= (array[a]) * (scalar);
	}
	return multiplied_array;
}


//Divides all elements of an array by a scalar
function divideArraybyScalar(array, scalar) {
	divided_array=newArray(lengthOf(array));
	for (a=0; a<lengthOf(array); a++) {
		divided_array[a]=array[a]/scalar;
	}
	return divided_array;
}


//Returns the indices at which a value occurs within an array
function indexOfArray(array, value) {
	count=0;
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			count++;
		}
	}
	if (count>0) {
		indices=newArray(count);
		count=0;
		for (a=0; a<lengthOf(array); a++) {
			if (array[a]==value) {
				indices[count]=a;
				count++;
			}
		}
		return indices;
	}
}


//Returns the first index at which a value occurs in an array
function firstIndexOfArray(array, value) {
	for (a=0; a<lengthOf(array); a++) {
		if (array[a]==value) {
			break;
		}
	}
	return a;
}


//Returns the index of the maximum of an array
function maxIndexOfArray(array) {
	Array.getStatistics(array, min, max, mean, stdDev);
	index = indexOfArray(array, max);
	return index[0];
}


//Converts an array into a string, elements separated by 'separator'
function arrayToString(array, separator) {
	outputString = "";
	for (i = 0; i < array.length; i++) {
		outputString += toString(array[i]) + separator;
	}
	return substring(outputString, 0, outputString.length - separator.length);
}


//Returns the average of all elements of an arrays, neglecting NaNs
function averageArray(array) {
	sum=0;
	nans=0;
	for (a=0; a<lengthOf(array); a++) {
		if(!isNaN(array[a])) sum=sum+array[a];
		else nans+=1;
	}
	return sum/(array.length-nans);
}


function occursInArray(array, value) {
	for(i=0; i<array.length; i++) {
		if(array[i] == value) return true;
	}
	return false;
}


//Returns a binary array (0-1) where all elements smaller than cutoff are set to 1 and the rest to 0
function createBinaryArrayWhereSmallerThanCutoff(array, cutoff) {
	for(i=0; i<array.length; i++) {
		if(array[i] < cutoff) array[i] = 1;
		else array[i] = 0;
	}
	return array;
}


//Insert an element with value at a certain position
function insertElementIntoArrayAtPosition(value, array, position) {
	if (position<lengthOf(array)) {
		Array.rotate(array, -position);
		Array.reverse(array);
		array[array.length]=value;
		Array.reverse(array);
		Array.rotate(array, position);
	}
	else array[array.length]=value;
	return array;
}


//Returns an array where the elements are randomly shuffled
function shuffle_array(array) {
	nrElements = array.length;
	randomArray = newArray(nrElements);
	for (i = 0; i < nrElements; i++) {
		randomArray[i] = random;
	}
	rankArray = Array.rankPositions(randomArray);
	shuffledArray = newArray(nrElements);
	for (i = 0; i < nrElements; i++) {
		shuffledArray[i] = array[rankArray[i]];
	}
	return shuffledArray;
}


//Appends the value to the array
function appendToArray(value, array) {
	temparray=newArray(lengthOf(array)+1);
	for (i=0; i<lengthOf(array); i++) {
		temparray[i]=array[i];
	}
	temparray[lengthOf(temparray)-1]=value;
	array=temparray;
	return array;
}


//Reverses a table
function reverse_table(inputTable){
	headings = split(Table.headings(inputTable), "\t");
    for (col=0; col<headings.length; col++) {
    	if(col==0 && headings[col]==" ") continue;
    	else {
			col_values = Table.getColumn(headings[col], inputTable);
			Array.reverse(col_values);
			Table.setColumn(headings[col], col_values, inputTable);
    	}
    }
	Table.update(inputTable);
}

Table.update;


//Get the persistent value of the script parameter 'param' in class. N.B. This returns 'null' when the parameter is set to the default value!
function getPref(class, param) {
	return eval("js",
		"var ctx = Packages.ij.IJ.runPlugIn(\"org.scijava.Context\", \"\");" +
		"var ps = ctx.service(Packages.org.scijava.prefs.PrefService.class);" +
		"var " + param + " = ps.get(" + class + ".class, \"" + param + "\", \"<null>\");" +
		param + ";"
	);
}


//Gets a persistent Script Parameter value. N.B. This returns 'null' when the parameter is set to the default value!
function getScriptParameterValue(string_parameter) {
	return eval("js",
		"var ctx = Packages.ij.IJ.runPlugIn(\"org.scijava.Context\", \"\");" +
		"var ps = ctx.service(Packages.org.scijava.prefs.PrefService.class);" +
		"var " + string_parameter + " = ps.get(org.scijava.script.ScriptModule.class, \"" + string_parameter + "\", \"<null>\");" +
		string_parameter + ";"
	);
}


//Sets a persistent Script Parameter value. N.B. This returns 'null' when the parameter is set to the default value!
function setScriptParameterValue(string_parameter, string_value) {
	eval("js",
		"var ctx = Packages.ij.IJ.runPlugIn(\"org.scijava.Context\", \"\");" +
		"var ps = ctx.service(Packages.org.scijava.prefs.PrefService.class);" +
		"var " + string_parameter + " = ps.put(org.scijava.script.ScriptModule.class, \"" + string_parameter + "\", \"" + string_value + "\");"
	);
}

