# FAST-HIPPOS: FLIM Analysis of Single-cell Traces for Hit Identification of Phenotypes in Pooled Optical Screening

FAST-HIPPOS is a collection of Fiji scripts to analyze and visualize multi-cell time-lapse experiments, detect hit cells based on user-set criteria and output their stage coordinates for e.g. photoactivation.
For non-screening applications it can function as a valuable tool for single-cell trace analysis, visualization and inspection.

## Simplified workflow
1. *[Loading input images](https://github.com/Jalink-lab/dynamic-pooled-screening/#1-input-images)*
2. *[Pre-processing](https://github.com/Jalink-lab/dynamic-pooled-screening/#2-pre-processing)*
3. *[Cell segmentation](https://github.com/Jalink-lab/dynamic-pooled-screening/#3-cell-segmentation)*
4. *[Measuring single-cell (FLIM/intensity) traces](https://github.com/Jalink-lab/dynamic-pooled-screening/#4-measuring-single-cell-flimintensity-traces)*
5. *[Visualization](https://github.com/Jalink-lab/dynamic-pooled-screening/#5-visualization)*
6. *[Hit selection](https://github.com/Jalink-lab/dynamic-pooled-screening/#6-screening-hit-selection)*
7. *Output files*

<hr>

## Installation & requirements
Required activated [Fiji Update Sites](https://imagej.net/update-sites/following):
- Fiji Update sites:
  - ImageScience
  - ...

In order to run the macro you need the following scripts:
- several scripts/macros -? find out which:
  - inspect and select traces (macro tool, now in startup macros)
  - helper scripts
  - stitching
  - drift
  - bidir correction


## 1. Input images
Input images are multi-channel `.tif` files, or any proprietary microscopy format that is [supported by Bio-Formats](https://bio-formats.readthedocs.io/en/v8.2.0/supported-formats.html), including files containing multiple images.
FLIM is supported in several ways (mostly for Leica images):
  - Confocal TCSPC: lifetime component images, fitted with a double-exponential and exported from LASX (as 'raw ImageJ tif'). 
  - TauSeparation, directly from the `.lif` file
  - Fast FLIM
  - TauContrast
  - Lambert Instruments Frequency Domain FLIM (`.fli` files)

  Supported non-FLIM images are:
  - Ratio Imaging (2-channel `.tif` files)
  - Intensity (single-channel `.tif` files)

The macro can batch-process multiple files as well.

## 2. Pre-processing
Optionally, image corrections can be performed as preprocessing steps:
1. Stitching of mylti-tiled experiments, with the separate command `Stitch tiles` (tested for Leica `.lif` files). Fast 'no-calculation' stitching is done using the [Grid/Collection Stitching plugin](https://imagej.net/plugins/grid-collection-stitching) to ensure correct pixel size[^1] and to prevent image distortion. XY-Stage coordinates of image tiles are extracted from the OME metadata of the `.lif` file and stored in the stitched `.tif` file, together with the grid layout, the tile size and the tile overlap percentage.
2. Drift correction. A cross-correlation (cc) of image frames is performed with either the first frame, the previous frame or the last frame. CLIJ2 image projection functions are used to quickly determine the coordinates of the maximum pixel in the cross-correlation image, its distance from the center representing the shift. Drift correction is performed on the calculated intensity image (i.e. the addition of the two components, if applicable).
3. Bidirectional scanning phase mismatch correction. The even and odd lines are split and turned into two separate images. Cross-correlating these images and locating the peak yields the shift between the phases. The two images are both shifted half this distance and interleaved again:
[^1]: The exported pixel size is not entirely correct, causing an increasing deviation in the cell coordinates for increasing x and y.

![image](https://github.com/user-attachments/assets/66408493-ec41-4c4b-9413-3d6ae136e932)

After these optional steps, a 'weighted lifetime' image is created. For single-component FLIM images ('FAST FLIM', 'TauContrast', Frequency-domain FLIM) this is simply the lifetime; for a two-component FLIM image this is equal to the *intensity-weighted* lifetime. (Note that the intensities here are the intensities of the lifetime components!)
Additionally, an 'RGB overlay image' is created, where a (x,y,t) smoothed version of this lifetime image is multiplied (overlayed) by the total intensity image, yielding a denoised visualization of the experiment.  

![Cos7H250_ADRB2KO_1 (weighted lifetime)](https://github.com/user-attachments/assets/c66ad486-1d79-4c6b-9b3c-e95d7b4d6683)
![Cos7H250_ADRB2KO_1 (lifetime   intensity RGB overlay)](https://github.com/user-attachments/assets/3aba7cb8-d64b-4d65-964f-4af346cb6eb5)

## 3. Cell segmentation
cell segmentation can be performed with [CellPose](https://github.com/MouseLand/cellpose) (2 or 3), operated from Fiji using [a wrapper](https://github.com/BIOP/ijl-utilities-wrappers). The image stack is first 'collapsed' using a summed-intensity projection of a subset or all of the time-lapse images, resulting in a single imageto segment. This procedure works well if the imaging is short enough that the cells do not move (a lot).
The user can select the segmentation model (pretrained or custom) and needs to provide a few key parameters, e.g. cell diameter and [flow threshold](https://cellpose.readthedocs.io/en/v3.1.1.1/settings.html#flow-threshold). Additionally, restrictions on cell size and circularity can be imposed.

![Cos7H250_ADRB2KO_1 (lifetime   intensity RGB overlay)](https://github.com/user-attachments/assets/eb3099c0-e9e6-4d8c-8000-a7d5b41892c9)

## 4. Measuring single-cell (FLIM/intensity) traces
After cell segmentation ROIs are created from the obtained label image, after which intensities and average fluorescence lifetimes are computed for every cell, at every time point. This average lifetime is the weighted lifetime, where each pixel of a cell is linearly weighted with its intensity fraction.

The script automatically determines the time points of a(nta)gonist stimulation and calibration by detecting peaks in the second derivative of the average trace of all cells. If the peaks are higher than a set number of times the stddev of the signal it is picked up. If successful, cell traces are divided into three parts: *baseline*, *response*, and *calibration*. If not, manual input of the time points is also possible. These three partitions are used for detection of hit cells when screening for dynamic phenotypes.

<img src="https://github.com/user-attachments/assets/52832275-d256-4162-8a86-7ee22ab6f2df" title="first and second derivative trace " width="500">

## 5. Visualization
The data is visualized in various graphs and images:
### Time traces plot
<img src="https://github.com/user-attachments/assets/a23f8818-b754-4b01-93ab-b5e1a6dff63e" title="lifetime traces plot" width="510">

### Timelapse Lifetime histogram and scatterplots
<img src="https://github.com/user-attachments/assets/87c92c02-cf5d-4e7f-b391-40c858a6ffa9" width="500">
<img src="https://github.com/user-attachments/assets/d750ebd6-a1d4-4687-83cc-bfae930cc858" width="500">

### Kymographs
This is an image with time as Y-coordinate, Cell ID as X-coordinate and cell lifetime as value. Additionally, a 'sorted kymograph' is created, where the X-axis is sorted on the average response lifetime.

![Cos7H250_ADRB2KO_1 (kymograph)](https://github.com/user-attachments/assets/5c64cad6-374c-481d-8b0c-273819a7bcce)
![Cos7H250_ADRB2KO_1 (kymograph sorted)](https://github.com/user-attachments/assets/bd32872a-b9af-4a52-b63e-8e685f17850e)


# 6. Screening: hit selection

### Hit criteria panel
![image](https://github.com/user-attachments/assets/4d2f37da-727e-4e86-a424-0fa101b05ba6)

### Time traces plot of hits
The time window used in the hit selection is highlighted in blue.
![image](https://github.com/user-attachments/assets/10db5f6a-a2f5-4ba7-9ba2-50f826b496a7)
![image](https://github.com/user-attachments/assets/be2f6ab5-ea0c-4646-97f4-77cfa6519d20)


