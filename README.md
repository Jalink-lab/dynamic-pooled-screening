# FATSHOPS - FLIM Analysis of Traces and Selection of Hits for Optical Pooled Screening

FATSHOPS is a versatile Fiji macro collection to analyze and visualize multi-cell time-lapse experiments, detect hit cells based on user-set criteria and output their stage coordinates for e.g. photoactivation.
For non-screening applications it can function as a valuable tool for single-cell trace analysis, visualization and inspection.

## Simplified workflow
- *[Loading images](https://github.com/Jalink-lab/dynamic-pooled-screening/blob/main/README.md#input-images)*
- *[Pre-processing](https://github.com/Jalink-lab/dynamic-pooled-screening/blob/main/README.md#pre-processing)*
- *[Cell segmentation](https://github.com/Jalink-lab/dynamic-pooled-screening/blob/main/README.md#cell-segmentation)*
- *Measuring single-cell (FLIM/intensity) traces*
- *Visualization*
- *Hit selection*
- *Output files*

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


## Input images
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


## Pre-processing
Optionally, image corrections can be performed as preprocessing steps:
1. Drift correction. A cross-correlation (cc) of image frames is performed with either the first frame, the previous frame or the last frame. CLIJ2 image projection functions are used to quickly determine the coordinates of the maximum pixel in the cross-correlation image, its distance from the center representing the shift. Drift correction is performed on the calculated intensity image (i.e. the addition of the two components, if applicable).
3. Bidirectional scanning phase mismatch correction. The even and odd lines are split and turned into two separate images. Cross-correlating these images and locating the peak yields the shift between the phases. The two images are both shifted half this distance and interleaved again:

  ![image](https://github.com/user-attachments/assets/66408493-ec41-4c4b-9413-3d6ae136e932)


## Cell segmentation
cell segmentation can be performed with [CellPose](https://github.com/MouseLand/cellpose) (2 or 3), operated from Fiji using [a wrapper](https://github.com/BIOP/ijl-utilities-wrappers). The image stack is first 'collapsed' using a summed-intensity projection of a subset or all of the time-lapse images, resulting in a single imageto segment. This procedure works well if the imaging is short enough that the cells do not move (a lot).
The user can select the segmentation model (pretrained or custom) and needs to provide a few key parameters, e.g. cell diameter and [flow threshold](https://cellpose.readthedocs.io/en/v3.1.1.1/settings.html#flow-threshold). Additionally, restrictions on cell size and circularity can be imposed.

## Visualization
### Time traces plot
![image](https://github.com/user-attachments/assets/0e74c287-5c79-4c21-8a04-5f3ba5db2ff9)

### Scatter plots
Create movies
![image](https://github.com/user-attachments/assets/795e8728-2183-4817-8853-4252df1c7b67)

### Screening: hit selection


### Hit criteria panel
![image](https://github.com/user-attachments/assets/4d2f37da-727e-4e86-a424-0fa101b05ba6)

### Time traces plot of hits
The time window used in the hit selection is highlighted in blue.
![image](https://github.com/user-attachments/assets/10db5f6a-a2f5-4ba7-9ba2-50f826b496a7)
![image](https://github.com/user-attachments/assets/be2f6ab5-ea0c-4646-97f4-77cfa6519d20)

### derivative trace
![image](https://github.com/user-attachments/assets/52832275-d256-4162-8a86-7ee22ab6f2df)
