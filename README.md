# Chromosome-analyze
## Overview

Chromosome-analyze is a set of tools to segment and trace the axes of pachytene chromosomes in three-dimensional images of meiotic nuclei acquired by fluorescence microscopy. It was designed for images of immunofluorescently stained germline tissue in C. elegans, but is applicable to images acquired using other imaging modalities (e.g., confocal microscopy) in other organisms, provided that the axes of individual chromosomes are resolvable in the raw images. 

## Data format and organization

The basic unit of raw data is a set of three-dimensional grayscale images, one for each fluorescent channel, corresponding to a single field of view (FOV). Usually, two or three fluorescent channels are present: one for chromosome axis labeling, and one or two for other features of interest (e.g., DAPI staining, FISH probes, DSB markers, etc). 

The raw images may be stored in a variety of proprietary formats. The code here was designed to operate on raw images stored in DeltaVision (DV) format, but the first step is to export these images as TIFF stacks, so other raw image formats can be accommodated as long as they are somehow re-saved as TIFF stacks.

In C. elegans, each FOV is of a portion of a single gonad from a single animal. The FOVs from a single animal are collected together into a single subdirectory; all of the subdirectories corresponding to all of the animals imaged from a single sample (i.e., batch of animals stained together but possibly mounted on multiple slides) are in turn collected into one directory which here is called the 'data directory.' 

The naming of this data directory is arbitrary, but here, the following convention is used:
    
    \\[year]_[month]_[day]_[genotype]_[experiment or staining description]

For example, for images of RAD-51 staining in wild-type (N2) captured on 7/22/2014, the data directory would be:

    \\2014_07_22_N2_RAD51\\ 

then the images from each animal from this sample are found in 

    \\2014_07_22_N2_RAD51\\01
    \\2014_07_22_N2_RAD51\\02

And so on. Within each of these subdirectories, the raw images (DV files) are named with the data directory name, followed by the subdirectory number (i.e., animal/gonad number), followed by the FOV number within that animal, followed by microscope-specific postfixes. So, for example:

    \\2014_07_22_N2_RAD51\\04\\2014_07_22_N2_RAD51_04_01_SIR_ALX.dv

This is FOV 01 in animal 04 (in subdirectory \\04) from the RAD-51-stained N2 strain imaged on 7/22/2014. The _SIR_ALX tag identifies the image as the final output from a structured illumination microscope (such images are used as examples here, but images acquired from other microscopes can be substituted as discussed above). There can be an arbitrary number of FOVs in each subdirectory.

## Writing TIFFs 

Given this data organization, the first step in MATLAB, assuming that one is starting from raw DV files, is to run the command

    batchBatchWriteForSubdirs(rootDir, filter_string, write_files_flag)

where rootDir is the data directory; e.g., `rootDir = E:\\Data\\2014_07_22_N2_RAD51_A02\\`

and filter_string is specific to the raw DV files; e.g., `filter_string = 'SIR_ALX'`

The function `batchBatchWriteForSubdirs` looks in each subdirectory and runs two functions:

    batchWriteTIFFfromOMX(rootDir, writeProjection, writeStack, filter_string)

This function opens each of the DV files (which it identifies using `filter_string`) in the data directory and writes each channel of each DV file as both a TIFF projection or TIFF stack. Each channel is assigned a single-digit single-character number corresponding to the order in which the channels appear in the DV file. This order is, in turn, usually determined by wavelength, but may also be microscope-specific. In practice channel identities are usually obvious by observing the image itself. 

Next, `batchBatchWriteForSubdirs` runs

    batchWriteIsotropicTIFFfromOMXTIFF(rootDir, z_range, filter_string)

This function acts on each subdirectory one by one and, in each, loads the TIFF stacks written by `batchWriteTIFFfromOMX`. It then makes a new sub-subdirectory for each FOV (named after the corresponding DV file) and writes isotropically resampled TIFF stacks for each channel from each image in this sub-subdirectory. These are simply named im1.tif, im2.tif, etc, where the number corresponds to the channel number assigned above.

The isotropic TIFF stacks are resampled only in z in order to yield isotropic voxel dimensions. In the images used as examples here, the raw voxel depth is 125nm, while the x-y pixel size is 39.6nm. The TIFF stacks are therefore resampled in z so that the number of slices is increased to `ceil([number of z slices]*125/39.6)`. (The images are usually about 4um or 32 slices deep, so this introduces an error of at most `1/(32*125/39.6) = 1%`.) 

***VERY IMPORTANT***: The dimensions of the raw voxels (that is, the x-y pixel size and the z-slice spacing) are hard-coded in the beginning of `batchWriteIsotropicTIFFfromOMXTIFF` in the line

    PX_SZ = [0.0792/2 0.0792/2 0.125]*1000;

**These numbers are specific to the example images used here. They are STRONGLY microscope-dependent and MUST BE CHANGED accordingly.** They are hard-coded, rather than passed as an argument, to ensure that they cannot easily or accidentally be changed---but this conversely requires that the user be aware of their importance. 

The path to the isotropic TIFF stacks is then (using the same example filename as before):

    \\2014_07_22_N2_RAD51\\04\\2014_07_22_N2_RAD51_04_01_SIR_ALX\\im1.tif
    \\2014_07_22_N2_RAD51\\04\\2014_07_22_N2_RAD51_04_01_SIR_ALX\\im2.tif

For example, here `im1.tif` might correspond to chromosome axis staining, and `im2.tif` to RAD-51 staining. These sub-subdirectories (corresponding to single FOVs from single animals) are the 'units' of analysis for the rest of the analysis.

Finally, `batchBatchWriteForSubdirs` prints a list of OMX image directories in the format

   dds(nn).name = E:\\Data\\2014_07_22_N2_RAD51_A02\\04\\2014_07_22_N2_RAD51_04_01_SIR_ALX\\

This list can be regenerated at any time by rerunning the function with `write_files_flag = 0`. Its use is explained below.

## Choosing parameters and other settings

The next step, now that isotropic three-dimensional images, in the form of TIFF stacks, have been obtained, is to make two new data-directory-specific files, both named after the data directory itself. The nomenclature is again arbitrary but here adheres to the following:

    batch_[root directory].m
    batch_[root directory]_settings.m

For example:

    batch_2014_07_22_N2_RAD51.m
    batch_2014_07_22_N2_RAD51_settings.m

The purpose of the first .m file is only to return a cell array of image sub-subdirectory names. This is done by simply copying and pasting the list of directories returned by `batchBatchWriteForSubdirs`, so that `batch_[root directory].m` looks like this:

    function dds = batch

    dds(1).name = 'E:\Data\2014_07_22_N2_RAD51\01\2014_07_22_N2_RAD51_01_SIR_ALX';
    dds(2).name = 'E:\Data\2014_07_22_N2_RAD51\01\2014_07_22_N2_RAD51_02_SIR_ALX';

and so on. 

The `_settings.m` file contains all of the parameters necessary for the nucleus segmentation, axis segmentation, and focus finding (but NOT the axis tracing). These are all hard-coded and are copied from previous settings files and adjusted as necessary. 

Before this function can be called, another function is called first to initialize the MATLAB data file that will contain all of these settings in each FOV sub-subdirectory. This file, `initializeGonadOMX.m`, writes a file called `gonad.mat` in each FOV directory. These files store all of the settings that will be applied by the `_settings.m` file as well as basic, common settings that are hard-coded in `initializeGonadOMX.m` (many of which are overwritten by settings in the `_settings.m file`---their existence in `initializeGonadOMX.m` is at this point archaic and indeed some other settings in initializeGonadOMX are no longer used by the code).

Important parameters hard coded in `initializeGonadOMX` are the voxel size of the isotropically resampled images (this should simply be the microscope's x-y pixel size).

The function to run `initializeGonadOMX` on each directory is

    batchInitializeGonadOMX(dds, num_chan, seg_chan)

where `dds` is the output of `batch_[root directory].m` and `num_chan` is the total number of channels in each FOV (usually only two or three) and `seg_chan` is the channel in which the axes are imaged (i.e., 1 for im1.tif, 2 for im2.tif, etc). 

Now, we can call the `_settings.m file`, which also takes the `dds` variable as its input. It opens the `gonad.mat` files created by `initializeGonadOMX` and adds/overwrites the hard-coded settings. 

The important settings it sets are:

    % segmentation settings
    gonad.SEG_SMOOTHING    = 3;
    gonad.SEG_GAMMA        = .6;
    gonad.SEG_THRESH       = 1.7;
    gonad.useGamma         = 1;
    gonad.useOtsuThresh    = 1;
    gonad.MIN_REG_VOL      = 25;
    gonad.MIN_SURF_AREA    = 1;

    % focus-finding settings
    fATable = zeros(6,3);
    gonad.focusAxisTable    = fATable;
    gonad.focusChannels     = [0,   2];
    gonad.focusSmoothing    = [0,   1];
    gonad.focusThresh       = [0,  10];

The segmentation settings determine how the axis image is manipulated before the initial mask is generated (upon which watershed boundaries are superimposed to generate regions and surfaces). 

`seg_smoothing` determines the radius of the Gaussian used to blur the image.

`seg_gamma` determines the gamma applied to the intensity (i.e., just image.^gamma). Usually 0.7 or 0.6 is best to retain the dimmer coverslip-distant axes without blowing out the bright coverslip-proximal axes. 

`seg_thresh` is not used if useOtsuThresh is enabled, which for OMX images it always should be---this option lets the segmentation code use Otsu threshholding to determine the right threshhold value (which is always a multiple of the mean image intensity). 

`Min_reg_vol` and `min_surf_area` are used in the segmentation code and should not need to be changed.

The focus-finding settings are a bit confusing. The `fATable` is archaic at this point and should always be set to zero. The `focusChannels` vector lists the channel numbers of the channels in which there are foci to be found---usually just one channel in OMX images, so that the vector is either [0,1] or [0,2]. There is no difference between [0,1] and [1,0], except that the order has to be consistent between `focusChannels`, `focusSmoothing`, and `focusThresh`. These vectors describe the smoothing radius and relative threshold---always 1 for smoothing and between 7 and 10 for threshold. 

Finally, the `_settings.m` file for OMX images also contains a few lines that use `estimateIlluminationProfile.m` to attempt to compensate for the uneven illumination present in most OMX images. This is not very successful but also in the worst case doesn’t seem to hurt either. 

## Nucleus segmentation
Once the `_settings.m` file has run, we are finally ready to segment the nuclei. The first step is to run the `makeBackgroundMask` function on all of the images, which is accomplished by 

    batchMakeBackgroundMask(dds);

This generates a bandpass filtered image, which is saved in each FOV subdirectory. Note that an estimated illumination profile compensation image (just a wide Gaussian that, multiplied by the OMX images, brightens the dimmer periphery of the field of view relative to the center) is hard-coded in this function. 

The nucleus mask itself is made by thresholding this filtered image using a threshold chosen by the user in segSupervisor. The command is

    segSupervisor(dds(ii).name)

This call a GUI and will also generate the nucleus mask for that FOV using `makeNukesFromMask4`, once an appropriate threshold is chosen by the user (the threshold can be previewed; usually a threshold between 18 and 24 is good, but there is a lot of variation between images). The ‘segment nuclei’ button calls this function, which identifies each isolated region in the masked image and assigns it to a nucleus. It uses a kind of k-means clustering. Running `makeNukesFromMask4` takes some time, so it is best to use `segSupervisor` to set the threshold for each image (this threshold is stored in gonad.mat) and then run `makeNukesFromMask4` on all of the images at once by calling

    batchMakeNukesFromMask4(dds);

There are invariably errors in the nucleus segmentation, so it is necessary to run `segSupervisor` again to correct these errors using the split, merge, and resegment buttons. 

Once a suitable nucleus mask is generated, adjust the min and max volumes, and the maximum image-edge-adjacent nucleus surface area, so that all of the “good” nucleus masks will be included (have green numbers, not red). Inclusion is defined by `isNuke2.m`.

Once this procedure has been completed for all of the OMX images, the function `batchMakeNukeData` runs the function `makeNukeData` on each OMX image directory. This function constructs the `nd.mat` data file that contains cropped images of each nucleus based on the nucleus mask from the `segSupervisor` output. It uses `isNuke2.m` to decide whether each region in this mask is in fact a nucleus. One `nd.mat` file is generated for each nucleus and is stored in a nucleus-specific subdirectory (also generated by `batchMakeNukeData`). These subdirectories are all in the /nukes/ subdirectory of each FOV subdirectory. Finally, run patchGonad_2013_07_24(dds) to update the nucleus properties in `gonad.mat`.

## Chromosome segmentation and focus finding

### Choosing parameters
The cropped image data for individual nuclei in a given FOV can be visualized using a GUI called nukeViewer. The command is the same as for segSupervisor:

    nukeViewer(dds(ii).name);

The code that attempts to segment the axes without supervision is contained within two functions. The first, `segAxis3DV`, generates an axis mask, segments the mask using the watershed transform, and then cleans up the resulting 3D regions and surfaces. The second function, `optiAxis3SmartFast2v3`, attempts to merge regions until there are six regions that correspond to the six synapsed chromosome axes. It does this by using a region and surface score to determine which surfaces between adjacent regions are most likely to be spurious. This procedure is more effective on confocal images (for which it was developed) than on DV or OMX images, which suffer from a greater range of “real” intensities and also from greater heterogeneity in intra-axis staining intensity. Fortunately, a good OMX image usually has relatively well-differentiated axes and needs little optimization. 

Both of these functions are run on all nuclei associated with the images in the `dds` file list by the function `batchAnalyzeOMX`. But, before doing this, it’s important to use `nukeViewer` to determine the right segmentation parameters. 

The first of these is the gamma used to generate the axis mask in `segAxis3DV`. It is necessary to adjust the raw image using a gamma < 1 in order to reduce the intensity gap between bright and dim axes in the same nucleus. Usually something between 0.6 and 0.8 works well---but experimentation in `nukeViewer` is necessary using the “view background mask” button.

The second of these parameters is the choice to use either an absolute threshold or Otsu thresholding to generate the chromosome axis background mask. Usually, an absolute threshold works well, but the right value has to be determined by inspection (usually in coordination with the gamma). 
 
Finally, if there are foci in the image (for example, RAD-51 foci), these can be identified automatically by another function, `spotAxis3`, that is also run on each nucleus by `batchAnalyzeOMX`. This function requires a few parameters that should also be determined empirically in nukeViewer. These are the relative intensity threshold (how bright a spot has to be relative to the background) and the smoothing radius. Usually a threshold of around 10 works, but it’s important to check using the ‘Find Foci’ button in nukeViewer. The smoothing radius is almost always set to 3. 

Once these parameters - the gamma and threshold for the axis image and the intensity cutoff for the focus finding - are determined in nukeViewer and checked for a few images in the data set, they must be applied to each image’s `gonad.mat` settings file. To do this, adjust the gamma and focus intensity threshold accordingly in the `_settings.m` file and run it again. Note that this will not overwrite any custom nucleus segmentation settings applied by segSupervisor (though they are not useful at this point). 

Finally, the function `batchAnalyzeOMX(dds)` can be called. Make sure that the `DO_SEG` and `DO_FOCUS_FIND` flags, and only those flags, are set to 1. This code will take quite a long time to run - the time seems to be on the order of tens of seconds per nucleus, so for 1000 nuclei, at least five hours will be required.

### Correction of segmentation in nukeViewer
The estimated axis segmentation needs to be verified by hand, and also, for RAD-51 analysis, “good” axes that are both correctly segmented and appear to have real foci need to be selected. This is accomplished in nukeViewer. The segmentation is corrected by selecting/deselecting surface ids in the textbox (always while holding down CTRL). Good axes are selected by clicking on the colored boxes. After finishing with each nucleus, be sure to click the ‘finalize nucleus’ button to save the changes to the nd.mat file. 

### Concatenation of correctly segmented axes
Once all of the nuclei in a data set have been checked and ‘good’ axes selected, the axis data is concatenated. This is accomplished by 

    axisData = batchCatAxesForRAD51OMX(dds);

The output, `axisData`, is a cell array with entries for each axis. This code uses `makeAxisDataForFISH` to construct the entry for each axis. Its most important function is to construct a mask and skeleton for each axis. 

## Chromosome tracing
The axes in `axisData` are traced using the Snakes algorithm by running

    axisData = batchSnakeTraceRAD51OMX(axisData);

This uses `snakeTrace3v1` with parameters chosen for OMX images. These parameters are hard coded in `defineSnakeOptsForRAD51OMX` and are resolution-dependent. This code takes a few seconds per axis to run. 

## Calculating focus positions
The final step in the automated analysis is to calculate the axial positions of the foci, if there were any. This code calculates the distance between each axis trace and all of the foci in that axis’s nucleus. It also generates a unique ID for each spot, which is used later eliminate double-counting of the same spot (since all of the code operates on axisData, and not a list of spot data, it’s possible for a spot that’s sufficiently close to two axes to be counted twice). The line is

    axisData = calcFocusPositionsRAD51(axisData, 0, 0);

The 0s are no longer used and present only for backwards compatibility. This function adds new fields to `axisData` that give the positions of each spot, its axial position along the axis trace generated by `batchSnakeTrace`, and intensity. Note again that all of the spots in each axis’s nucleus are included here, so that whenever two or more axes in `axisData` are from the same nucleus, this spot data is repeated in the `axisData`. 



