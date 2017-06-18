# EyeStalker: robust video-based eye tracking using recursive estimation of pupil characteristics

Video-based eye tracking is a valuable technique in many research fields. Numerous open-source eye tracking algorithms have been developed in recent years, primarily designed for general application with many different camera types. However, these algorithms do not capitalize on the high frame rate of eye tracking cameras often employed in psychophysical studies. We present a pupil detection method that utilizes this high-speed property to obtain reliable predictions through recursive estimation about certain pupil characteristics in successive camera frames. These predictions are subsequently used to carry out novel image segmentation and classification routines to improve pupil detection performance. Based on results from hand-labelled eye images, our approach was found to have a greater detection rate, accuracy and speed compared to other recently published open-source pupil detection algorithms.

<b>Paper</b>

Will be added soon. 

Novel hand-labelled pupil detection dataset of 5000 images from 12 participants, as mentioned in paper, can be downloaded from [here](https://drive.google.com/open?id=0Bw57olSwQ4EbUWV5ajNKeG93NEk).

<b>GUI Layout</b>

![alt text](https://cloud.githubusercontent.com/assets/10850074/26767383/ef50d9d2-499f-11e7-858e-4c08660d4b82.png)

<b>Quick start</b>

Follow the install instructions on "Releases" page (https://github.com/tbrouns/eyestalker/releases/latest). After starting up the application, a GUI similar to the one shown above should appear, but without the displayed eye image. To see the eye tracking algorithm in action, we will test it on a sample data set. You can download it from here (*data.zip*):

[Download data set here](https://drive.google.com/open?id=0Bw57olSwQ4EbOG9kVTAzUjBtNTA)

This data set consists of a 1500 ms recording of a saccadic eye movement recorded at 250 Hz (375 images in total). After extracting the ZIP file, your directory should look as follows:

```
data
│
└───images
│   │   
│   └───trial_0
│       │   
│       └───raw
|           | 0.png
|           | 1.png
|           | ...
```

You can give the *data* directory any name you wish, but the subdirectories and filenames must not be altered. 

In the GUI, press *Load session* and select the *data* directory. We will now perform pupil detection on the whole data set by clicking *All frames*. The program will run through every image and draw a teal ellipse on the pupil-iris boundary together with a cross which marks the pupil centre, if it has successfully detected the pupil.  The result is immediately visible in the display frame. Once detection is completed, you can move through all the images with the slider (directly to the right of  *Combine*) to see the result for each individual frame. 

In the *trial_0* directory, a new subdirectory will have been created called *processed*, which contains all the processed images for display purposes. Furthermore, there will be a DAT file called *tracking_data.dat* that contains the eye tracking measurements. The DAT file contains a single row of data. The first value gives the number of samples, which is 375 for the sample data set. This is followed by 5 concatenated data vectors, each containing 375 elements. These are:
1. The first vector comprises of ones and zeroes, indicating whether the pupil was detected in the corresponding camera frame (Yes = 1, No = 0).
2. Pupil centre x-position in image coordinates (pixels)
3. Pupil centre y-position in image coordinates (pixels)
4. Pupil circumference (pixels)
5. Pupil aspect ratio (ratio between major and minor axes)

You can play around with the various parameters in the *Eye tracking* tab. You can press *One frame* to see the effect of a change in parameter value on pupil detection in the current camera frame. 

If you select *Box* and *Edges*, the Haar-like feature detector and Canny edges will also be drawn in the procesed image, respectively. 

Pressing the *Quit* button will ensure that a *config_user.ini* file is saved in the same directory as the AppImage, which contains your current parameter configuration. This INI file is automatically loaded next time you start the application. To reset all parameters you can remove the INI file and restart the application or press the *Reset parameters* button. 

More details will be added soon.

<b>Citation</b>

If you use this software, please find the citation here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.810494.svg)](https://doi.org/10.5281/zenodo.810494)

Or directly for BibTex: https://zenodo.org/record/810494/export/hx

<b>Third-party libraries</b>

EyeStalker has been built using:

* Qt, cross-platform application framework (The Qt Company, https://www.qt.io/, version 5)
* OpenCV, open-source computer vision library (Intel Corporation, Itseez Inc., http://opencv.org/)
* Eigen, C++ template library for linear algebra (http://eigen.tuxfamily.org)
* Boost, peer-reviewed portable C++ source libraries (http://www.boost.org/).
* QDarkStyleSheet, a dark stylesheet for Qt applications (https://github.com/ColinDuquesnoy/QDarkStyleSheet)

EyeStalker been deployed using *linuxdeployqt*, available from: https://github.com/probonopd/linuxdeployqt

Copyright (C) 2017 Terence Brouns 
