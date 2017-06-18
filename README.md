# EyeStalker: robust video-based eye tracking using recursive estimation of pupil characteristics

Video-based eye tracking is a valuable technique in many research fields. Numerous open-source eye tracking algorithms have been developed in recent years, primarily designed for general application with many different camera types. However, these algorithms do not capitalize on the high frame rate of eye tracking cameras often employed in psychophysical studies. We present a pupil detection method that utilizes this high-speed property to obtain reliable predictions through recursive estimation about certain pupil characteristics in successive camera frames. These predictions are subsequently used to carry out novel image segmentation and classification routines to improve pupil detection performance. Based on results from hand-labelled eye images, our approach was found to have a greater detection rate, accuracy and speed compared to other recently published open-source pupil detection algorithms.

<b>GUI</b>

![alt text](https://cloud.githubusercontent.com/assets/10850074/26767383/ef50d9d2-499f-11e7-858e-4c08660d4b82.png)

If you use this software, please cite:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.810495.svg)](https://doi.org/10.5281/zenodo.810495)

or:

@misc{terence_brouns_2017_810495,
  author       = {Terence Brouns},
  title        = {tbrouns/eyestalker: Linux release v.1.0.1},
  month        = jun,
  year         = 2017,
  doi          = {10.5281/zenodo.810495},
  url          = {https://doi.org/10.5281/zenodo.810495}
}

Copyright (C) 2017 Terence Brouns 
