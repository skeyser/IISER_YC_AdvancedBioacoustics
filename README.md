# IISER_YC_AdvancedBioacoustics

Repository for the IISER-Yang Center Advanced Bioacoustics workshop in Tirupati, India.

## Instructors

### Yang Center for Conservation Bioacoustics, Cornell Lab of Ornithology

- [Rebecca Cohen, Ph.D.](https://www.birds.cornell.edu/ccb/rebecca-cohen-ph-d/)
- [Irina Tolkova, Ph.D.](https://www.birds.cornell.edu/ccb/irina-tolkova-ph-d/)
- [Spencer Keyser, Ph.D.](https://www.birds.cornell.edu/ccb/spencer-r-keyser-ph-d/)
- [Meghan A. Beatty, Ph.D.](https://www.birds.cornell.edu/ccb/meghan-a-beatty-ph-d/)

### IISER

- V.V. Robin, Ph.D.
- Isha Bopardikar, Ph.D.
- Chiti Arvind

## Workshop Overview

This workshop is intended as a practical, hands-on opportunity for bioacoustics practitioners to gain advanced analytic skills that enable efficient handling of large bioacoustics data sets and their use for addressing behavioral and ecological research questions. 
Some familiarity with basic concepts in acoustics, bioacoustics, and digital signal processing is expected, though an overview of key fundamentals will be provided. An opening symposium will introduce the workshop organizers and instructors and their research, and allow participants to get to know one another. 
Participants will move through several modules together, each of which will include lecture, case studies, and hands on analysis of provided demonstration data sets. A mid-week field trip will provide opportunity for participants to visit some local field sites and address practical considerations of passive acoustic data collection.
There will also be opportunity for participants to apply analytic concepts covered in the workshop to their own data sets with guidance and troubleshooting by workshop instructors, and a closing session in which participants can present their works-in-progress. 
Analyses will be carried out in R, and all code and demonstration data sets will be made publicly available following the workshop to provide reference and teaching materials.

## Modules

### *Basic Acoustics Review & Data Management*

A review of basic concepts in bioacoustics and digital signal processing, including: physics of sound and acoustic propagation; the SONAR equation; sound levels; decibels; analog vs digital signals; sampling considerations (sampling rate, bit depth, recording schedule); calibration; Fourier transforms; sound visualization (spectrograms & LTSAs); signal measurements. 
Best practices in acoustic data management will be covered, including metadata, file formats (uncompressed vs compressed), file naming, version controlling, and data storage.

### *Detection/Classification and Machine Learning*

A defining characteristic of passive acoustic monitoring is the ability to collect large data sets (10’s GB - 100s TB) with a minimum of field effort. While these rich data sets can be used to address myriad research questions, their size and complexity necessitate specialized analytic approaches to efficiently extract meaningful insights. 
Automated detection and classification techniques are commonly employed to identify signals of interest in large passive acoustic data sets much more quickly and objectively than could be achieved via manual review and annotation. Recently, machine learning algorithms have gained popularity for detection and classification tasks, as well as for addressing needs such as dimensionality reduction and clustering. 
In this module we will review some popular automated detectors/classifiers available in popular acoustic data analysis software, delve into the basics of machine learning, and train and evaluate a custom machine learning detector/classifier. 

### *Localization*

If a vocalization is heard across multiple recording units, it may be possible to localize signals in space, enabling analysis beyond detection and classification. This module will describe the fundamentals of acoustic localization and direction-of-arrival estimation,  demonstrate these methods through several case studies, and describe connections to abundance estimation. Finally, we will discuss the associated challenges, such as synchronization, call association, and quantifying uncertainty.

### *Occupancy Modelling and Density Estimation*

Robust estimates of species geographic distributions, habitat associations, abundances, and population trajectories are critical to basic and applied ecology and conservation. However, many animals are cryptic, rare, and otherwise difficult to detect such that our own observation process is imperfect. Occupancy models provide a robust, hierarchical modelling framework that explicitly accounts for the imperfect observation process by leveraging repeated surveys, which bioacoustic survey methods are optimized for. 
Similarly, the occupancy modelling framework can be adapted to integrate continuous data on species counts providing a natural extension to estimating abundance that is critical for conservation decision making. In this module, we will explore the fundamentals of occupancy modelling from theory to application, provide hands-on practice for running basic occupancy models (i.e., single-species, single-season) in R, 
and provide examples of extensions to relax assumptions of the basic occupancy model case to explore advanced applications (e.g., multi-species occupancy models; dynamic occupancy models). Next, we will focus on extending concepts from species detection/non-detection data to species count data for estimating abundance after accounting for imperfect detection, or so-called N-mixture models. We will focus on data derived from bioacoustic surveys and facilitate group discussions on the benefits, limitations, and applications of bioacoustic data in the context of occupancy and abundance modelling. 

## Software Installations

The following programs are used in the demos and hands-on sections of the various workshop modules. Workshop participants should download and install these prior to the beginning of the workshop.

### BirdNET Analyzer
BirdNET Analyzer is a GUI for processing audio data, with functionality including supervised deep learning detection/classification using the BirdNET model, custom model development using BirdNET as a foundation, detection review and classifier performance evaluation, etc. BirdNET and BirdNET Analyzer are developed and maintained by the K. Lisa Yang Center for Conservation Bioacoustics at Cornell University, and Chemnitz University of Technology.
  1. Download BirdNET Analyzer v2.4.0 installer from the [BirdNET GitHub repository releases page](https://github.com/birdnet-team/BirdNET-Analyzer/releases/tag/v2.4.0)
  2. Double click the downloaded installer file to run it; you may need to tell your computer that it is safe to run. Click **More info** and **Run anyway**.
  3. Accept the license agreement and click **Next**. Check the box to create a desktop shortcut and click **Next**. Click **Install**.

### Raven Annotate
Raven Annotate is a new [Raven Workbench](https://www.birds.cornell.edu/ccb/raven-pro/) module that provides functionality for loading and visualizing raw sound files and annotating acoustic events. It is developed and maintained by the K. Lisa Yang Center for Conservation Bioacoustics at Cornell University.
  1. Use one of the following links to download the appropriate Raven Annotate installer, according to your operating system:
    Windows: https://updates.ravensoundsoftware.com/updates/workbench/raven_annotate/win_x64/Raven-Annotate-1.0.1.exe
    MacOS: https://updates.ravensoundsoftware.com/updates/workbench/raven_annotate/mac_arm64/Raven-Annotate-1.0.1.dmg
    Linux: https://updates.ravensoundsoftware.com/updates/workbench/raven_annotate/linux_x64/raven-annotate_1.0.1-1_amd64.deb 
  2. Locate the installer in your **Downloads** folder and double click to run.
  3. Follow the prompts to complete installation by accepting the license agreement and clicking **Install**.

