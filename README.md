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
