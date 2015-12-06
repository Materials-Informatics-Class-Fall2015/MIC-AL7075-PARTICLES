---
layout:     	post
title:      	Final Post
date:       	2015-12-06 15:00
author:     Paul Kern and Chris Shartrand
tags:        result 
---
<!-- Start Writing Below in Markdown -->

**Introduction and Background**
Aluminum Alloy 7075 (AA 7075) is an aluminum alloy whose composition is composed largely of 5.1-6.1% zinc, 2.1-2.9% magnesium, 1.2-2.0% copper and <0.5% Iron, Silicon, Manganese, Chromium, Zirconium and Titanium. During the material processing and hot rolling stages of the alloy, it is common for μm-sized particles to become ingrained in the material matrix. Rollett Et al. have stated that these particles contribute up to 2% of the total volume of the matrix. Deformation of the alloy can cause these particles to either crack or become fully detached from the matrix. When the alloy is eventually put under stress, these cracked particles and the areas where the particles have detached can cause the matrix itself to crack. A crack in the matrix combined with continued stress can cause a failure within the section of the alloy.
![An example of a cracked particle](/MIC-AL7075-PARTICLES/img/crackedParticle.png)
**Motivation**
This research will be the first investigations into both shear and unaxial loadings for AA-7075 for the largely anisotropic nature of the particles. This project would be the first to use a multistep approach with a low cost informatics framework for high cycle fatigues.
![Shear Lives](/MIC-AL7075-PARTICLES/img/shear_lives.png)
**Problem Statement/Objectives**
We wanted to be able to predict the life-limiting FIP for a given particle cluster, load amplitude, and load biaxiality. We will be using the Fatemi-Socie FIP to quantify fatigue lives, defined by $$FIP_{FS} = \frac{\gamma}{2} ( 1 + k\frac{\sigma}{\sigma_y} )$$. In this equation $\sigma$ is the tension on the maximum shear plane, $\sigma_y$ is the cyclic yield stress, $\frac{\gamma}{2}$ is the maximum plastic shear and k is a calibration constant for multi-axial fatigue (usually taken as 0.5). In this particular project, we have investigated the non-local averaged FIP surrounding the particles. An averaging volume of ~5% particle volume is used. The maximum of these FIP determines the life-limiting crack extension into the matrix.

**Data**
We have 6 microstructure images taken from multiple different samples. A 21x21x21 microstructure created from delta microstructures was used to train the MKS model for the purpose of the project. The six microstructure images were used in the creation of 10 200x200x200 microstructure reconstructions.
![An example of one of the microstructure images](/MIC-AL7075-PARTICLES/img/Presentation_Images/refined-4.png)
![Reconstruction](/MIC-AL7075-PARTICLES/img/Presentation_Images/3D_reconstruction.png)

**Workflow**
In general, our modeling process began with training our model using delta microstructures. We then tested our model under the different loading conditions and created quantifications of error to understand issues with the model. These error quantifications were used to improve the model before we predicted strain fields for the large microstructure reconstructions. Finally we used this information to preform the life predictions from FIPs.
![Workflow](/MIC-AL7075-PARTICLES/img/Presentation_Images/workflow.png)

**Project Challenges**
One of the first challenges we encountered was simply upon how to define distance measured to the nearest particle. Our original plan involved using Manhattan distance measurements, however our error quantifications found this to be a non-desirable way to measure distance. Therefore we ended up deciding on using radial distance for our measurements.
When we proceeded to the training section using delta microstructures in the MKS method, we initially used different Poisson ratios between the particle and matrix phases to a negative outcome on our predictions. After changing to using solely the Poisson ratio for the matrix, our results improved drastically. However we are still not entirely certain if this change was correct.
In order to have a demonstrative example of the difference between the reconstruction microstructure that was trained on the microstructure images and a reconstruction microstructure with randomized dispersal of particles, we created 2-point statistical images to compare the two. However in our first run, we found that while the 2-point statistics for the randomized microstructure was not capturing the stringers, as we would expect, the error quantifications between the two indicated that they were identical, which they were not. As a result, we changed from matrix auto-correlation to particle auto-correlation, which showed a slight improvement. From there we truncated based on the average stringer length and compared the two point statistics within that range. In doing this, we were able to see a massive difference between the reconstruction and the randomized reconstruction.
The MSE for the randomized microstructure is seen below
![Random](/MIC-AL7075-PARTICLES/img/Presentation_Images/MSE_random_stringers.png)
And the MSE for the reconstruction
![Reconstruction](/MIC-AL7075-PARTICLES/img/Presentation_Images/MSE_recon_stringer.png)

**Results**
Seen below are the fatigue predictions for both shear and uniaxial loadings. 
Several interesting features become readily apparent in the plots here. The first is that both loading conditions do indeed quickly approach a far-field condition. For this reason, both plots are restricted to the first 6 spatial bins. Another feature is that both loading conditions appear to produce a bimodal distribution. The uniaxial distribution's peak at the higher FIP value represents a larger relative intensification than that of the shear FIP distribution. The FIP values are still found to be higher in the shear loading condition, however, so while the intensification aspect is proved due to the anisotropic particle distribution, it is still insufficient to capture the discrepancy in our current fatigue model when compared to the literature results.
![Shear](/MIC-AL7075-PARTICLES/img/Presentation_Images/shear002.png)
![Uniaxial](/MIC-AL7075-PARTICLES/img/Presentation_Images/shear002.png)

**Conclusions**
Our work does seem to provide some evidence that the anisotropic distribution of particles contributes to the knock down in fatigue lives for uniaxial loadings.
There is definitely some physical phenomenon that our model is not fully capturing due to the assumptions that we made. In order to rectify this issue, we would need to move away from MKS, as we need to model the crack itself which is a contrast issue that MKS currently cannot capture.

**Acknowledgements**
Dr. Kalidindi
Dr. McDowell
David Brough
Yuksel Yabansu
NSF FLAMEL Program





