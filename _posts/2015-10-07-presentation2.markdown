---
layout:     	post
title:      	Presentation 2 for Comments
date:       	2015-10-07 19:23
author:     	Paul Kern and Chris Shartrand
tags:         result
---
<!-- Start Writing Below in Markdown -->
## Presentation Overview ##
 - Transition of project plan​
 - Localization with MKS​
 - Data formats​
 - Workflow​
 - Preliminary Analysis
![](/MIC-AL7075-PARTICLES/tree/gh-pages/img/Presentation_Images/Pres2_Img3)
## Where We've Been and Where We are Going ##
Goal still to use meta-model to predict life damaging nature of particle clusters in multiaxial fatigue. We are not using techniques talked about so far in class. We are using the localization technique briefly introduced in the workshop where we quantify uncertainties associated with MKS prediction of multiaxial fatigue.
We are interested in stability of predictions and measures of dispersion between small test microstructures and MKS predictions.

 - Microstructure – array(#loads,#ms,x,y,z)
 - Strains – array(#loads,6,#ms,x,y,z)
 - Models – array(#loads,6)
 - Base Folder
 -DoE.csv mapping folder # to load conditions
 - Elem_grain_#.txt
 - Strain_#.txt
 ![](/MIC-AL7075-PARTICLES/tree/gh-pages/img/Presentation_Images/Pres2_Img1)

Training
 - Simulate delta microstructures (2-phase so same as example online), Predict full tensor necessary because of the multiaxial nature, Store
   influence coefficients for later use
Testing
 - Same size MS, Compare to FEM, Currently test at same loads
Predicting
 - Slicing of larger microstructure, Expand coefficients, Use strains to
   predict local plasticity and recover FIP
![](/MIC-AL7075-PARTICLES/tree/gh-pages/img/Presentation_Images/Pres2_Img2)
## Source ##
https://commons.wikimedia.org/wiki/File:Stress_Strain_Ductile_Material.pdf