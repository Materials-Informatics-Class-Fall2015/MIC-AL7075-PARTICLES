---
layout: post
title: Proposed Methodology
date: 2015-09-09 01:30
author: Paul Kern
tags: background information
---
<!-- Start Writing Below in Markdown -->
**The Structure**
Al 7075-T6 contains many very stiff particles (2-3 times matrix stiffness) whose behavior bounds the fatigue lives [1]. These particles tend to cluster and form "stringers" along the rolling direction [1].
Rollett et al. have successfully used pair correlation functions to describe the spatial distribution of particles within the matrix [2]. Pair correlation functions (PCF) describe the likelihood of finding a like-phase (in this case particle) at a given distance and radial direction from the current point. This inherently 2D (or higher dimensional) function can be very useful to quantify microstructure.

![A sample image of the PCF from the rolling plane](/MIC-AL7075-PARTICLES/img/pcf.png) [2]

We hope to use the Representative Volume Element (RVE) of particles matching the experimental PCF used by Rollett et al. to instantiate several smaller Statistical Volume Elements (SVE) containing clusters of 3-5 particles. We then will process the SVEs and compute the PCF for each SVE in all planes.

**The Property**
Fatigue Indicator Parameters (FIP) are used to quantify damage, especially in metals. The most common FIP used in the McDowell lab to quantify fatigue lives is known as the Fatemi-Socie FIP and is defined as follows
$$FIP_{FS} = \frac{\gamma}{2} ( 1 + k\frac{\sigma}{\sigma_y} )$$
In this equation $\sigma$ is the tension on the maximum shear plane, $\sigma_y$ is the cyclic yield stress, $\frac{\gamma}{2}$ is the maximum plastic shear and k is a calibration constant for multi-axial fatigue (usually taken as 0.5) [3].
Since the FIP depends upon the plastic deformation, at the minimum an elastic-plastic constitutive model must be defined. In this case we are using a linear-elastic model for the stiff particles with a Young's Modulus of 169 GPa and Poisson Ratio of 0.3
The elastic-plastic model defined for the aluminum matrix has a Young's Modulus of 69 GPa and Poisson Ratio of 0.345. The combined kinematic-isotropic hardening model uses the stresses and plastic strains from a calibrated crystal plasticity model.

In this particular project, we will be investigating the non-local averaged FIP surrounding the particles. An averaging volume of ~5% particle volume will be used. The maximum of these FIP will determine the life-limiting crack extension into the matrix.

**The Linkages**
Ideally, we want to be able to predict the life-limiting FIP for a given particle cluster, load amplitude, and load biaxiality. As previously stated, the variation in the particle clusters will be due to taking SVE from the larger RVE.
The load amplitude and biaxiality will be varied parametrically within the region of interest ([0.2,0.6]% strain and [-1,1] respectively)
To create the linkage, lowering of the PCF dimensionality will be necessary via Principal Component Analysis (PCA). Once this is accomplished a Neural Network will be trained on the inputs and outputs to generate a meta-model applicable for predictions of life-limiting clusters.

[1] Xue, Et al, 2007, “Micromechanisms of multistage fatigue crack growth in a high-strength aluminum alloy”, Acta Materialia, Vol. 55, pp 1975-1984.
[2] Rollett Et al, 2006, “Three Dimensional Microstructures: Statistical Analysis of Second Phase Particles in AA7075-T651”, Materials Science Forum, Vols. 519-521, pp 1-10.
[3] Fatemi, Ali, and Darrell F. Socie. "A Critical Plane Approach to Multiaxial Fatigue Damage Including out‐of‐Phase Loading." Fatigue & Fracture of Engineering Materials & Structures 11.3 (1988): 149-165.


