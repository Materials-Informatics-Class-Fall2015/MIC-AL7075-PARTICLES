---
layout:     	post
title:      	Poisson Ratio Effects
date:       	2015-10-31 23:00
author:     Paul Kern
tags:      result
---
<!-- Start Writing Below in Markdown -->
**Review**

Based on the previous discussion of MKS prediction errors and initial attempts to resolve this, the Poisson Ratio contrast was investigated. Additionally, instead of letting the sides contract naturally, the strain contraction was imposed in the loading.
MKS predictions were performed with simple superposition of the single local tensor component to the respective macroscopic tensor component.
The Poisson ratio for the particle and matrix phase are both now set as 0.3595

**Success!**
It appears that the contrast between the Poisson Ratio is detrimental to the predictions of the phase with the smaller volume fraction based on the below plot. We are now able to predict loading from a superposition without the discontinuity in the particle phase predictions.
![A](/MIC-AL7075-PARTICLES/img/Paul/ChangePoisson.png)

Additional simulations will be conducted to compare the predictions with free surfaces to isolate the effect of the imposition of loading on the contracting faces in predictions.
Unfortunately, it appears that the contrast between Poisson Ratio is not captured by the first order influence coefficients used in pyMKS. This means that we could be forced to choose between what effects we want to capture and how to accomplish it.
Compromising the material system representation further is required to accurately predict all strain values, it may be advisable to revert to the contrasting Poisson Ratios predictions. The reason is that the matrix values were accurately predicted and ultimately we want to quantify damage in the matrix.
Further discussion and final decision will be made in the coming months.

**Plan of Attack for November**
Ultimately our finalized product will be a trained series of models to predict multiaxial fatigue damage in dual phase matrix-particle systems for Al 7075-T6.
We will quantify error across the range of strains under consideration (0.2-0.4%)
Additionally we will quantify the damage intensification under specific loading conditions found in literature to verify discrepancies in fatigue lives.
In the next week we will determine which set of material system/MKS prediction methods we will be used moving forward.
Immediately following this, additional error quantification will be performed if necessary.
Verification of the reconstructed particle "stringers" will be conducted using the 2-point statistics discussed in class.
Predictions of the intensification under various loading conditions will finally be performed during the remainder of the class.