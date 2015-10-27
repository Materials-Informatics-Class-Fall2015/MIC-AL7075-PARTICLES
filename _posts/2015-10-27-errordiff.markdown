---
layout:     	post
title:      	Error Magnitude between Finite Element and PyMKS
date:       	2015-10-27 13:36
author:     Chris Shartrand
tags:         result
---
<!-- Start Writing Below in Markdown -->
**Review**
Recall from our previous post that after failing to observe the predicted phenomenon with the binning for Manhattan distance, it was decided to purely look into the magnitude of error between the Finite Element Method and the PyMKS method.
**Loading Conditions**
For our project, we have been working with both uniaxial and shear loading conditions in our investigation of our microstructures strain values. For uniaxial, we set a baseline strain of .002 and conducted two other computations at a 50% and 200% increase in the baseline value. Under shear, the set baseline was .004 with two other computations conducted at a 25% and 50% increase.
Under uniaxial with a 200% increase in the baseline strain, the following selection of error magnitude was found:
![A](/MIC-AL7075-PARTICLES/img/Presentation_Images/uniaxial200percent.png)
Where clearly these orders of magnitude are quite large.
Now when considering shear with a 50% increase in the baseline strain, the following selection of error magnitude was found:
![A](/MIC-AL7075-PARTICLES/img/Presentation_Images/shearbaseline.png)
Whose error is much smaller than that of the uniaxial 200% increase.
Further investigation yielded clear results that showed the order of magnitude to be small under uniaxial baseline and 50% and shear baseline, 25% and 50%. Therefore we concluded that the uniaxial 200% would be outside the scope of our model's predictive capabilities. 