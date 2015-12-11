---
layout:     	post
title:      	Last Thoughts
date:       	2015-12-11 18:00
author:     Paul Kern
tags:         
---
<!-- Start Writing Below in Markdown -->

**Model Integration**

Some questions were asked about how we plan to use this moving forward. The current model in the McDowell group works by calculating local (simulation element) fatigue response and propagating a crack radially (and on varying slip planes).
The first step in the process, however, is actually "cracking" a particle based on the highest local plasticity response. We can use our new FIP responses to provide the intensification factor assuming that the far-field FIPs that we get from MKS are proportional to the local response in the homogenized model.
Since it's hard (computationally inefficient) to model these very small particles at the same time as we're modelling the local grain structure, we can instead overlay virtual particles simulated with MKS and use the intensification responses from these cheaper predictions.
Furthermore, during crack propagation we can use these virtual particles and crack propagation to influence the local crack growth rate as the crack tip interacts with areas of particle driven intensification. This part will require some more studies with high fidelity simulation to train some interaction factors between the crack tip and particle distributions.

**Further MKS Improvements**

The MKS method itself can be improved to incorporate crystal orientation and higher order terms to better capture the actual plasticity.
Local stress responses can then be driven by the appropriate anisotropic stiffness tensor to recover the local responses more accurately.
This approach will still be limited to fully reversed simulations, but will be a much better approximation than we are currently achieving.
Noah Paulson has done a lot of great work on both of these things with a Ti64 material system, so we will be reaching out to him to explore opportunities for improvements in the MKS prediction side of things.
We also want to investigate more fundamentally what happens with the differing Poisson ratios in the MKS predictions. We understand that the changing them both to be the same improved the prediction, but are uncertain of why this introduced such strong bias in the first place.