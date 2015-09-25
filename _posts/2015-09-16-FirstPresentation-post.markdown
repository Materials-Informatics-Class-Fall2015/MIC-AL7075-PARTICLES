---
layout:     	post
title:     		Copy of Presentation for Comments
date:      	2015-09-16 03:00 
author:     	Paul Kern and Chris Shartrand

theme:		moon # default/beige/blood/moon/night/serif/simple/sky/solarized
trans:		default # default/cube/page/concave/zoom/linear/fade/none

horizontal:	</section></section><section markdown="1" data-background="http://matin-hub.github.io/project-pages/img/slidebackground.png"><section markdown="1">
vertical:		</section><section markdown="1">
---
<section markdown="1" data-background="http://matin-hub.github.io/project-pages/img/slidebackground.png"><section markdown="1">
## {{ page.title }}

<hr>

#### {{ page.author }}

#### {{ page.date | | date: "%I %M %p ,%a, %b %d %Y"}}

<!-- Start Writing Below in Markdown -->

## Background Knowledge ##

 - Aluminum Alloy consisting of Zinc,  Magnesium, Copper, Iron, Silicon, Manganese, Chromium, Zirconium, and Titanium.
 - During processing it is common for μm-sized particles to become ingrained in the matrix
 - Deformation can cause particles to become cracked or detached
 - Cracks combined with stress can lead to failure in the alloy


 ![A Cracked Particle](/MIC-AL7075-PARTICLES/img/crackedParticle.png)

## Importance ##

 - AA7075 is commonly used in the construction of airframes
 - Product failure = Costly
 - Improved modeling can lead to prediction of microstructures that have a higher probability of cracking.
 - Quantification of product quality.
 

## The Microstructure ##
 - Particles are very stiff
 - Tend to cluster and form "stringers" along the direction of the roll
 - Use Pair Correlation Functions (PCF) to quantify probability of finding another cluster from a given point
![A sample image of the PCF from the rolling plane](/MIC-AL7075-PARTICLES/img/pcf.png)


## Microstructure Property ##

 - Fatigue Indicator Parameters (FIP) are used to quantify damage
 - We will use the Fatemi-Socie FIP: $FIP_{FS} = \frac{\gamma}{2} ( 1 + k\frac{\sigma}{\sigma_y} )$
 - An elastic-plastic constitutive model must be defined
 - We are using a linear-elastic model for the stiff particles with: Young’s Modulus of 169 GPa and Poisson Ratio of 0.3
 - Investigating the non-local averaged FIP surrounding the particles
 - Averaging volume of ~5% particle volume will be used
  
## Simulation ##

![3D-Simulation of Particles](/MIC-AL7075-PARTICLES/img/Presentation_Images/Particles_D3D.png)

## Two-Point Statistics ##
![Two Point Statistics](/MIC-AL7075-PARTICLES/img/Presentation_Images/2_point.png)

<!-- End Here -->

#[Print]({{ site.url }}{{ site.baseurl }}{{ page.url }}/?print-pdf#)

#[Back]({{ site.url }}{{ site.baseurl }})

</section></section>
