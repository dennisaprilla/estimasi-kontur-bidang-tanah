# Introduction

 This repository contains a quick implementation of the Direction Cosine Matrix (DCM) algorithm in Matlab found in [1]. Two functions are provided `reset_fusion` and `dcm_algorithm`. The first one is used for initializing the DCM matrix, and the second one is used to update it each time interval. In the `test` folder you can see an example on how to use it.

# Hardware setup
TODO:

# Testing
The Matlab implementation has been tested and verified against the C implementation in [1], using an Arduino with a SparkFun 9 DoF Sensor Stick [2].
 
 
<p align="center">
  <img src="https://github.com/alrevuelta/dcm-implementation/blob/master/img/matlab_vs_c.png">
</p>

# References
* [1] https://github.com/Razor-AHRS/razor-9dof-ahrs/tree/master/Arduino/Razor_AHRS
* [2] https://www.sparkfun.com/products/retired/10724
* [3] http://www.mdpi.com/1424-8220/15/3/7016/pdf
* [4] http://www.ti.com/lit/an/slaa518a/slaa518a.pdf
* [5] http://www.starlino.com/dcm_tutorial.html
* [6] http://www.academia.edu/11778144/Direction_Cosine_Matrix_IMU_Theory
