# Control_Theory repo

Since I have a lot of free time, in this repo I will sum up some aspects of Control Theory as I studied it at the electrical engineering school of NTUA.

It ain't much but it's honest work.

## 1) root_locus_calculation.m and root_locus.m

Both functions use as input a matlab continuous-time transfer function tf.

* *root_locus_calculation.m*: This function uses matlab's rlocus to find the trajectories and respective gains. Then it calculates either analytically or numerically the basic characteristics of the rlocus: Departure and Arrival angles from the poles and zeros respectively (γωνίες αναχώρησης-αφιξης), breakaway/in points (σημεία θλάσης) , asymptotes and approximate intersection points with the y'y axis (with the respective gain K of borderline stability of the closed loop).

* *root_locus.m*: It uses *root_locus_calculation.m* to find the characteristics of the root locus and plots them nicely.

The example.m shows a use of the latter functions, enjoy :) 


Simple rlocus: | root_locus result:
------------ | -------------
![rlocus_image](/Example_Images/rlocus.png) | ![root_locus_image](/Example_Images/root_locus.png) 



