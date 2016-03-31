# ptrackr

3D-tracking of particles from the surface to the seafloor works. The basic function for that is "trackit_3D". Examples work as well.

2D-tracking of particles with an option to stop particles when floor currents are low is NOT WORKING yet.
The basic function for that is "trackit_2D". 
It's based on "trackit_3D" with the additional lines marked quite obvious.

For a large number of points, looping over short run trackit-functions makes things a lot faster.
Thats the aim of "loopit_2D3D" which is NOT WORKING yet for 2D, only for 3D.

I wasn't able to make "loopit_2D3D" run with the basic "trackit" functions. 
So I altered them slightly and named those functions called "loopittrackit2D" and "loopittrackit_3D".
Again, the 2D version is not working yet, and I've altered the basic "trackit_2D" function since I created the loopit-version, so I'd need to recreate the loopit-version once the basic function works.



## other things to do: 
 - write a function to map the output of loopit

