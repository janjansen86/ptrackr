# ptrackr

Both tracking in 2D and 3D work well. 

For a large number of points, looping over short trackit-functions and removing the particles that have settled makes things a lot faster. That's the purpose of the function loopit_2D3D.

## other things to do: 
 - Merge the two trackit functions into a single function as most of the code is redundant (the respective lines have comments).
 - write a function for mapping the output of loopit
 - finish vignette

