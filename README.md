## FENRIR Ray-Tracing Suite ##

Clean codebase for raytracing around thick accretion disks in AGN illuminated by a lamppost corona. Further instructions on how to use each file are in their headers. There are three sections, creating emissivity ray tracing, disk to observer ray tracing, and then a few helper functions and slurm scripts to combine the outputs and plot everything.  

Uses:

Corona Ray Tracing: /corona_new_clean

- Create emissivity plots: run jobs_clean/emissivity_profile.sb
- Create lamppost .fits files: run 

Disk to Observer Ray Tracing: /diskimg_clean

- Create redshift maps: run redshift_maps.sb
- Create transfer .fits files: run 

Combining the outputs: 

- Intensity Maps: run intensity_maps.sb
- Line profiles (one at a time): run line_prof.sb
- Create a large set of line profiles: line_prof_lots.sb

Analyzing the outputs: 

Other files/Information: 
- test_corona_params.txt and disk_params.txt hold the spin, height, and inclination parameters to loop through when doing sets of ray tracing.
- At the moment, all of the shell scripts are hardcoded with paths, but all of the python and cpp files should work if those shell scripts are edited.
- 
