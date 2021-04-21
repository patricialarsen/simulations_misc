## Check halo catalogs 

This inputs quantities from a halo catalog output and makes histograms of distributions. These are output as txt files in the outputs folder. 

### Implemented Tests 

- Position histograms
- Velocity histograms
- Spin/ angular momentum histograms
- Mass/ angular momentum 2d histogram 


### To implement

- Triaxiality histogram
- SOD checks 
- vmax checks 

### How to run
Source the hacc environment scripts, edit the makefile to link to the latest version of trunk then type make, the syntax for running this is 

```
mpirun -np 1 ./make_hist $FILENAME $OUTPUT_FOLDER
```

You can then plot these histograms using 

```
python plot_outputs.py $OUTPUT_FOLDER 
```

This will plot to the display by default. If you want to suppress this, set display=False in the code. In both cases it will save the plots to the output folder. 

### Last Journey tests 

We may want to implement some of these here. They were

- compare M200 values to the Bolshoi Planck ones

checked that minimum fof count is as expected (20 particles)
checked that maximum fof mass sounds sensible (10^14.5)
checked that fof mass and fof count match (with particle mass of 2.79x10^9)
checked that fof_count distribution looks sensible (plot in folder)
checked uniqueness of fof tags (3936 duplicates, ~0.07%)
checked mass distribution of duplicate halos (plot in folder)
checked that fof and sod masses roughly line up (plot in folder - note these are log10 scaled)

checked ke distribution is similar between high mass fof and sod values (plot in folder)
checked ke distribution looks like a peaked distribution (plot in folder)
checked ke is linked to fof mass
checked ke orders of magnitude seem roughly right compared to velocity dispersions

checked that all sod masses for fof_count<500 are -101 (fill value)
checked that no sod masses for fof_count>=500 have fill values

checked that center positions are distributed sensibly
checked that mean and center positions roughly line up (maximum distance is roughly ~1 Mpc)
sod_halo_min_pot_x etc. are exact duplicates of the center positions (note redundancy)
sod_halo_mean_x etc. roughly line up with mean positions

ang_momentum has similar distributions in all directions (sod and fof) - these values are very high
sod values are biased lower than fof values
checking that angular momentum goes as M^5/3 - virial equilibrium for massive halos
- not checking overall amplitude but the mass scaling is correct.
- checking values against bolshoi planck - assuming *a^2 is the correct conversion from comoving to physical we agree pretty well.


velocity dispersion information:
note eq. 1 in https://arxiv.org/pdf/astro-ph/0702241.pdf - 1d velocity dispersion.
Specific thermal energy in dark matter halo should scale with GM/R, R goes as M1/3 and kinetic energy and M^2/3.
- relation with mass has roughly right scatter.


