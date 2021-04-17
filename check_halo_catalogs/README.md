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

