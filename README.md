# Egocentric Value maps
This contains the scripts necessary to reproduce the simulations and analyses for the Egocentric Value maps manuscript

*For any questions, you can contact me rory.bufacchi@gmail.com*

## System Requirements
All plots and analyses were performed on Matlab 2022b, but any version of Matlab from 2022a onward should do.

The analyses require only a standard computer with enough RAM to support the in-memory operations.

## Installation
Just unzip the EgocentricValueMaps folder wherever you find appropriate and run
```
addpath(genpath('chosen_filepath\EgocentricValueMaps'));
```

## Usage Examples

### Demo
Demo.m shows how to generate and visualise bodypart-centred Q-fields using tabular learning (SECTION 1), and simple neural networks (SECTION 2).

### Run analyses for figures and stats
The scripts necessary to create each figure and the corresponding stats are in the 'Figures' folder. If you donwload the precomputed data from https://doi.org/10.5281/zenodo.16408688, you will be able to run all the "*Figures*" .m files directly. Just change the line
```
dataPath                    = '';
```
at the top of the .m file to be whatever folder you have unzipper the 'Results' folder into. Alternatively, select the appropriate folder when prompted by the pop-up. 

### Generate data
To generate the simulated data for both ANN and tabular simulations, run (sections of) *CreateDataForPlots_Final.m*

*! ----------------WARNING:-----------------!*

This will take VERY long to run. I *highly recommend* downloading the precomputed data from https://doi.org/10.5281/zenodo.16408688. If you want to follow the full analyses step by step, I  recommended running the sub-sections of *CreateDataForPlots_Final.m* individually, rather than just running the whole .m file as that willtake forever.


### MATLAB Workspaces
The necessary MATLAB workspaces for this project can be downloaded from the following link:

https://doi.org/10.5281/zenodo.16408688

- For recreating figures and stats, the appropriate workspace for each Figure is loaded at the top of the corresponding .m file.
- In CreateDataForPlots_Final.m, the appropriate workspace is loaded at the top of each sub-section of the .m file.
