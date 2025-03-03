# EEG-DOT
Codes used for paper "Diffuse optical tomography spatial prior for EEG source localization in human visual cortex"

Note that these are NOT the codes used for our previous simulation work "Enhanced spatiotemporal resolution imaging of neuronal activity using joint electroencephalography and diffuse optical tomography". The codes for that paper can be found [here](https://github.com/JiamingCao/NIRS-EEG).

However, we do recommend using the codes here, as they are more up-to-date.

If you use these codes, please make sure to also cite our paper:
- Cao, Jiaming, et al. "Diffuse optical tomography spatial prior for EEG source localization in human visual cortex." NeuroImage 277 (2023): 120210. https://doi.org/10.1016/j.neuroimage.2023.120210

## How to use
- The routine starts with ```ForwardModels.m```, which does exactly what the name suggests. The differences from the codes you have are 1) EEG forward model is calculated using BEM, and 2) the cortical source locations are assumed to be the pial surface (so the brain looks more "folded" in the new results)
 
- Then ```processEEG_all.m``` and ```processNIRS_reconall.m```. They loop through all the subjects and run all the preprocessing, and the core functions that do all the heavy lifting are ```func_eegblocks.m``` and ```func_nirsrecon.m```.

   Quite a bit of preprocessing was involved (also reported in the paper), and you may wish to add, remove, or change the parameters of some steps, depending on your own data. Reconstruction of NIRS data using classical Tikhonov regularization is also performed in this stage.
 
- After processing each individual, the group average is done using ```GroupEEG_leftright.m``` and ```GroupAverage_recon_leftright.m```.
 
- Finally, you will need to run ```EEG_NIRS.m```, which runs the ReML-based fusion at each EEG timestep. You may wish to change the parameters here and there according to your data, although from my experience, the difference may not be too dramatic.

## Some notes
- You will need the following extra toolboxes
  - [FieldTrip](https://github.com/fieldtrip/fieldtrip)
  - [NIRFASTer](https://github.com/nirfaster/NIRFASTer)
  - [iso2mesh](https://github.com/fangq/iso2mesh)
- Some of the supporting functions are adapted from the [nirs-toolbox](https://github.com/huppertt/nirs-toolbox), which is marked in the preamble parts of the files. The ```iso2mesh_plotmesh``` function is, as the name suggests, the ```plotmesh``` function that shipped with the [iso2mesh](https://github.com/fangq/iso2mesh) toolbox. It is renamed simply to avoid confusion with a NIRFASTer function with the same name
- The functions adapted from the other toolboxes were based on the previous versions. Therefore, it is possible that they do not resemble or match the performance of the latest versions of the corresponding functions anymore.

