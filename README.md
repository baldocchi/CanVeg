# CanVeg
CANVEG is a one-dimensional, multi-layer biosphere-atmosphere gas exchange model to compute water vapor, CO2 and isoprene flux densities. The model consists of coupled micrometeorological and eco-physiological modules.  The micrometeorological modules compute leaf and soil energy exchange, turbulent diffusion, scalar concentration profiles and radiative transfer through the canopy.  Environmental variables, computed with the micrometeorological module, in turn, drive the physiological modules that compute leaf photosynthesis, stomatal conductance, transpiration and leaf, bole and soil/root respiration, as well as isoprene emission rates. The micrometeorological turbulence diffusion is based on Lagrangian theory that produces a dispersion matrix.

Canveg.m is the main program. It calls specific parameter subroutines for the canopy being studied. Here we call input met files, set up latitude and longitude, set canopy, vegetaion, and soil conditions, etc

Details about the model are in the file Canveg_documentation
