# 3D Single Molecule Localisation Microscopy Data Analysis 
## applied to 3D PSD95 localisations


<img src=Sample.png width="600" height="300">

AnalyseExperiment:
- separates the localisations in a monomer and nanocluster component
- based on envelope tests checks if the momers are spatially random (CSR) or better modelled by a Thomas process. Similarly, if nanoclusters are better described by a Thomas or double Thomas process
- computes the parameters of the above models, based on pair correlation function fits.
The R scripts allows the simulation of 3D Thomas and double Thomas point processes.
