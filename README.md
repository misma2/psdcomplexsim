# psdcomplexsim
 
 
 
 # PSD protein complex simulator

See online tool at: http://psdcomplexsim.cytocast.com/ calculates protein complexes of seven major postsynaptic proteins (two membrane receptors, NMDAR and AMPAR, as well as seven scaffold proteins) as a function of their individual abundance. This site is intended for demonstrative/test purposes, our results obtained using 524 experimentally obtained data sets is described in Miski et al.: Diversity of synaptic protein complexes as a function of the abundance of their constituent proteins: a modeling approach (submitted for publication). The abundance of the constituent proteins can be set using the form below. The (simplified) interactions considered between these proteins are shown below. By default, Homer proteins are modeled as dimers and Shank polymerization (via its SAM domain) is not considered, you can change these settings by ticking the appropriate boxes. The default length of the simulation run with this service is 1. The output contains the protein complexes formed along with their abundance.

Please refer to results obtained by this tool as:

Marcell Miski,Bence Márk Keömley-Horváth,Dorina Rákóczi Megyeriné, Attila Csikász-Nagy,Zoltán Gáspári; Diversity of synaptic protein complexes as a function of the abundance of their constituent proteins: a modeling approach (2021) Submitted to PLOS Computational Biology

# Main scripts

## run.sh
It runs the Cytocast application for several times for each brain region.
### for cycle j: 
sets how many times a simulation should be run with the same setup - this will be averaged during postprocessing
### for cycle i:
set the input data directory, runs through each different setup (brain region)
### /cytocastCoreNoOsm
is the command to start Cytocast see the Cytocast documantations for further information or -h for help

## generateInput/brainspan.py
Gets the preferred data from the mRNA data downloaded from: https://brainspan.org/static/download.html 
Generates files for SiComPre
## generateInput/create_input.py
Generates input json files for Cytocast (run.sh) from SiComPre abundance data.

## outputs/runPostProcessing.py
Runs the postprocessing calculation.
See the details and functions in the code - there is minimalist GUI
