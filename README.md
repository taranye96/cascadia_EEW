This repository contains the codes referenced in the following paper:

*Nye, T. A., V. J. Sahakian, A. Schlesinger, D. Melgar, A. Babaie Mahani, A. Williamson, E. Ferguson, A. Lux, and B. Pirenne (2026). Validation framework for semi-stochastic simulations in Cascadia earthquake early warning. Seismica, 5(2). doi: 10.26443/seismica.v5i1.1411*

**cascadia_EEW/figures** contains the scripts used to make the figures in the paper. 

**cascadia_EEW/src** contains the scripts used to generate the simulated data and perform subsequent analyses. 
- The scripts in this directory are organized numerically by the order in which they are run. For example, scripts starting with "0" should be run before those starting with "1". Scripts with a letter following the number should be run in order. 
- Scripts starting with "0" are used to run the simulations and analyze real noise from offshore sensors.
- Scripts starting with "1" are used to process the simulated waveforms and calculate intensity measures.
- Scripts starting with "2" are used to calculate P-arrivals and source-to-site azimuths and run the STA/LTA detection algorithms.
- The script starting with "3" is used to calculate the P-wave displacement amplitude.
- Scripts without a number are called within the other scripts.

The data for this project can be found (include Borealis reference once published)
