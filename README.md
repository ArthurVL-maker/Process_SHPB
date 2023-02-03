Process_SHPB, an open-source algorithm for stress wave dispersion correction in split-Hopkinson pressure bar experiments

When processing signals from split-Hopkinson pressure bar (SHPB) experiments, it is frequently presumed that longitudinal stress waves in the pressure bars travel one-dimensionally at a common velocity of c0. Hence, measurements recorded at the strain gauges are commonly simply translated to the end of the bar using a suitable time delay. In reality, stress waves travel at a certain phase velocity, cp, which varies with frequency, bar diameter, one-dimensional wave speed and Poisson’s ratio. As the frequency of a wave rises, the phase velocity drops, resulting in signal dispersion as it propagates down the bar. The dispersion of the stress pulse is followed by a frequency-dependent fluctuation in stress and strain throughout the bar cross-section. Therefore, a signal recorded on the surface of the bar at some distance from the specimen will not accurately represent the stresses the specimen was subjected to, and hence cannot be used to objectively measure the specimen response.

Process_SHPB includes shpb_main.py, the main algorithm required to process all the SHPB test data, from which the subroutines dispersion.py and dispersion_factors.py are called. It was inspired by a MATLAB script created by Barr (2016).

The inputs for the key function named shpb_main are: raw file, which is the raw file collected from the SHPB test; sample data, which is a list containing the length, mass and dry mass of the sample tested; confinement, which specifies the confinement of the sample tested using the SHPB apparatus either ‘”None”, “Ring” or “Reservoir”; and correction, which is either a dispersion correction or a simple time shift analysis. In the Processed_Data folder, all the necessary key results are saved as .csv files.

The subroutine dispersion.py, in the main function shpb_main.py, uses an adaptation of Tyas and Pope’s (2005) dispersion-correction approach to verify that the inferred axial stress and strain data appropriately represent the specimen behaviour, as specified below:

1.	The fast Fourier transform (FFT) is used to convert the time-domain strain signal to the frequency domain.
2.	Using Bancroft’s equation, a correction is made to the phase angle of each frequency component to account for dispersion over the distance between the strain gauge and the bar end. 
3.	The amplitude of each frequency component is corrected using the factors m1, m2, v_ratio and norm_freq, which account for strain and Young’s modulus fluctuation over the bar cross section, respectively. Davies’ investigation of radial effects in a cylindrical pressure bar yielded these results.
4.	The inverse FFT is used to transform the signal back into the time domain.

To save time, dispersion.py utilises a subroutine dispersion_factors.py, which includes a precalculated, normalised look-up table of phase velocity, m1, m2, v_ratio and norm_freq. Process_SHPB, includes a lookup table for a Poisson’s ratio of 0.29. More tables can be generated using the calculation method outlined in Tyas and Pope (2005).
