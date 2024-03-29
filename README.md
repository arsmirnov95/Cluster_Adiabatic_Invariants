# Adiabatic Invariants Calculations for Cluster Mission: codes.

### Codes for Phase Space Density (PSD) conversion of Cluster/RAPID/IES electron data.

The codes include:
* matlab script to compute the second invariant K
* Python script to compute Lstars with 1-min resolution
* set of Python functions to convert Cluster/RAPID/IES electron flux to PSD


The result product is available through Cluster science archive ('LSTAR' product). 

##### When using the 'LSTAR' product please add the DOI for this code to references: 10.5281/zenodo.3519999.


For computations, one needs the information from cdf files taken from Cluster science
archive and Kp index from OMNI data (with 1-hour resolution).

Generated by A.G. Smirnov, E.A. Kronberg, P.W. Daly, N.A. Aseev, Y.Y. Shprits and A. Kellerman.

In case of any questions, contact Artem Smirnov via arsmirnov95@gmail.com.


To run the code, it is necessary to have IRBEM library installed on the computer.
Other Python dependencies:
* numpy 
* scipy
* sys 
* datetime
* spacepy
* pandas
