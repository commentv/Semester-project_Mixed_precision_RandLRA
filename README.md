# Presentation
This project aims at studying the impact of mixed precision in the approximation quality of different randomized low rank approximation algorithms
such as the randomized nyström approximation or the randomized SVD. The code provided can be used to reproduce the figures in ```Report_semester_project.pdf```
or derive other figures. The code only implements the Nyström approximation using 2 different methods to derive the pseudoinverse:
The cholesky factorisation [1] or the epsilon pseudoinverse [2].

# Explanation
The code of the project is organized as follows :

-```/src/create_example``` is a matlab function file to create the examples used in Figures 1 to 5.

-```/src/Nystom```, ```/src/Nystom_single```, ```/src/Nystom_eps_pinv``` and ```/src/Nystom_eps_pinv_single``` are matlab function files used to compute
the nystrom approximation of a given matrix, using the different ways to produce the pseudoinverse, using mixed precision or by having all the steps in single precision.

-```/src/chop-master``` is a folder which contains all necessary tools to simulate the impact of the use of half precision arithmetic. For more details, see [3].

-```/src/Script_nystrom_chol```, ```/src/Script_nystrom_chol_all_single```, ```/src/Script_nystrom_eps_pinv``` and ```/src/Script_nystrom_eps_pinv_all_single``` are matlab scripts used to
plot the absolute error of the approximation when the nystrom approximation is computed with all steps in single precision or with mixed precision, depending on how we compute the pseudoinverse.
Those scripts are used to produce Figures 1 to 5.

-```/src/Script_nystrom_eps_pinv_all_single``` is a matlab script used to compare the absolute error when the approximation is computed by the nystrom algorithm with all the steps in single
precision for the cholesy factorisation and the epsilon pseudoinverse. This script is used to produce Figures 4 and 5.

-```/src/Script_abalone``` and ```/src/Script_uniform``` is a matlab script that produces Figures 6 and 7 when both the mixed precision nystrom algorithms are applied to the RBFK matrix derived 
from the abalone data or a uniform distribution. Note that the abalone data is located in the ```/data``` folder. To get more information on the abalone data, see [4].

Note that the files used to produce the different results were executed on Matlab R2022b and do not need any special libraries.

# References
[1] Erin Carson and Ieva Dauzickait ˇ e. Single-pass nystr ˙ om approximation in mixed precision. 05 2022.<br /><br />

[2] Yuji Nakatsukasa. Fast and stable randomized low-rank matrix approximation, 2020.<br /><br />

[3] Nicholas Higham and Srikara Pranesh. Simulating low precision floating-point arithmetic. SIAM Journal
on Scientific Computing, 41(5), 2019.<br /><br />

[4] A. Asuncion and D.J. Newman. UCI Machine Learning Repository. ´ http://www.ics.uci.edu/
˜mlearn/MLRepository.html, 2007.<br /><br />
