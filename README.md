# Imaginary-time-Time-Dependent-DFT-for-Quantum-ESPRESSO v7.1
A modification to pw.x of the DFT package Quantum ESPRESSO (QE) v7.1 that enables the user to imaginary-time propagate the wave function.  https://iopscience.iop.org/article/10.1088/1361-648X/abbe7e
For Linux.

Imaginary-time propagation is an alternative means of producing the Kohn-Sham ground state or delta-SCF excited states. Compared to the standard SCF method, the main benefit of imaginary-time propagation is that it usually converges in cases where SCF does not. However, it is slower than SCF.

This modification can automatically switch from SCF to imaginary-time propagation when SCF is unable to converge.

This modification is GPU accelerated, however, combining GPU acceleration with ultra-soft pseudo potentials my be very slow because the inverse S-matrix used by the propagator is not GPU accelerated. Gamma only and multiple k-point calculations are possible. DFT+U calculations do not seem to work.

To setup:

  -first install QE v7.1  https://www.quantum-espresso.org/quantumespresso-v-7-1/ or use an existing version
  
  -download the modified files in this repository then place the file input_parameters.f90 in the folder MODULES and the rest of the files in the folder PW/src
  
  -open a terminal in the main folder and type './configure'  then 'make all'
  

Input variables for this modification are part of the "system" name list of the QE input file. Descriptions of the input variables are as follows:

run_ITDDFT
   A logical value that will run imaginary time propagation if .true. and will run standard QE if .false. 
   Defalt: .false.
   
      
g_dependant_dtau

   A logical value. If true, g-vectors that have less kinetic energy will be propagated with a greater time step. This results in faster convergence.
   Default: .false.
   
dtau

   The time step length is equal to dtau multiplied by the inverse of the maximum g-vector kinetic energy.  Except when g_dependant_dtau=.true., where The time  step length is less than or equal to dtau multiplied by the inverse of the g-vector kinetic energy.
   Default: 1   
   
ITDDFT_backup 

   A logical value. If true, pw.x will switch to imaginary time propagation automatically if SCF fails to converge. run_ITDDFT must be set to false.

Default: .false.


It is recommended that the user use the following settings:

ITDDFT_backup = .true.
g_dependant_dtau = .true.
dtau = 1
run_ITDDFT = .false.
