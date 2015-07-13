# 3 site code for Anderson-Hubbard systems
<html>
<body>
<h2>Summary</h2>
<p> This code finds the density of states (DOS) and the general inverse participation (GIPR) for a three site system for half filling. The program makes the hamiltonian matrices for each combination of up and down electrons and diagonalizes them using LAPACK (linear algebra library). It compares the eigenenergies to find the lowest grand potential (the ground state). 
</p>
<p>The program using a pre made lookup table (found in routines.f90) to find the states that the ground state vector can photo emmit or inverse photo emmit too. These tables are called IPES and PES (inverse photo emmision and photo emmision spectrum). </p>
<p>It makes contributions to the DOS at the frequency (energy) equal to the difference of that state to the ground grand potential. This calculates the local density of states (LDOS). The LDOS for each of the sites is averaged and added to the total DOS. The DOS is made smooth by using energy (frequency) bins.  </p>
</body>
</html>
