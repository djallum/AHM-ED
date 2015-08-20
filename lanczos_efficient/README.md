<html>
<body>
<h2>Lanczos Efficient</h2>
<p> This program was orginally written by Dr. Atkinson. It has since been slightly modified and comments have been added with the end goal of using it for Anderson-Hubbard model simulations with large system sizes due to it's ability to solve the many-body eigenstates and eigenvalues of large systems quickly (up to 16 sites). <p>
<h3>General Structure</h3>
<p>The program is split into two segments: Zgroundstate and Zspectrum. The first find the lowest many-body eigenstate of a particular Hamiltonian submatrix and the second uses this information to contruct the denstiy of states (DOS).</p>
<p>The second half of the program has not been properly investigated however there is a flow chart for the first program (zgroundstate.html) which can be opened at this website: https://www.draw.io as well as considerable commenting of the code.</p>
<h3> Current Status</h3>
<p> The program is currently able to solve the Hamiltonian submatrices of fock states with  0 < n_up < nsites  0 < n_dn < nsites. For example for a 4 site system it can solve all the matrices that don't have either 0 or 4 up electrons or 0 or 4 down electrons.</p>
<p> A script (<em>main.sh</em>) has been written which will run the Zgroundstate portion for all the possible submatrices (except ones that currently have issues as mentioned above) and compares their values to determine the amount of electrons in the groundstate. It then sends this information to the second portion (Zspectrum) for it to solve the DOS in terms of the proper ground state.<p>
<p> The first program needs to be modified in order to solve all possible Hamiltonian submatrices.<p>
<p> An additonal error has also found, when the site potentials are all the same the lowest grand potential is not accurate. If the values are change slightly (1,1,1,1 to 1,1,1,1.01) the lowest grand potential changes by a lot and returns to the proper value. A simple fix to this problem would be a <em>if</em> statement that could test if the error would occur and if so slightly modify one of the site potentials.</p> 
</body>
</html>
