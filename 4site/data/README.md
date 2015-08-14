<html>

<h2>Data Naming Conventions</h2>
<p>All the data files begin with "<em>4-</em>" indicating it is an ensemble of 4 site systems. This is followed by "<em>dos+ipr_</em>" signifying that it contains the density of states and the generalized inverse participation ratio. The rest of the filename is assigned based on the parameters of the simulation that created it. The parameter are in the order t,U,W,mu  where t is the hopping integral, U is the on-site interactions, W is the width of the disorder potential and mu is the chemical potential. The name of the variable comes first and is immediately followed by it's value.</p>
<p><strong>Example:</strong> <em>4-dos+ipr_t-1U4W8mu4.0.dat</em> is for a simulation run with t=-1, U=4, W=8 and mu=4.0.</p>
</html>
