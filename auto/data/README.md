<html>

<h2>Data Naming Conventions</h2>
<h3>Basic</h3>
<p>All the data files begin with "<em>a</em>" indicating it was created by the auto code. This is followed by a number to indicate the size of the systems (ex. 4). This is followed by "<em>dos+ipr_</em>" signifying that it contains the density of states and the generalized inverse participation ratio. The rest of the filename is assigned based on the parameters of the simulation that created it. The parameter are in the order t,U,W,mu  where t is the hopping integral, U is the on-site interactions, W is the width of the disorder potential and mu is the chemical potential. The name of the variable comes first and is immediately followed by it's value.</p>
<p><strong>Example:</strong> <em>a4-dos+ipr_t-1U4W8mu4.0.dat</em> is for a ensemble of 4 site systems run with t=-1, U=4, W=8 and mu=4.0.</p>
<h3>Additional Features</h3>
<p>To avoid data being over written if a file with a certain name already exists the data file will have a "<em>_1</em>" add the end. The number "1" can be replaced by any number less then 10. This allows 10 different data files for a simulations run with identical parameters.</p>
<p><strong>Example:</strong> <em>a4-dos+ipr_t-1U4W8mu4.0&#95 1.dat</em></p>
<p>Multiple data files can also be added together to create a larger ensemble. When this occurs the resulting data file with have the same name as the first input file with a "<em>c</em>" appended at the end. The width of the energy bining can also be increased and when this occurs a "<em>b</em>" is appended to the end.<p>
</html>
