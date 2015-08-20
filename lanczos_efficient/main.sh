#!/bin/bash
site_potentials="$(./site_pot.e)"

for n_up in {1..7};
do
	for n_dn in {1..7};
	do
		echo -e "$n_up \n $n_dn \n $site_potentials"| ./Zgroundstate.e
	done
done

MIN="$(./reader.e)"

echo -e "$MIN" | ./Zspectrum.e

echo "$MIN"