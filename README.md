fish-conditions
===

This project contains code to select optimal uniFISH 1st and 2nd hybridization conditions for probes composed of oligonucleotides with the following structure:

- color sequence (C)
- forward sequence (F)
- target sequence (T)
- reverse sequence (R)

I will refer to the longest flap comprising R and C as L (i.e., L=R+C) ans as the whole oligo as O (i.e., O=C+F+T+R).

### Considerations on the two hybridization steps

Default hybridization conditions are, in general: 300 mM \[Na+] (corresponding to 2xSSC), oligo concentration, 37 &deg;C and 25% formamide (FA).

For the first hybridization (H1), the script will identify the optimal hybridization temperature and formamide concentration that will maximize the ratio between the number of hybridized oligos without secondary structure and the total number of oligonucleotides.

For the second hybridization (H2), the script will identify the optimal hybridization temperature and formamide concentration that will maximize the ratio between the number of hybridized flaps without secondary structure and the total number of flaps, keeping into account the original melting temperature of the target to avoid oligo detachment.

### Parameters:

- RNA/DNA FISH
- Probe concentration, default at 1 &micro;M
- Universal oligo concentration, default at 1 &micro;M
- Sodium concentration, default at 300 mM
- For H1 and H2, separately:
    + Temperature, default at 37 &deg;C
    + Temperature step, default at 0.5 &deg;C (increase and decrease)
    + FA concentration, default at 25%
    + FA concentration step, default at 5% (only increase)

### Scripts:

- `find_fish_conditions.sh` main bash script. Runs OligoArrayAux to analyze secondary structures and uses `oligo-melting-temperature` for first and second hybridization melting temperature calculations.