fish-conditions
===

This project contains code to select optimal uniFISH 1st and 2nd hybridization conditions for probes composed of oligonucleotides with the following structure: color sequence (C), forward sequence (F), target sequence (T), and reverse sequence (R), from 5' to 3'. Below, I will refer to the longest flap comprising R and C as L (i.e., L=R+C) and to the whole oligo as O (i.e., O=C+F+T+R).

### FISH overview

<img src="https://github.com/ggirelli/fish-conditions/blob/master/images/fish.png" alt="octocat" />

Our FISH protocol comprises a first hybridization round in which and unlabeled oligo is hybridized to its target, and a second hybridization round where a conjugate oligo is hybridized to the color flap.

#### Considerations on the two hybridization steps

Default hybridization conditions are, in general: 300 mM \[Na+] (corresponding to 2xSSC), 1 uM probe concentration, 37 &deg;C and 25% formamide (FA).

For the **first hybridization** (H1), the script will identify the optimal hybridization temperature, i.e., the temperature that will maximize the fraction of **hybridized** (T, as opposed to **dissociated**) and **unfolded** (as opposed to **folded**) oligos (O).

For the **second hybridization** (H2), the script will identify the optimal hybridization temperature, i.e., the temperature that will maximize the fraction of **hybridized** conjugated oligos (C) over **unfolded** color-forward (CF) flaps, while avoiding whole oligo (O) dissociation.

In both cases, a formamide concentration is suggested to shift the optimal temperature to the default one, i.e., the user will be free to decide whether to change the incubation temperature or the formamide concentration in the buffers.

### Requirements

- GNU parallel, tested with version 20161222.
- OligoArrayAux, tested with version 3.8.

### Scripts:

- `find_fish_conditions.sh`: main bash script. Runs OligoArrayAux to analyze secondary structures and uses `oligo-melting` for first and second hybridization melting temperature calculations. Can be used to iterate over a multi-probe fasta, or to harmonize the conditions of multiple probes to be used in the same experiment.
- `find_fish_conditions.single_probe.sh`: bash script for single-probe evaluation.

### Parameters:

Both scripts can explore a range of temperature, with specified increments, around the default hybridization temperature.

```
 Mandatory arguments:
  -i fasta        Input fasta file.
  -o outdir       Output folder.

 Optional arguments:
  -h, --help      Show this help page.
  -a, --harmonize Find optimal condition for the provided probes to be used in
                  the same experiment. Thus, only one condition output.
                  Specific for multi-probe script only.
  -v, --verbose   Verbose mode.
  -y              Do not ask for settings confirmation.
  --dtype type    Duplex type: DNA:DNA, RNA:RNA, DNA:RNA, RNA:DNA.
                  Default: DNA:RNA
  --famode mode   Formamide correction mode: 'mcconaughy' (classic) or 'wright'.
                  Default: mcconaughy.
  --mvalue m      Formamide m-value, either single x value or xL+y format.
                  Only used with wrighte famode. Default: 0.522.
  --t1 temp       Default temperature for 1st hybridization. Default: 37 degC
  --t1step step   Step for 1st hyb. temp. exploration. Default: 0.1 degC
  --t1min tmin    Lower boundary for temperature exploration (1st hybrid.).
                  Default: 32
  --t1max tmax    Upper boundary for temperature exploration (1st hybrid.).
                  Default: 42
  --fa1 conc      Default formamide conc. for 1st hyb. Default: 25 %
  --na1 conc      Monovalent ion conc for 1st hyb. Default: 0.300 M
  --t2 temp       Default temperature for 2nd hybridization. Default: 37 degC
  --t2min tmin    Lower boundary for temperature exploration (2nd hybrid.).
                  Default: 32
  --t2max tmax    Upper boundary for temperature exploration (2nd hybrid.).
                  Default: 42
  --t2step step   Step for 2nd hyb. temp. exploration. Default: 0.1 degC
  --fa2 conc      Default formamide conc. for 2nd hyb. Default: 25%
  --na2 conc      Monovalent ion conc for 2nd hyb. Default: 0.300 M
  -n pname        Probe name. Default: 'probe'
                  Specific for single-probe script only.
  -p conc         Probe concentration. Default: 1e-6 M
  -r pattern      Regular expression for probe name identification.
                  Note: the regular expression output should be a substring of
                        the input fasta header.
                  Default: s/^[>\ ]*([^:]*):.*/\\1/
  -s struct       Comma separated color,forward,target,reverse length in nt.
                  Default: 20,20,30,20
  -t nthreads     Number of threads for parallelization. GNU parallel is
                  required for parallelization to occur.
  -u conc         Universal (labeled) oligo concentration. Default: 1e-6 M
```

### Background

Tabulated thermodynamic estimates for duplex hybridization, without mismatches, used in the [oligo-melting](http://github.com/ggirelli/oligo-melting/) package are available for: DNA:DNA (Allawi&SantaLucia, Biochemistry(36), 1997), DNA:RNA (or RNA:DNA; Sugimoto N, et al., Biochemistry(34), 1995), and RNA:RNA (Freier SM, et al., PNAS(83), 1986).

Mono- (Na+) and di-valent (Mg2+) ion concentration effects
are well characterized on the melting temperature. Formulas used in [oligo-melting](http://github.com/ggirelli/oligo-melting/) come from the following literature:

* Na+: Owczarzy R, et al., Biochemistry(43), 2004
* Mg2+: Owczarzy R, et al., Biochemistry(47), 2008

Formamide (denaturant) effect is well characterized on the free energy
(based on *m-value* estimates; Wright ES, et al., Appl Environ Microb(80), 2014). Classical model instead predicted -0.72 · [FA%] degC (McConaughy BL, et al., Biochemistry(8), 1969). A nice paper on *m-value*s: Myers JK, et al., Protein Science(4), 1995.

OligoArrayAux (part of UNAFold) is used to study the secondary structure (Markham NR, et al., Bioinformatics, 2008). The software already includes mono- (Na+) and di-valent (Mg2+) ion concentration effects  at the level of free energy.

As formamide effects on secondary structure are not well described in literature, only the classical correction model (McConaugy) is allowed in the current implementation.

### Score interpretation

The current implementation uses as score the **sum** of the fraction of *good* oligos (unfolded & hybridized) in a probe. As such, dividing the score by the number of oligos provides the average fraction of *good* oligos at given hybridization conditions.

For example, a score of 18.56 for a prove of 24 oligonucleotides means that on average only 77.33% of each oligo will be properly hybridized.

```
18.56 / 24 = 0.7733
```

Previous implementation focused on the **minimum** score of a probe, instead of the sum. The old approach resulted in lower averages, triggering the implementation of a new approach that could take all oligos into account at the same time.

### Usages

#### Single probe evaluation

...

#### Multi-probe evaluation

...

#### Multi-probe harmonization

...

### Setup

```
...
```