fish-conditions
===

This project contains code to select optimal uniFISH 1st and 2nd hybridization conditions for probes composed of oligonucleotides with the following structure: color sequence (C), forward sequence (F), target sequence (T), and reverse sequence (R), from 5' to 3'. Below, I will refer to the longest flap comprising R and C as L (i.e., L=R+C) and to the whole oligo as O (i.e., O=C+F+T+R).

<img src="https://github.com/ggirelli/fish-conditions/blob/master/images/fish.png" alt="octocat" />

### Installation

1. Install OligoArrayAux if not already available. To check if the software is already available run `command -v melt.pl` or `melt.pl --version`. The first will give no output if the software is not installed, and the second will show the version of the currently installed software, if any. OligoArrayAux is available [here](http://unafold.rna.albany.edu/OligoArrayAux.php).
2. Install GNU parallel if not already available. To check if the software is already available run `command -v parallel` or `parallel --version`. The first will give no output if the software is not installed, and the second will show the version of the currently installed software, if any. GNU parallel is available [here](https://www.gnu.org/software/parallel/).
3. Download `fish-conditions` repository and submodules.

```
git clone https://github.com/ggirelli/fish-conditions
cd fish-conditions
git submodule init
git submodule update
```

4. Test the installation with `./test/test.sh` from within the main code folder (e.g., `fish-conditions`).

#### For Mac users

Homebrew package manager is needed to install required software. To install homebrew, check out [this page](https://brew.sh/).

If not yet available, please install `gnu-getopt` as follows:

```
brew install gnu-getopt
brew link --force gnu-getopt
```

Also, gnu-parallel (step 2 above) can be easily install with:

```
brew install parallel
brew link --force parallel
```

Please, notice that `fish-conditions v1.1.0` has been tested only with `parallel v20161222`.

## Usage

Before starting the analysis, the setting are showed and, after pressing "**q**", the user has the option to either review the settings again (**s**), abort the analysis (**a**) or confirm the settings and start the run (**y**).

* **Consult help page**. Run `find_fish_conditions.single_probe.sh -h` or `find_fish_conditions.sh -h` for single- or multi-probe manual, respectively.

* **Single probe evaluation**. Use `find_fish_conditions.single_probe.sh -i input.fa -o output_dir/` to run the default parameters over a single probe fasta file.

* **Multi-probe evaluation**. Use `find_fish_conditions.sh -i input.fa -o output_dir/` to run the default parameters over a single probe fasta file.

* **Multi-probe harmonization**. Use `find_fish_conditions.sh -i input.fa -o output_dir/ --harmonize` to run the default parameters over a single probe fasta file.

### First run

Plot generation is implemented in R and requires the R packages `argparse` and `ggplot2`. The script automatically checks for the packages, and tries to automatically install them, providing help in manually installing them if the automatic installation fails.

## License

```
MIT License
Copyright (c) 2017 Gabriele Girelli
```
