#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 20171016
# Project: FISH probe condition picking
# Description: select optimal uniFISH 1st and 2nd hybridization conditions
#   for probes composed of oligonucleotides with the following structure:
#   forward-target-reverse-color.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# PARAMS =======================================================================


# Help string
helps="
usage: ./find_fish_conditions.sh [-h|--help][-v|--verbose]
    [--t1 temp][--t1step step][--fa1 conc][--fa1step step][--na1 conc]
    [--t2 temp][--t2step step][--fa2 conc][--fa2step step][--na2 conc]
    [-p conc][-u conc][-r pattern]
    -i fasta -o outdir

 Description:
  This script selects optimal uniFISH 1st and 2nd hybridization conditions for
  probes composed of oligonucleotides with the following structure:
   forward-target-reverse-color.
  Takes a fasta file in input, where the probe ID is specified in the header and
  identified through the provided -p pattern option.

 Mandatory arguments:
  -i fasta        Input fasta file.
  -o outdir       Output folder.

 Optional arguments:
  -h, --help      Show this help page.
  -v, --verbose   Verbose mode.
  --t1 temp       Default temperature for 1st hybridization. Default: 37 degC
  --t1step step   Step for 1st hyb. temp. exploration. Default: 0.5 degC
  --fa1 conc      Default formamide conc. for 1st hyb. Default: 25 %
  --fa1step step  Step for FA conc. exploration. Default: 5 %
  --na1 conc      Monovalent ion conc for 1st hyb. Default: 0.300 M
  --t2 temp       Default temperature for 2nd hybridization. Default: 37 degC
  --t2step step   Step for 2nd hyb. temp. exploration. Default: 0.5 degC
  --fa2 conc      Default formamide conc. for 2nd hyb. Default: 25%
  --fa2step step  Step for FA conc. exploration. Default: 5%
  --na2 conc      Monovalent ion conc for 2nd hyb. Default: 0.300 M
  -p conc         Probe concentration. Default: 1e-6 M
  -u conc         Universal (labeled) oligo concentration. Default: 1e-6 M
  -r pattern      Regular expression for probe name identification.
                  Default: ...
"

# Default values
t1=37
t1step=0.5
fa1=25
fa1step=5
na1=0.3
t2=37
t2step=0.5
fa2=25
fa2step=5
na2=0.3
probe_conc=0.000001
uni_conc=0.000001
pregexp=""
verbose=false

# Set option parsing strings
opt_name="find_fish_conditions.sh"
opt_short="hvi:o:p:u:r:"
opt_long="help,verbose,t1:,t1step:,fa1:,fa1step:,na1:,"
opt_long=$opt_long"t2:,t2step:,fa2:,fa2step:,na2:"

# Parse options
TEMP=`getopt -o $opt_short --long $opt_long -n $opt_name -- "$@"`
eval set -- "$TEMP"
while true ; do
    case "$1" in
        -h| --help) # Print help page
            echo -e "$helps"
        exit 0 ;;
        --t1) # 1st hybr. temp.
            if [ $2 -ge 0 ]; then t1="$2"; else
                msg="$helps\n!!!ERROR! --t1 cannot be lower than 0 degC."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        --t1step) # 1st hybr. temp. step
            if [ $2 -gt 0 ]; then t1step=$2; else
                msg="$helps\n!!!ERROR! --t1step must be higher than 0 degC."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        --fa1) # 1st hybr. formamide conc.
            if [ $2 -gt 0 ]; then fa1=$2; else
                msg="$helps\n!!!ERROR! --fa1 must be higher than 0 %."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        --fa1step) # 1st hybr. formamide conc. step
            if [ $2 -gt 0 ]; then fa1step=$2; else
                msg="$helps\n!!!ERROR! --fa1step must be higher than 0 %."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        --na1) # 1st hybr. monovalen ion conc.
            if [ $2 -gt 0 ]; then fa1=$2; else
                msg="$helps\n!!!ERROR! --na1 must be higher than 0 M."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        --t2) # 2nd hybr. temp.
            if [ $2 -ge 0 ]; then t2="$2"; else
                msg="$helps\n!!!ERROR! --t2 cannot be lower than 0 degC."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        --t2step) # 2nd hybr. temp. step
            if [ $2 -gt 0 ]; then t2step=$2; else
                msg="$helps\n!!!ERROR! --t2step must be higher than 0 degC."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        --fa2) # 2nd hybr. formamide conc.
            if [ $2 -gt 0 ]; then fa2=$2; else
                msg="$helps\n!!!ERROR! --fa2 must be higher than 0 %."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        --fa2step) # 2nd hybr. formamide conc. step
            if [ $2 -gt 0 ]; then fa2step=$2; else
                msg="$helps\n!!!ERROR! --fa2step must be higher than 0 %."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        --na2) # 2nd hybr. monovalen ion conc.
            if [ $2 -gt 0 ]; then fa2=$2; else
                msg="$helps\n!!!ERROR! --na2 must be higher than 0 M."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        -p) # Probe concentration
            if [ $2 -gt 0 ]; then probe_conc=$2; else
                msg="$helps\n!!!ERROR! -p must be higher than 0 M."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        -u) # Universal oligo concentration
            if [ $2 -gt 0 ]; then uni_conc=$2; else
                msg="$helps\n!!!ERROR! -u must be higher than 0 M."
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        -r) # Probe name regular expression
            pregexp="$2"
        shift 2 ;;
        -i) # Input fasta
            if [ -e "$2" ]; then fain_path="$2"; else
                msg="$helps\n!!!ERROR! Invalid -i option."
                  msg="$msg\n          File not found: $2\n"
                echo -e "$msg"; exit 1
            fi
        shift 2 ;;
        -o) # Output folder
            outdir="$2"
        shift 2 ;;
        -v| --verbose)
            # Verbose mode on
            verbose=true
            shift
        ;;
        --) shift ; break ;;
        *) echo "Internal error!" ; exit 1 ;;
    esac
done

# Check mandatory options
if [ -z "$fain_path" ]; then
    echo -e "$helps\n!!! ERROR! Missing mandatory -i option.\n"
    exit 1
fi
if [ -z "$outdir" ]; then
    echo -e "$helps\n!!! ERROR! Missing mandatory -o option.\n"
    exit 1
fi

# Print options
opt_string="
#------------ GENERAL ------------#

     Input fasta : $fain_path
   Output folder : $outdir
         Verbose : $verbose

#------- 1st HYBRIDIZATION -------#

         [probe] : $probe_conc M
     Temperature : $t1 degC
      Temp. step : $t1step degC
            [FA] : $fa1 %
       [FA] step : $fa1step %
           [Na+] : $na1 M

#------- 2nd HYBRIDIZATION -------#

           [uni] : $uni_conc M
     Temperature : $t2 degC
      Temp. step : $t2step degC
            [FA] : $fa2 %
       [FA] step : $fa2step %
           [Na+] : $na2 M

#---------------------------------#
"
echo -e "$opt_string"

# RUN ==========================================================================

# Read fasta and identify probes -----------------------------------------------
# Iterate from default conditions ----------------------------------------------
# 2nd structure and tm & FA
# 

# END ==========================================================================

################################################################################
