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
usage: ./find_fish_conditions.single_probe.sh [-h|--help][-v|--verbose][-y]
    [--t1 temp][--t1step step][--fa1 conc][--fa1step step][--na1 conc]
    [--t2 temp][--t2step step][--fa2 conc][--fa2step step][--na2 conc]
    [--famode mode][-p conc][-u conc][-r pattern][-n pname][-t nthreads]
    -i fasta -o outdir

 Description:
  This script selects optimal uniFISH 1st and 2nd hybridization conditions for
  probes composed of oligonucleotides with the following structure:
   color-forward-target-reverse.
  Takes a fasta file in input, where each sequence is an oligo in a probe. This
  script works on a single probe and considers the provided fasta file as such.

 Notes:
  Monovalent ion (e.g., Na+) concentration must be > 0 M for 2ndary structure
  calculations using OligoArrayAux.

 Mandatory arguments:
  -i fasta        Input fasta file.
  -o outdir       Output folder.

 Optional arguments:
  -h, --help      Show this help page.
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
  -p conc         Probe concentration. Default: 1e-6 M
  -u conc         Universal (labeled) oligo concentration. Default: 1e-6 M
  -s struct       Comma separated color,forward,target,reverse length in nt.
                  Default: 20,20,30,20
  -n pname        Probe name. Default: 'probe'
  -t nthreads     Number of threads for parallelization. GNU parallel is
                  required for parallelization to occur.
"

# Default values
t1=37
t1step=0.1
t1min=32
t1max=42
fa1=25
na1=0.3
t2=37
t2step=0.1
t2min=32
t2max=42
fa2=25
na2=0.3
probe_conc=0.000001
uni_conc=0.000001
pregexp=""
verbose=false
dtype="DNA:RNA"
fa_mode="mcconaughy"
fa_mvalue=0.522
struct="20,20,30,20"
probe_name="probe"
nthreads=1
ask=true

# Set option parsing strings
opt_name="find_fish_conditions.sh"
opt_short="hvyi:o:u:p:s:n:t:"
opt_long="help,verbose,dtype:,famode:,mvalue:,t1:,t1step:,t1min:,t1max:,fa1:,"
opt_long=$opt_long"na1:,t2:,t2step:,t2min:,t2max:,fa2:,na2:"

# Parse options
TEMP=`getopt -o $opt_short --long $opt_long -n $opt_name -- "$@"`
eval set -- "$TEMP"
while true ; do
  case "$1" in
    -h| --help) # Print help page
      echo -e "$helps"
    exit 0 ;;
    -n) # Probe name
      probe_name="$2"
    shift 2 ;;
    --dtype) # Duplex type
      dtype="$2"
    shift 2 ;;
    --famode) # FA correction mode
      if [ "wright" == "$2" -o "mcconaughy" == "$2" ]; then fa_mode="$2"; else
        msg="$helps\n!!!ERROR! --famode possible values are"
        msg="$msg 'mcconaughy' and 'wright'."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --famvalue) # FA m-value
      fa_mvalue="$2"
    shift 2 ;;
    --t1) # 1st hybr. temp.
      if (( $(bc <<< "$2 >= 0") )); then t1="$2"; else
        msg="$helps\n!!!ERROR! --t1 cannot be lower than 0 degC."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --t1step) # 1st hybr. temp. step
      if (( $(bc <<< "$2 > 0") )); then t1step=$2; else
        msg="$helps\n!!!ERROR! --t1step must be higher than 0 degC."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --t1min) # 1st hybr. temp. lower boundary of exploration space
      if (( $(bc <<< "$2 > 0") )); then t1min=$2; else
        msg="$helps\n!!!ERROR! --t1min must be higher than 0 degC."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --t1max) # 1st hybr. temp. lower boundary of exploration space
      if (( $(bc <<< "$2 > 0") )); then t1max=$2; else
        msg="$helps\n!!!ERROR! --t1max must be higher than 0 degC."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --fa1) # 1st hybr. formamide conc.
      if (( $(bc <<< "$2 >= 0") )); then fa1=$2; else
        msg="$helps\n!!!ERROR! --fa1 must be higher than or equal to 0 %."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --na1) # 1st hybr. monovalen ion conc.
      if (( $(bc <<< "$2 > 0") )); then na1=$2; else
        msg="$helps\n!!!ERROR! --na1 must be higher than 0 M."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --t2) # 2nd hybr. temp.
      if (( $(bc <<< "$2 >= 0") )); then t2="$2"; else
        msg="$helps\n!!!ERROR! --t2 cannot be lower than 0 degC."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --t2step) # 2nd hybr. temp. step
      if (( $(bc <<< "$2 > 0") )); then t2step=$2; else
        msg="$helps\n!!!ERROR! --t2step must be higher than 0 degC."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --t2min) # 2nd hybr. temp. lower boundary of exploration space
      if (( $(bc <<< "$2 > 0") )); then t2min=$2; else
        msg="$helps\n!!!ERROR! --t2min must be higher than 0 degC."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --t2max) # 2nd hybr. temp. lower boundary of exploration space
      if (( $(bc <<< "$2 > 0") )); then t2max=$2; else
        msg="$helps\n!!!ERROR! --t2max must be higher than 0 degC."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --fa2) # 2nd hybr. formamide conc.
      if (( $(bc <<< "$2 >= 0") )); then fa2=$2; else
        msg="$helps\n!!!ERROR! --fa2 must be higher than or equal to 0 %."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --na2) # 2nd hybr. monovalen ion conc.
      if (( $(bc <<< "$2 > 0") )); then na2=$2; else
        msg="$helps\n!!!ERROR! --na2 must be higher than 0 M."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    -p) # Probe concentration
      if (( $(bc <<< "$2 > 0") )); then probe_conc=$2; else
        msg="$helps\n!!!ERROR! -p must be higher than 0 M."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    -u) # Universal oligo concentration
      if (( $(bc <<< "$2 > 0") )); then uni_conc=$2; else
        msg="$helps\n!!!ERROR! -u must be higher than 0 M."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    -s) # Probe structure
      struct="$2"
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
    -t) # Number of thread sfor parallelization
      if (( $(bc <<< "0 < $2") )); then nthreads=$2; else
        msg="$helps\n!!!ERROR! -t must be higher than 0."
        echo -e "$msg"; exit 1
      fi
    shift 2;;
    -v | --verbose)
      # Verbose mode on
      verbose=true
    shift ;;
    -y) # Ask for confirmation
    ask=false
    shift ;;
    --) shift ; break ;;
    *) echo "Internal error!" ; exit 1 ;;
  esac
done

# Check that metl.pl is installed
if [ ! -x "$(command -v melt.pl)" ]; then
  echo -e "fish-conditions requires OligoArrayAux melt.pl script to work."
  exit 1
fi
if [ 1 -lt $nthreads -a ! -x "$(command -v parallel)" ]; then
  echo -e "fish-conditions requires GNU parallel to parallelize."
  exit 1
else
  maxthreads=$(parallel --number-of-cores)
  if [ $nthreads -gt $maxthreads ]; then
    echo -e "Cannot use more than $maxthreads threads."
    exit 1
  fi
fi

# Check mandatory options
if [ -z "$fain_path" ]; then
  echo -e "$helps\n!!! ERROR! Missing mandatory -i option.\n"
  exit 1
fi
if [ -z "$outdir" ]; then
  echo -e "$helps\n!!! ERROR! Missing mandatory -o option.\n"
  exit 1
fi

# Additional manipulations
IFS=',' read -r -a astruct <<< "$struct"
if [ ${#astruct[@]} -ne 4 ]; then
  msg="""$helps\n !!! ERROR!"
  msg="$msg This script expect a structure of the form:"
  msg="$msg color-forward-target-reverse."
  msg="$msg The provided structure does not match the expected one: $struct"
  echo -e "$msg"
  exit 1
fi
exp_size=0
frag_start=(0)
for len in ${astruct[@]}; do
  if [ $len -le 0 -o -z "$len" ]; then
    msg="$helps\n!!! ERROR! Every structure part must be longer than 0: $struct"
    echo -e "$msg\n"
    exit 1
  fi
  exp_size=$(bc <<< "$exp_size + $len")
  frag_start+=($exp_size)
done

# Make paths absolute
if [ "/" != ${fain_path:0:1} ]; then fain_path=$(pwd)/$fain_path; fi
if [ "/" != ${outdir:0:1} ]; then outdir=$(pwd)/$outdir/; fi

if $ask; then
  # Print options
  opt_string="
  #------------ GENERAL ------------#

        Probe name : $probe_name
       Input fasta : $fain_path
     Output folder : $outdir
           Verbose : $verbose
     FA correction : $fa_mode
         Structure : $struct
     Expected size : $exp_size
           Threads : $nthreads

  #------- 1st HYBRIDIZATION -------#

           [probe] : $probe_conc M
       Temperature : $t1 degC
        Temp. step : $t1step degC
       Temp. range : $t1min - $t1max degC
              [FA] : $fa1 %
         [FA] step : $fa1step %
             [Na+] : $na1 M

  #------- 2nd HYBRIDIZATION -------#

             [uni] : $uni_conc M
       Temperature : $t2 degC
        Temp. step : $t2step degC
       Temp. range : $t2min - $t2max degC
              [FA] : $fa2 %
         [FA] step : $fa2step %
             [Na+] : $na2 M

  #---------------------------------#
  "

  # Ask confirmation
  settings_confirm=`echo -e "$opt_string" | sed 's/^/ /'`
  settings_confirm="
   ##############################################
   #                                            #
   #  PLEASE, DOUBLE CHECK THE SETTINGS BELOW   #
   # (press 'q' to continue when you are ready) #
   #                                            #
   ##############################################

  $settings_confirm"

  echo -e "$settings_confirm" | less

  msg="$msg\nRun the analysis?\nYes (y), Abort (a), Show again (s)"
  clear; echo -e $msg; read -e ans

  end=0
  while [[ 0 -eq $end ]]; do
    if [[ -z $ans ]]; then
      echo -e $msg
      read -e ans
    elif [[ 'a' == $ans ]]; then
      end=1
      echo "Aborted."
      exit 1
    elif [[ 'y' == $ans ]]; then
      echo -e "\n"
      end=1
    elif [[ 's' == $ans ]]; then
      echo -e "$settings_confirm" | less
      clear; echo -e $msg; read -e ans
    else
      echo -e $msg
      read -e ans
    fi
  done
fi

# RUN ==========================================================================

# Read fasta and identify probes -----------------------------------------------
# Iterate from default conditions ----------------------------------------------
# 2nd structure and tm & FA
# 

# END ==========================================================================

################################################################################
