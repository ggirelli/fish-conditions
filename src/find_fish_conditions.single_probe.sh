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
#   color-forward-target-reverse.
#   
# TODO:
#  - Check if melt.pl is available
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# PARAMS =======================================================================

moddir="`dirname $(pwd)/${BASH_SOURCE}`/../lib/"
srcdir="`dirname $(pwd)/${BASH_SOURCE}`/"
curdir="$(pwd)/"

# Help string
helps="
usage: ./find_fish_conditions.single_probe.sh [-h|--help][-v|--verbose]
    [--t1 temp][--t1step step][--fa1 conc][--fa1step step][--na1 conc]
    [--t2 temp][--t2step step][--fa2 conc][--fa2step step][--na2 conc]
    [--famode mode][-p conc][-u conc][-r pattern][-n pname]
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
  --dtype type    Duplex type: DNA:DNA, RNA:RNA, DNA:RNA, RNA:DNA.
                  Default: DNA:RNA
  --famode mode   Formamide correction mode: 'mcconaughy' (classic) or 'wright'.
                  Default: mcconaughy.
  --mvalue m      Formamide m-value, either single x value or xL+y format.
                  Only used with wrighte famode. Default: 0.522.
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
  -s struct       Comma separated color,forward,target,reverse length in nt.
                  Default: 20,20,30,20
  -n pname        Probe name. Default: 'probe'
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
dtype="DNA:RNA"
fa_mode="mcconaughy"
fa_mvalue=0.522
struct="20,20,30,20"
probe_name="probe"

# Set option parsing strings
opt_name="find_fish_conditions.sh"
opt_short="hvi:o:u:p:s:n:"
opt_long="help,verbose,dtype:,famode:,mvalue:,t1:,t1step:,fa1:,fa1step:,na1:,"
opt_long=$opt_long"t2:,t2step:,fa2:,fa2step:,na2:"

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
    --fa1) # 1st hybr. formamide conc.
      if (( $(bc <<< "$2 >= 0") )); then fa1=$2; else
        msg="$helps\n!!!ERROR! --fa1 must be higher than or equal to 0 %."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --fa1step) # 1st hybr. formamide conc. step
      if (( $(bc <<< "$2 > 0") )); then fa1step=$2; else
        msg="$helps\n!!!ERROR! --fa1step must be higher than 0 %."
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
    --fa2) # 2nd hybr. formamide conc.
      if (( $(bc <<< "$2 >= 0") )); then fa2=$2; else
        msg="$helps\n!!!ERROR! --fa2 must be higher than or equal to 0 %."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    --fa2step) # 2nd hybr. formamide conc. step
      if (( $(bc <<< "$2 > 0") )); then fa2step=$2; else
        msg="$helps\n!!!ERROR! --fa2step must be higher than 0 %."
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
    -v| --verbose)
      # Verbose mode on
      verbose=true
    shift ;;
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

# Create output directory if missing
mkdir -p "$outdir"

# Copy input fasta to output directory
cp $fain_path $outdir/input.fa

# Read fasta and print a summary -----------------------------------------------
fain=$(cat "$fain_path")
fain_id=$(echo -e "$fain" | grep ">")
fain_seq=$(echo -e "$fain" | grep -v ">")

# Number of oligos
n_oligo=$(echo -e "$fain_seq" | wc -l)

# Oligo length(s)
olen=($(echo -e "$fain_seq" | awk '{ print length($0) }' | sort | uniq))
if [ 1 -ne ${#olen[@]} ]; then
  msg="!!!ERROR! Inconsistent oligo length detected: ${olen[@]}"
  echo -e "$helps\n$msg"
  exit 1
fi
if [ $olen -ne $exp_size ]; then
  msg="!!!ERROR! Oligo length does not match the provided structure."
  msg="$msg Expecting $exp_size nt oligos, found $olen nt oligos instead."
  echo -e "$helps\n$msg"
  exit 1
fi

echo -e " · Found $n_oligo $olen nt oligos."

# Extract only target sequences from the fasta file ----------------------------
targs=$(echo -e "$fain_seq" | \
  awk -v start=${frag_start[2]} -v len=${astruct[2]} \
  '{ print substr($0, start, len) }')
paste <(echo -e "$fain_id") <(echo -e "$targs") | tr '\t' '\n' \
  > "$outdir/targets.fa"









# Iterate at different temperature =============================================

ct=$t1











# Log current conditions
echo -e " · Working on: [FA] $fa1%; [Na] $na1 M; t $ct degC" 

# Generate condition folder
cond_string=$probe_name"_FA"$fa1"p_Na"$na1"M_t"$ct"degC"
cond_dir="$outdir/$cond_string"
mkdir -p "$cond_dir"

# Calculate hybridization of target portions -----------------------------------
$moddir/oligo_melting/melt_duplex.py "$outdir/targets.fa" -FC \
  -o $probe_conc -n $na1 -f $fa1 --fa-mode $fa_mode -t $dtype \
  --fa-mvalue $fa_mvalue --t-curve 30 0.5 \
  --out-curve "$cond_dir/targets.melt_curve.$ct.FA"$fa1"p.tsv" \
  > "$cond_dir/targets.melt.$ct.FA"$fa1"p.tsv"
if [ ! -e "$cond_dir/targets.melt.$ct.FA"$fa1"p.tsv" ]; then exit 1; fi

# Plot melting curves
$moddir/oligo_melting/scripts/plot_melt_curves.R -n "$cond_string : Targets" \
  "$cond_dir/targets.melt_curve.$ct.FA"$fa1"p.tsv" \
  "$cond_dir/targets.melt_curve.$ct.FA"$fa1"p.pdf"

# 2nd structure and tm & FA ----------------------------------------------------
cd $cond_dir
cp $outdir/input.fa $cond_dir/input.fa

# Use OligoArrayAux for 2nd structure calculation
melt.pl -n DNA -t $ct -N $na1 -C $probe_conc input.fa >> "oligo_melt.tsv.tmp"

# Re-format data for easy manipulation
melt_id="oligo_name\n$(cat "oligo_melt.tsv.tmp" | grep "Calculating")"
melt_data=$(cat "oligo_melt.tsv.tmp" | grep -v "Calculating")
paste <(echo -e "$melt_id") <(echo -e "$melt_data") | \
  sed 's/^Calculating for//' | sed 's/ //g' | \
  paste - <(echo -e "Seq\n$fain_seq") > "second.melt.$ct.tsv"
rm "oligo_melt.tsv.tmp"

# FA correction
$moddir/oligo_melting/melt_second.py -f $fa1 --t-curve 60 0.5 \
  --out-curve "$cond_dir/second.melt_curve.$ct.FA"$fa1"p.tsv" -C \
  $cond_dir/"second.melt.$ct.tsv" > $cond_dir/"second.melt.$ct.FA"$fa1"p.tsv"

# Plot secondary structure melting curves
$moddir/oligo_melting/scripts/plot_melt_curves.R \
  -n "$cond_string : Secondary structure" \
  "$cond_dir/second.melt_curve.$ct.FA"$fa1"p.tsv" \
  "$cond_dir/second.melt_curve.$ct.FA"$fa1"p.pdf"

# Score function ---------------------------------------------------------------
cscore=$($moddir/score_temp.py -d "$dtype" -t $ct -o $probe_conc -n $na1 \
  -f $fa1  --fa-mode "$fa_mode" --fa-mvalue "$fa_mvalue" \
  --out-single "$cond_dir/oligo.scores.$ct.FA"$fa1"p.tsv" \
  "$cond_dir/targets.melt.$ct.FA"$fa1"p.tsv" \
  "$cond_dir/second.melt.$ct.FA"$fa1"p.tsv")

rm input.fa*
cd $curdir

# Identify FA concentration for optimal temperature ============================

# Plot per-oligo coupled melting curves ========================================
# $moddir/oligo_melting/scripts/plot_melt_curves_coupled.R -n $probe_name -t $t1 \
#   "$cond_dir/targets.melt_curve.$ct.FA"$fa1"p.tsv" \
#   "$cond_dir/second.melt_curve.$ct.FA"$fa1"p.tsv" \
#   "$cond_dir/oligo.melt_curve.$ct.FA"$fa1"p.pdf"

# !!!Only for optimal condition!!!

# END ==========================================================================

################################################################################
