#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.0
# Date: 20171016
# Project: FISH probe condition picking
# Description: select optimal uniFISH 1st and 2nd hybridization conditions
#   for probes composed of oligonucleotides with the following structure:
#   color-forward-target-reverse.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

moddir="`dirname ${BASH_SOURCE}`/../lib/"
if [ "/" != ${moddir:0:1} ]; then moddir="$(pwd)/$moddir"; fi
srcdir="`dirname ${BASH_SOURCE}`/"
if [ "/" != ${srcdir:0:1} ]; then srcdir="$(pwd)/$srcdir"; fi
curdir="$(pwd)/"

# DEPENDENCIES =================================================================

source $moddir/find_fish_conditions.lib.sh

# PARAMS =======================================================================

# Help string
helps="
usage: ./find_fish_conditions.single_probe.sh [-h|--help][-v|--verbose][-y]
    [--dtype dtype][--famode mode][--mvalue m][--version][--noplot]
    [--t1 temp][--t1step step][--t1min tmin][--t1max tmax]
    [--fa1 fa][--na1 na][--mg1 mg]
    [--t2 temp][--t2step step][--t2min tmin][--t2max tmax]
    [--fa2 fa][--na2 na][--mg2 mg]
    [-p conc][-u conc][-r pattern][-s struct][-n pname][-t nthreads]
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
  --noplot        Do not produce any plot.
  -v, --verbose   Verbose mode.
  --version       Print version and stop.
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
  --mg1 conc      Divalent ion conc for 1st hyb. Default: 0 M
  --t2 temp       Default temperature for 2nd hybridization. Default: 37 degC
  --t2min tmin    Lower boundary for temperature exploration (2nd hybrid.).
                  Default: 32
  --t2max tmax    Upper boundary for temperature exploration (2nd hybrid.).
                  Default: 42
  --t2step step   Step for 2nd hyb. temp. exploration. Default: 0.1 degC
  --fa2 conc      Default formamide conc. for 2nd hyb. Default: 25%
  --na2 conc      Monovalent ion conc for 2nd hyb. Default: 0.300 M
  --mg2 conc      Divalent ion conc for 1st hyb. Default: 0 M
  -n pname        Probe name. Default: 'probe'
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
"

# Default values
t1=37
t1step=0.1
t1min=32
t1max=42
fa1=25
na1=0.3
mg1=0
t2=37
t2step=0.1
t2min=32
t2max=42
fa2=25
na2=0.3
mg2=0
probe_conc=0.000001
uni_conc=0.000001
probe_regexp="s/^[>\ ]*([^:]*):.*/\\1/"
verbose=false
dtype="DNA:RNA"
fa_mode="mcconaughy"
fa_mvalue=0.522
struct="20,20,30,20"
probe_name="probe"
nthreads=1
ask=true
doplot=true

# Set option parsing strings
opt_name="find_fish_conditions.sh"
opt_short="hvyi:o:u:p:s:n:t:r:"
opt_long="help,verbose,dtype:,famode:,mvalue:,t1:,t1step:,t1min:,t1max:,fa1:,"
opt_long=$opt_long"na1:,mg1:,t2:,t2step:,t2min:,t2max:,fa2:,na2:,mg2:,"
opt_long=$opt_long"version,noplot"

# Parse options
TEMP=`getopt -o $opt_short --longoptions $opt_long -n $opt_name -- "$@"`
eval set -- "$TEMP"
while true ; do
  case "$1" in
    -h| --help) # Print help page
      echo -e "$helps"
    exit 0 ;;
    --noplot) # No plot mode
      doplot=false
    shift ;;
    --version) # Print version and stop
      echo "find_fish_conditions.single_probe.sh v1.1.0"
    exit 0 ;;
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
    --mvalue) # FA m-value
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
    --mg1) # 1st hybr. divalent ion conc.
      if (( $(bc <<< "$2 >= 0") )); then mg1=$2; else
        msg="$helps\n!!!ERROR! --mg1 must be higher than or equal to 0 M."
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
    --mg2) # 2nd hybr. divalent ion conc.
      if (( $(bc <<< "$2 >= 0") )); then mg2=$2; else
        msg="$helps\n!!!ERROR! --mg2 must be higher than or equal to 0 M."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    -i) # Input fasta
      if [ -e "$2" ]; then fain_path="$2"; else
        msg="$helps\n!!!ERROR! Invalid -i option."
        msg="$msg\n          File not found: $2\n"
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    -n) # Probe name
      probe_name="$2"
    shift 2 ;;
    -o) # Output folder
      outdir="$2"
    shift 2 ;;
    -p) # Probe concentration
      if (( $(bc <<< "$2 > 0") )); then probe_conc=$2; else
        msg="$helps\n!!!ERROR! -p must be higher than 0 M."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    -r) # Probe ID pattern
      probe_regexp="$2"
    shift 2 ;;
    -s) # Probe structure
      struct="$2"
    shift 2 ;;
    -t) # Number of thread sfor parallelization
      if (( $(bc <<< "0 < $2") )); then nthreads=$2; else
        msg="$helps\n!!!ERROR! -t must be higher than 0."
        echo -e "$msg"; exit 1
      fi
    shift 2;;
    -u) # Universal oligo concentration
      if (( $(bc <<< "$2 > 0") )); then uni_conc=$2; else
        msg="$helps\n!!!ERROR! -u must be higher than 0 M."
        echo -e "$msg"; exit 1
      fi
    shift 2 ;;
    -v | --verbose)
      # Verbose mode on
      verbose=true
    shift ;;
    -y) # Ask for confirmation
    ask=false
    shift ;;
    --) shift ; break ;;
    *) echo "Internal error! Unrecognized option '$1'." ; exit 1 ;;
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
              Plot : $doplot

  #------- 1st HYBRIDIZATION -------#

           [probe] : $probe_conc M
       Temperature : $t1 degC
        Temp. step : $t1step degC
       Temp. range : $t1min - $t1max degC
              [FA] : $fa1 %
             [Na+] : $na1 M
            [Mg2+] : $mg1 M

  #------- 2nd HYBRIDIZATION -------#

             [uni] : $uni_conc M
       Temperature : $t2 degC
        Temp. step : $t2step degC
       Temp. range : $t2min - $t2max degC
              [FA] : $fa2 %
             [Na+] : $na2 M
            [Mg2+] : $mg2 M

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

# Create output directory if missing
mkdir -p "$outdir"

# Copy input fasta to output directory
if [ "$fain_path" != "$outdir/input.fa" ]; then
  cp $fain_path $outdir/input.fa
fi

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

echo -e "Found $n_oligo $olen nt oligos."

# Extract only target sequences from the fasta file ----------------------------

targs=$(echo -e "$fain_seq" | \
  awk -v start=${frag_start[2]} -v len=${astruct[2]} \
  '{ print substr($0, start, len) }')
paste <(echo -e "$fain_id") <(echo -e "$targs") | tr '\t' '\n' \
  > "$outdir/targets.fa"

# 1st HYBRIDIZATION ============================================================

echo -e "
# 1st HYBRIDIZATION
# ================="

# Iterate at different temperature ---------------------------------------------

# Prepare output
fname="$outdir/H1.temp.score.tsv"
if [ -e $fname ]; then
  echo -e "t\tscore\tnorm" > $fname
fi

if [ 1 == $nthreads ]; then # SINGLE THREAD #

  # Run default condition
  echo -e "Checking default temperature"
  run_single_condition1 $outdir $probe_name $fa1 $na1 $mg1 $probe_conc \
    $fa_mvalue $fa_mode $dtype $t1 \
    "$moddir" "$srcdir" "$outdir/input.fa" 0 $doplot
  echo -e "$ct\t$cscore\t$nscore" >> $fname

  # Default best values
  best_score=$cscore
  best_t=$t1
  best_cond="H1_"$probe_name"_FA"$fa1"p_Na"$na1"M_t"$best_t"degC"

  # Explore lower temperatures
  echo -e "\nExploring lower temperatures"
  ct=$(bc <<< "$t1 - $t1step")
  while
    run_single_condition1 $outdir $probe_name $fa1 $na1 $mg1 $probe_conc \
      $fa_mvalue $fa_mode $dtype $ct \
      "$moddir" "$srcdir" "$outdir/input.fa" 0 $doplot
    (( $(bc <<< "$ct >= $t1min" ) ));
    #(( $(bc <<< "$cscore >= $best_score") ));
  do
    if (( $(bc <<< "$cscore >= $best_score") )); then
      best_score=$cscore
      best_cond=$cond_string
      best_t=$ct
    fi
    echo -e "$ct\t$cscore\t$nscore" >> $fname
    ct=$(bc <<< "$ct - $t1step")
  done
  echo -e " · Reached lower boundary, stop.\n"

  # Explore higher temperatures
  echo -e "Exploring higher temperatures"
  ct=$(bc <<< "$t1 + $t1step")
  while
    run_single_condition1 $outdir $probe_name $fa1 $na1 $mg1 $probe_conc \
      $fa_mvalue $fa_mode $dtype $ct \
      "$moddir" "$srcdir" "$outdir/input.fa" 0 $doplot
    (( $(bc <<< "$ct <= $t1max" ) ));
    #(( $(bc <<< "$cscore >= $best_score") ));
  do
    if (( $(bc <<< "$cscore >= $best_score") )); then
      best_score=$cscore
      best_cond=$cond_string
      best_t=$ct
    fi
    echo -e "$ct\t$cscore\t$nscore" >> $fname
    ct=$(bc <<< "$ct + $t1step")
  done
  echo -e " · Reached upper boundary, stop.\n"

else # MULTI-THREAD #

  # Explore lower temperatures
  echo -e "\nExploring temperature space"

  # Run in parallel
  export -f run_single_condition1
  pout=$(parallel -kj $nthreads run_single_condition1 ::: $outdir ::: \
    "'$probe_name'" ::: $fa1 ::: $na1 ::: $mg1 ::: $probe_conc ::: $fa_mvalue ::: \
    $fa_mode ::: $dtype ::: $(seq $t1min $t1step $t1max) ::: "$moddir" ::: \
    "$srcdir" ::: "$outdir/input.fa" ::: 1 ::: $doplot)
  
  # Reformat with temperature
  pout=$(paste <(seq $t1min $t1step $t1max) <(echo -e "$pout"))
  echo -e "$pout" >> $fname

  # Select best copndition
  cscore=$(echo -e "$pout" | datamash max 2)
  best_score=$cscore
  best_t=$(echo -e "$pout" | grep $best_score | cut -f1)
  best_cond="H1_"$probe_name"_FA"$fa1"p_Na"$na1"M_t"$best_t"degC"
  
fi

# Select -----------------------------------------------------------------------

echo -e " · Best condition: $best_cond"
echo -e " >>> Score: $best_score\n"

# Focus on best condition
cond_dir="$outdir/$best_cond"
ct=$best_t

# Identify FA concentration for optimal temperature ----------------------------

optimal_fa1=$(bc <<< "scale = 2; $fa1 + (($ct - $t1) / 0.72)")
echo -e "For optimal 1st hybridization at $t1 degC, use $optimal_fa1% FA.\n"

# Plot per-oligo coupled melting curves ----------------------------------------

if $doplot; then
  echo -e " · Plotting single-oligo description for optimal condition..."
  $moddir/oligo_melting/scripts/plot_melt_curves_coupled.R -n $probe_name \
    -t $best_t "$cond_dir/targets.melt_curve.$ct.FA"$fa1"p.tsv" \
    "$cond_dir/second.melt_curve.$ct.FA"$fa1"p.tsv" \
    "$cond_dir/oligo.melt_curve.$ct.FA"$fa1"p.pdf" & pid=$!
  wait $pid
  mv "$cond_dir/oligo.melt_curve.$ct.FA"$fa1"p.pdf" \
    "$outdir/H1.oligo.melt_curve.optimal.pdf"
fi

# Save conditions --------------------------------------------------------------

echo -e "probe\tscore\th%FA\tTh\to%FA\tTo\tNa\tprobe_conc" \
  > "$outdir/H1.picked.tsv"
line="$probe_name\t$best_score\t$optimal_fa1\t$t1\t$fa1"
line=$line"\t$best_t\t$na1\t$probe_conc"
echo -e "$line" >> "$outdir/H1.picked.tsv"

# Move to subfolder
if [ -d "$outdir/H1" ]; then rm -r $outdir/H1; fi
mkdir -p $outdir/H1
mv $outdir/H1_* $outdir/H1/

# 2nd HYBRIDIZATION ============================================================

echo -e "
# 2nd HYBRIDIZATION
# ================="

if [ -d "$outdir/H2/" ]; then rm -r "$outdir/H2/"; fi
mkdir -p "$outdir/H2/"

# Extract color sequences
col_seq=($(echo -e "$fain_seq" | awk \
  -v s=${frag_start[0]} -v e=${frag_start[1]} \
  '{ print substr($1, s, e) }' | sort | uniq))
col_id=()
for s in ${col_seq[@]}; do col_id+=($(cat "$outdir/input.fa" | paste - - | \
  grep "$s" | cut -f1 | sed -E "$probe_regexp" | sort | uniq)); done
paste <(echo ${col_id[@]} | tr ' ' '\n') <(echo ${col_seq[@]} | tr ' ' '\n') | \
  sed 's/^/>/' | tr '\t' '\n' > "$outdir/color.fa"

# Extract color-forward sequences
colfor_seq=$(echo -e "$fain_seq" | awk \
  -v s=${frag_start[0]} -v e=${frag_start[2]} \
  '{ print substr($1, s, e) }' | sort | uniq)
for s in ${colfor_seq[@]}; do colfor_id+=($(cat "$outdir/input.fa" | \
  paste - - | grep "$s" | cut -f1 | sed -E "$probe_regexp" | sort | uniq));
done
paste <(echo ${colfor_id[@]} | tr ' ' '\n') <(echo ${colfor_seq[@]} | \
  tr ' ' '\n') | sed 's/^/>/' | tr '\t' '\n' > "$outdir/color.forward.fa"

# Calculate hybridization of target portions -----------------------------------
$moddir/oligo_melting/melt_duplex.py "$outdir/targets.fa" -FC \
    -o $probe_conc -n $na2 -m $mg2 -f $fa2 --fa-mode $fa_mode -t $dtype \
    --fa-mvalue $fa_mvalue --t-curve 30 0.5 \
    --out-curve "$outdir/H2/targets.melt_curve.$t2.FA"$fa2"p.tsv" \
    > "$outdir/H2/targets.melt.$t2.FA"$fa2"p.tsv" & pid=$!
wait $pid

if [ ! -e "$outdir/H2/targets.melt.$t2.FA"$fa2"p.tsv" ]; then exit 1; fi

if $doplot; then
  # Plot melting curves
  $moddir/oligo_melting/scripts/plot_melt_curves.R \
      -n "Harmonized targets" \
      "$outdir/H2/targets.melt_curve.$t2.FA"$fa2"p.tsv" \
      "$outdir/H2/targets.melt_curve.$t2.FA"$fa2"p.pdf" & pid=$!
  wait $pid
fi

# Iterate at different temperature =============================================

# Prepare output
fname="$outdir/H2.temp.score.tsv"
if [ -e $fname ]; then
  echo -e "t\tscore\tnorm" > $fname
fi

if [ 1 == $nthreads ]; then # SINGLE THREAD #

  # Run default condition
  echo -e "Checking default temperature"
  run_single_condition2 $outdir $probe_name $fa2 $na2 $mg2 $uni_conc \
    $fa_mvalue $fa_mode "DNA:DNA" $t2 "$moddir" "$srcdir" \
    "$outdir/color.forward.fa" 0 $t2 $doplot
  echo -e "$ct\t$cscore\t$nscore" >> $fname

  # Default best values
  best_score=$cscore
  best_t=$t2
  best_cond="H2_"$probe_name"_FA"$fa2"p_Na"$na2"M_t"$best_t"degC"

  # Explore lower temperatures
  echo -e "\nExploring lower temperatures"
  ct=$(bc <<< "$t2 - $t2step")
  while
    run_single_condition2 $outdir $probe_name $fa2 $na2 $mg2 $uni_conc \
      $fa_mvalue $fa_mode "DNA:DNA" $ct \
      "$moddir" "$srcdir" "$outdir/color.fa" 0 $t2 $doplot
    (( $(bc <<< "$ct >= $t2min" ) ));
    #(( $(bc <<< "$cscore >= $best_score") ));
  do
    if (( $(bc <<< "$cscore >= $best_score") )); then
      best_score=$cscore
      best_cond=$cond_string
      best_t=$ct
    fi
    echo -e "$ct\t$cscore\t$nscore" >> $fname
    ct=$(bc <<< "$ct - $t2step")
  done
  echo -e " · Reached lower boundary, stop.\n"

  # Explore higher temperatures
  echo -e "Exploring higher temperatures"
  ct=$(bc <<< "$t2 + $t2step")
  while
    run_single_condition2 $outdir $probe_name $fa2 $na2 $mg2 $uni_conc \
      $fa_mvalue $fa_mode "DNA:DNA" $ct \
      "$moddir" "$srcdir" "$outdir/color.fa" 0 $t2 $doplot
    (( $(bc <<< "$ct <= $t2max" ) ));
    #(( $(bc <<< "$cscore >= $best_score") ));
  do
    if (( $(bc <<< "$cscore >= $best_score") )); then
      best_score=$cscore
      best_cond=$cond_string
      best_t=$ct
    fi
    echo -e "$ct\t$cscore\t$nscore" >> $fname
    ct=$(bc <<< "$ct + $t2step")
  done
  echo -e " · Reached upper boundary, stop.\n"

else # MULTI-THREAD #

  # Explore lower temperatures
  echo -e "\nExploring temperature space"

  # Run in parallel
  export -f run_single_condition2
  pout=$(parallel -kj $nthreads run_single_condition2 ::: $outdir ::: \
    "'$probe_name'" ::: $fa2 ::: $na2 ::: $mg2 ::: $probe_conc ::: $fa_mvalue \
    ::: $fa_mode ::: "DNA:DNA" ::: $(seq $t2min $t2step $t2max) ::: "$moddir" \
    ::: "$srcdir" ::: "$outdir/color.fa" ::: 1 ::: $t2 ::: $doplot)
  
  # Reformat with temperature
  pout=$(paste <(seq $t2min $t2step $t2max) <(echo -e "$pout"))
  echo -e "$pout" >> $fname

  # Select best copndition
  cscore=$(echo -e "$pout" | datamash max 2)
  best_score=$cscore
  best_t=$(echo -e "$pout" | grep $best_score | cut -f1)
  best_cond="H2_"$probe_name"_FA"$fa2"p_Na"$na2"M_t"$best_t"degC"

fi

# Select -----------------------------------------------------------------------

echo -e " · Best condition: $best_cond"
echo -e " >>> Score: $best_score\n"

# Focus on best condition
cond_dir="$outdir/$best_cond"
ct=$best_t

# Identify FA concentration for optimal temperature ----------------------------

optimal_fa2=$(bc <<< "scale = 2; $fa2 + (($ct - $t2) / 0.72)")
echo -e "For optimal 2nd hybridization at $t2 degC, use $optimal_fa2% FA.\n"

# Plot per-oligo coupled melting curves ----------------------------------------

if $doplot; then
  echo -e " · Plotting single-oligo description for optimal condition..."
  $moddir/oligo_melting/scripts/plot_melt_curves_coupled.R -n $probe_name \
    -t $best_t "$cond_dir/color.melt_curve.$ct.FA"$fa2"p.tsv" \
    "$cond_dir/second.melt_curve.$ct.FA"$fa2"p.tsv" \
    --addit-tsv "$outdir/H2/targets.melt_curve.$t2.FA"$fa2"p.tsv" \
    "$cond_dir/h2.melt_curve.$ct.FA"$fa2"p.pdf" & pid=$!
  wait $pid
  mv "$cond_dir/h2.melt_curve.$ct.FA"$fa2"p.pdf" \
    "$outdir/H2.h2.melt_curve.optimal.pdf"
fi

# Save conditions --------------------------------------------------------------

echo -e "probe\tscore\th%FA\tTh\to%FA\tTo\tNa\tprobe_conc" \
  > "$outdir/H2.picked.tsv"
line="$probe_name\t$best_score\t$optimal_fa2\t$t2\t$fa2"
line=$line"\t$best_t\t$na2\t$probe_conc"
echo -e "$line" >> "$outdir/H2.picked.tsv"


# Move to subfolder
mv $outdir/H2_* $outdir/H2/

# END ==========================================================================

echo -e "\nDONE"

################################################################################
