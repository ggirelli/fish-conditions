#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 20171030
# Description: functions for FISH hybridization condition picking
# 
# ------------------------------------------------------------------------------



# FUNCTIONS ====================================================================

function run_single_condition1() {

    # First hybridization analysis
    outdir=${1}
    probe_name=$(echo "${2}" | sed -r "s/'//g")
    fa1=${3}
    na1=${4}
    probe_conc=${5}
    fa_mvalue=${6}
    fa_mode=${7}
    dtype=${8}
    ct=${9}
    moddir=${10}
    srcdir=${11}
    fain_path=${12}
    fain_seq=$(cat "$fain_path" | grep -v ">")
    parallel=${13}
    doplot=${14}

    # Log current conditions
    if [ 0 -eq $parallel ]; then
        echo -e " · Working on: [FA] $fa1%; [Na] $na1 M; t $ct degC"
    fi

    # Generate condition folder
    cond_string="H1_"$probe_name"_FA"$fa1"p_Na"$na1"M_t"$ct"degC"
    cond_dir="$outdir/$cond_string"
    mkdir -p "$cond_dir"

    # Calculate hybridization of target portions -------------------------------
    $moddir/oligo_melting/melt_duplex.py "$outdir/targets.fa" -FC \
        -o $probe_conc -n $na1 -f $fa1 --fa-mode $fa_mode -t $dtype \
        --fa-mvalue $fa_mvalue --t-curve 30 0.5 \
        --out-curve "$cond_dir/targets.melt_curve.$ct.FA"$fa1"p.tsv" \
        > "$cond_dir/targets.melt.$ct.FA"$fa1"p.tsv" & pid=$!
    wait $pid

    if [ ! -e "$cond_dir/targets.melt.$ct.FA"$fa1"p.tsv" ]; then exit 1; fi

    if [ "true" == "$dplot" ]; then
        # Plot melting curves
        $moddir/oligo_melting/scripts/plot_melt_curves.R \
            -n "$cond_string : Targets" \
            "$cond_dir/targets.melt_curve.$ct.FA"$fa1"p.tsv" \
            "$cond_dir/targets.melt_curve.$ct.FA"$fa1"p.pdf" & pid=$!
        wait $pid
    fi

    # 2nd structure and tm & FA ------------------------------------------------
    cd $cond_dir
    cp $outdir/input.fa $cond_dir/input.fa

    # Use OligoArrayAux for 2nd structure calculation
    melt.pl -n DNA -t $ct -N $na1 -C $probe_conc input.fa \
        >> "oligo_melt.tsv.tmp" & pid=$!
    wait $pid

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
        $cond_dir/"second.melt.$ct.tsv" \
        > $cond_dir/"second.melt.$ct.FA"$fa1"p.tsv" & pid=$!
    wait $pid

    if [ "true" == "$dplot" ]; then
        # Plot secondary structure melting curves
        $moddir/oligo_melting/scripts/plot_melt_curves.R \
            -n "$cond_string : Secondary structure" \
            "$cond_dir/second.melt_curve.$ct.FA"$fa1"p.tsv" \
            "$cond_dir/second.melt_curve.$ct.FA"$fa1"p.pdf" & pid=$!
        wait $pid
    fi

    # Score function -----------------------------------------------------------
    cscore=$($moddir/score_temp.py -d "$dtype" -t $ct -o $probe_conc -n $na1 \
        -f $fa1  --fa-mode "$fa_mode" --fa-mvalue "$fa_mvalue" \
        --out-single "$cond_dir/oligo.scores.$ct.FA"$fa1"p.tsv" \
        "$cond_dir/targets.melt.$ct.FA"$fa1"p.tsv" \
        "$cond_dir/second.melt.$ct.FA"$fa1"p.tsv")

    if [ 0 -eq $parallel ]; then
        echo -e " >>> Score: $cscore"
    else
        echo $cscore
    fi

    rm input.fa*
}

function run_single_condition2() {

    # First hybridization analysis
    outdir=${1}
    probe_name=$(echo "${2}" | sed -r "s/'//g")
    fa2=${3}
    na2=${4}
    probe_conc=${5}
    fa_mvalue=${6}
    fa_mode=${7}
    dtype=${8}
    ct=${9}
    moddir=${10}
    srcdir=${11}
    fain_path=${12}
    fain_seq=$(cat "$fain_path" | grep -v ">")
    parallel=${13}
    t2=${14}
    doplot=${15}

    # Log current conditions
    if [ 0 -eq $parallel ]; then
        echo -e " · Working on: [FA] $fa2%; [Na] $na2 M; t $ct degC"
    fi

    # Generate condition folder
    cond_string="H2_"$probe_name"_FA"$fa2"p_Na"$na2"M_t"$ct"degC"
    cond_dir="$outdir/$cond_string"
    mkdir -p "$cond_dir"

    # Calculate hybridization of color portions -------------------------------
    $moddir/oligo_melting/melt_duplex.py "$outdir/color.fa" -FC \
        -o $probe_conc -n $na2 -f $fa2 --fa-mode $fa_mode -t $dtype \
        --fa-mvalue $fa_mvalue --t-curve 30 0.5 \
        --out-curve "$cond_dir/color.melt_curve.$ct.FA"$fa2"p.tsv" \
        > "$cond_dir/color.melt.$ct.FA"$fa2"p.tsv" & pid=$!
    wait $pid

    if [ ! -e "$cond_dir/color.melt.$ct.FA"$fa2"p.tsv" ]; then exit 1; fi

    if [ "true" == "$doplot" ]; then
        # Plot melting curves
        $moddir/oligo_melting/scripts/plot_melt_curves.R \
            -n "$cond_string : Color" \
            "$cond_dir/color.melt_curve.$ct.FA"$fa2"p.tsv" \
            "$cond_dir/color.melt_curve.$ct.FA"$fa2"p.pdf" & pid=$!
        wait $pid
    fi

    # 2nd structure and tm & FA ------------------------------------------------
    cd $cond_dir
    cp $outdir/color.forward.fa $cond_dir/color.forward.fa

    # Use OligoArrayAux for 2nd structure calculation
    melt.pl -n DNA -t $ct -N $na2 -C $probe_conc color.forward.fa \
        >> "colfor_melt.tsv.tmp" & pid=$!
    wait $pid

    # Re-format data for easy manipulation
    melt_id="oligo_name\n$(cat "colfor_melt.tsv.tmp" | grep "Calculating")"
    melt_data=$(cat "colfor_melt.tsv.tmp" | grep -v "Calculating")
    paste <(echo -e "$melt_id") <(echo -e "$melt_data") | \
        sed 's/^Calculating for//' | sed 's/ //g' | \
        paste - <(echo -e "Seq\n$fain_seq") > "second.melt.$ct.tsv"
    rm "colfor_melt.tsv.tmp"

    # FA correction
    $moddir/oligo_melting/melt_second.py -f $fa2 --t-curve 60 0.5 \
        --out-curve "$cond_dir/second.melt_curve.$ct.FA"$fa2"p.tsv" -C \
        $cond_dir/"second.melt.$ct.tsv" \
        > $cond_dir/"second.melt.$ct.FA"$fa2"p.tsv" & pid=$!
    wait $pid

    if [ "true" == "$doplot" ]; then
        # Plot secondary structure melting curves
        $moddir/oligo_melting/scripts/plot_melt_curves.R \
            -n "$cond_string : CF secondary structure" \
            "$cond_dir/second.melt_curve.$ct.FA"$fa2"p.tsv" \
            "$cond_dir/second.melt_curve.$ct.FA"$fa2"p.pdf" & pid=$!
        wait $pid
    fi

    # Score function -----------------------------------------------------------

    cscore=$($moddir/score_temp.py -d "$dtype" -t $ct -o $probe_conc -n $na2 \
        -f $fa2  --fa-mode "$fa_mode" --fa-mvalue "$fa_mvalue" \
        --out-single "$cond_dir/oligo.scores.$ct.FA"$fa2"p.tsv" \
        --addit "$outdir/H2/targets.melt.$t2.FA"$fa2"p.tsv" \
        "$cond_dir/color.melt.$ct.FA"$fa2"p.tsv" \
        "$cond_dir/second.melt.$ct.FA"$fa2"p.tsv")

    if [ 0 -eq $parallel ]; then
        echo -e " >>> Score: $cscore"
    else
        echo $cscore
    fi

    rm color.forward.fa*
}

################################################################################