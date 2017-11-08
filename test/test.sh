#!/usr/bin/env bash

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.1.1
# Date: 20171102
# Project: FISH optimal conditions
# Description: test the code installation.
# 
# ------------------------------------------------------------------------------



# ENV VAR ======================================================================

export LC_ALL=C

# PARAMS =======================================================================

basedir="`dirname ${BASH_SOURCE}`"
if [ "/" != ${basedir:0:1} ]; then basedir="$(pwd)/$basedir"; fi

	version="v1.2.0"

# RUN ==========================================================================

# Test requirements ------------------------------------------------------------
echo "Checking required software..."

# Check if metl.pl is installed
if [ ! -x "$(command -v melt.pl)" ]; then
  echo -e "fish-conditions requires OligoArrayAux melt.pl script to work."
  exit 1
fi

# Check if parallel is installed
if [ ! -x "$(command -v parallel)" ]; then
  echo -e "fish-conditions requires GNU parallel to parallelize."
  exit 1
fi

echo "...checked!"

# Test single scripts ----------------------------------------------------------
echo "Checking that scripts are properly installed..."

test1=$($basedir/../lib/score_temp.py --version)
if [ "score_temp.py $version" != "$test1" ]; then
	echo "score_temp.py not properly installed."
	exit 1
fi

test3=$($basedir/../src/find_fish_conditions.sh --version)
if [ "find_fish_conditions.sh $version" != "$test3" ]; then
	echo "find_fish_conditions.sh not properly installed."
	exit 1
fi

test3=$($basedir/../src/find_fish_conditions.single_probe.sh --version)
if [ "find_fish_conditions.single_probe.sh $version" != "$test3" ]; then
	echo "find_fish_conditions.single_probe.sh not properly installed."
	exit 1
fi

echo "...checked!"

# Run single-probe script ------------------------------------------------------
echo "Checking single-probe script..."

test4=$($basedir/../src/find_fish_conditions.single_probe.sh -y \
	--t1min 37 --t2min 37 --t1max 37 --t2max 37 \
	-i "$basedir/test.fa" -o "$basedir/out_test4/" --noplot)

test4std=$(cat $basedir/tests/test4.H1.picked.tsv)
test4=$(cat $basedir/out_test4/H1.picked.tsv)
if [ "$test4std" != "$test4" ]; then
	echo "Unexpected output from single-probe script. Broken installation."
	exit 1
fi

test4std=$(cat $basedir/tests/test4.H2.picked.tsv)
test4=$(cat $basedir/out_test4/H2.picked.tsv)
if [ "$test4std" != "$test4" ]; then
	echo "Unexpected output from single-probe script. Broken installation."
	exit 1
fi

rm -rf "$basedir/out_test4/"

echo "...checked!"

# Run multi-probe script -------------------------------------------------------
echo "Checking multi-probe script..."

test5=$($basedir/../src/find_fish_conditions.sh -y \
	--t1min 37 --t2min 37 --t1max 37 --t2max 37 \
	-r "s/^>.*(C[0-9]*:[^:]*):.*$/\\1/" --noplot \
	-i "$basedir/test.harmonize.fa" -o "$basedir/out_test5/")

test5std=$(cat $basedir/tests/test5.H1.picked.tsv)
test5=$(cat $basedir/out_test5/single_probes/C1:AFDN_cds/H1.picked.tsv)
if [ "$test5std" != "$test5" ]; then
	echo "Unexpected output from multi-probe script. Broken installation."
	exit 1
fi

test5std=$(cat $basedir/tests/test5.H2.picked.tsv)
test5=$(cat $basedir/out_test5/single_probes/C1:AFDN_cds/H2.picked.tsv)
if [ "$test5std" != "$test5" ]; then
	echo "Unexpected output from multi-probe script. Broken installation."
	exit 1
fi

rm -rf "$basedir/out_test5/"

echo "...checked!"

# Run probe harmonizing script -------------------------------------------------
echo "Checking multi-probe harmonization script..."

test6=$($basedir/../src/find_fish_conditions.sh -y --harmonize \
	--t1min 37 --t2min 37 --t1max 37 --t2max 37 \
	-r "s/^>.*(C[0-9]*:[^:]*):.*$/\\1/" --noplot \
	-i "$basedir/test.harmonize.fa" -o "$basedir/out_test6/")

test6std=$(cat $basedir/tests/test6.H1.picked.tsv)
test6=$(cat $basedir/out_test6//H1.picked.tsv)
if [ "$test6std" != "$test6" ]; then
	msg="Unexpected output from probe harmonization script. "
	msg=$msg"Broken installation."
	echo -e "$msg"
	exit 1
fi

test6std=$(cat $basedir/tests/test6.H2.picked.tsv)
test6=$(cat $basedir/out_test6//H2.picked.tsv)
if [ "$test6std" != "$test6" ]; then
	msg="Unexpected output from probe harmonization script. "
	msg=$msg"Broken installation."
	echo -e "$msg"
	exit 1
fi

rm -rf "$basedir/out_test6/"

echo "...checked!"

# END ==========================================================================

echo "All checked!"

################################################################################
