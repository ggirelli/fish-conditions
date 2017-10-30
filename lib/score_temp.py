#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 1.0.0
# Date: 20171027
# Project: FISH probe condition picking
# Description: Provide score for temperature condition.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse
import os
import re
import sys

from oligo_melting.lib.meltlib import *

# PARAMETERS ===================================================================


# Add script description
parser = argparse.ArgumentParser(description = '''
Calculates score of a temperature condition.
''')

# Add mandatory arguments
parser.add_argument('target_tsv', type = str, nargs = 1,
    help = 'Path to target hybridization characterization tsv file.')
parser.add_argument('second_tsv', type = str, nargs = 1,
    help = 'Path to secondary structure characterization tsv file.')

# Add arguments with default value
parser.add_argument('-t', type = float, nargs = 1,
    metavar = 'temp', help = """Condition temperature in Celsius degree.""",
    default = [37])
parser.add_argument('-d', type = str, nargs = 1,
    help = '''Duplex type. Possible values: DNA:DNA, RNA:RNA, DNA:RNA or
    RNA:DNA. The first nucleic acid type indicates the provided sequence.''',
    choices = ['DNA:DNA', 'RNA:RNA', 'RNA:DNA', 'DNA:RNA'],
    default = ['DNA:DNA'])
parser.add_argument('-o', metavar = "oligo_conc",
    type = float, nargs = 1,
    help = '''Oligonucleotide concentration [M].
    Default: 0.25e-6 M''',
    default = [0.25e-6])
parser.add_argument('-n', '--naconc', metavar = "na_conc",
    type = float, nargs = 1,
    help = '''Na+ concentration [M].
    Default: 50e-3 M''',
    default = [50e-3])
parser.add_argument('-m', '--mgconc', metavar = "mg_conc",
    type = float, nargs = 1,
    help = '''Mg2+ concentration [M]. Note: Mg2+ correction overwrites Na+
    correction. Default: 0 M''',
    default = [0])
parser.add_argument('-f', '--faconc', type = float, nargs = 1,
    metavar = 'fa_conc', help = '''Formamide concentration in %%(v,v).''',
    default = [35])
parser.add_argument('--fa-mode', type = str, nargs = 1,
    metavar = 'fa_mode', help = '''Mode of formamide correction. "mcconaughy"
    for classical -0.72%%FA correction from ref. 7 (default), "wright" for
    single reaction model correction from ref.8.''',
    choices = ['mcconaughy', 'wright'], default = ['mcconaughy'])
parser.add_argument('--fa-mvalue', type = str, nargs = 1, metavar = 'm',
    help = '''Specify the formamide m-value to be used with the wright
    correction model. Use either a single value "x" or two values "xL+y" where
    L is the probe length. Default: 0.1734''', default = ["0.1734"])
parser.add_argument('--out-single', type = str, nargs = 1, metavar = 'outfile',
    help = '''Path to output tsv file to store single oligo scores.''',
    default = [None])

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
tin_path = args.target_tsv[0]
sin_path = args.second_tsv[0]
temp = args.t[0]
dtype = args.d[0]

# Concentration
oligo_conc = args.o[0]
na_conc = args.naconc[0]
mg_conc = args.mgconc[0]

# Formamide
fa_conc = args.faconc[0]
fa_mode = args.fa_mode[0]
fa_mval_s = args.fa_mvalue[0]

# Output path
outfile = args.out_single[0]
doSingleOut = type(None) != type(outfile)

# Prepare buffer pointer and write header
if doSingleOut:
    fout = open(outfile, 'w')
    fout.write("name\tscore\n")

# Additional checks ------------------------------------------------------------

# Evaluate formamide m-value
parsed_mval = False
mregex = ["(([+-]?[0-9\.]*)L([+-][0-9\.]*))", "([+-]?[0-9\.]*)"]

# Search for xL+y m-value format
msearch = re.search(mregex[0], fa_mval_s)
if not type(None) is type(msearch):
    mgroups = msearch.groups()
    if fa_mval_s == mgroups[0]:
        mvalue = lambda x: float(mgroups[1]) * x + float(mgroups[2])
        parsed_mval = True

# Search for x m-value format
msearch = re.search(mregex[1], fa_mval_s)
if not parsed_mval and not type(None) is type(msearch):
    mgroup = msearch.groups()[0]
    if fa_mval_s == mgroup:
        mvalue = lambda x: float(mgroup)
        parsed_mval = True

# Trigger error otherwise
if not parsed_mval:
    sys.exit("!!!ERROR! Unexpected formamide m-value format Check help page.")

# Fix m-value label
if not fa_mode in ["wright"]:
    fa_mval_s = "-"

# FUNCTIONS ====================================================================

def read_tsv(path):
    # Reads a tsv creating a dictionary per each row.
    # The first row fields are used as keys.
    # 
    # Args:
    #   path (string): path to tsv file.
    # 
    # Returns:
    #   list: list of dictionaries.
    
    # Default delimiter for tsv files
    delim = "\t"

    # Output parsed row list
    out = []
    
    # Check if file exists
    if not os.path.exists(path):
        return([])

    # Open buffer to file
    with open(path) as fin:
        # Setup for different first row parsing
        first_row = []

        # Iterate content rows
        for row in fin:
            # Parse first row
            if 0 == len(first_row):
                first_row = row.strip().split(delim)
                continue
            
            # Parse other rows
            row = row.strip().split()

            # Check number of columns
            if not len(row) == len(first_row):
                print("ERROR! Variable column number found in %s" % path)
                return([])
            else:
                # Store row
                rowd = {}
                for i in range(len(first_row)):
                    rowd[first_row[i]] = row[i]
                out.append(rowd)

    # Output
    return(out)

def reparse_tsv_list(l, k, delim):
    # Reparse a list of dictionaries into a dictionary.
    # Uses the specified field as aggregation key.
    # 
    # Args:
    #   l (list): dictionary list.
    #   k (string): aggregation key.
    #   delim (string): to remove additional parts of the key.
    # 
    # Returns:
    #   dict: aggregated dictionaries.
    
    # Empty dictionary
    outd = {}

    # Iterate through list items
    for d in l:
        # All list items must have the aggregation key
        if not k in d.keys():
            msg = "ERROR! Not all the items contain the aggregation key '%s'"
            sys.exit(msg % k)
            return({})

        # Split on delimiter and keep prefix
        if delim in d[k]:
            d[k] = d[k].split(delim)[0]

        # Store
        if d[k] in outd.keys():
            outd[d[k]].append(d)
        else:
            outd[d[k]] = [d]

    # Output
    return(outd)

# RUN ==========================================================================

# Read tables
t = reparse_tsv_list(read_tsv(tin_path), 'oligo_name', ',')
s = reparse_tsv_list(read_tsv(sin_path), 'oligo_name', ',')

# Save lowest score
min_score = 1

# Iterate through the oligos
for oligo_name in t.keys():

    # Target -------------------------------------------------------------------
    tdh = float(t[oligo_name][0]['dH'])
    tds = float(t[oligo_name][0]['dS'])
    seq = t[oligo_name][0]['Seq'].upper()
    fgc = (seq.count('G') + seq.count('C')) / (len(seq) + 0.)

    # Reset melting temperature to standard conditions
    ttemp = temp + 273.15
    ttemp = duMelt_ion_adj(ttemp, 1, 0, fgc, na_conc)
    ttemp = duMelt_fa_adj(ttemp, tdh, tds, seq, oligo_conc,
        0, fa_mode, mvalue, dtype, fa_conc)

    # Calculate fraction
    tf = duMelt_curve(seq, oligo_conc, na_conc, mg_conc,
        fa_conc, fa_mode, mvalue, fgc, tdh, tds, ttemp, 0, 0.1, dtype)
    tf = tf[0][1]

    # Secondary structure ------------------------------------------------------
    sdh = float(s[oligo_name][0]['dH'])
    sds = float(s[oligo_name][0]['dS']) / 1000

    # Reset melting temperature to standard conditions
    stemp = ssMelt_fa_adj(temp + 273.15, 0, fa_conc)

    # Calculate fraction
    sf = ssMelt_curve(sdh, sds, stemp, fa_conc, 0, 0.1)
    sf = sf[0][1]

    # Calculate score ----------------------------------------------------------
    
    unfold = sf             # Unfolded fraction
    folded = 1 - unfold     # Folded fraction
    dissoc = tf             # Dissociated fraction
    hybrid = 1 - dissoc     # Hybridized fraction

    # Probability of good/bas oligos
    good = unfold * hybrid     # Both unfolded AND hybridized
    bad = (folded + dissoc) - folded * dissoc    # Folded OR dissociated

    # Score to be maximized
    score = good

    # Output if needed ---------------------------------------------------------
    if doSingleOut:
        fout.write("%s\t%f\n" % (oligo_name, score))

    # Save lowest --------------------------------------------------------------
    if score <= min_score:
        min_score = score

# Close output buffer
if doSingleOut:
    fout.close()

# Print result
print(min_score)

# END ==========================================================================

################################################################################
