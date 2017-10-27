#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ------------------------------------------------------------------------------
# 
# Author: ...
# Email: ...
# Version: X.X.X
# Date: YYYYMMDD
# Project: ...
# Description: ...
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

import argparse

# PARAMETERS ===================================================================


# Add script description
parser = argparse.ArgumentParser(
	description = 'Description...'
)

# Add mandatory arguments
parser.add_argument('test', type = str, nargs = 1,
	help = 'Example.')

# Add arguments with default value
parser.add_argument('-e', type = int, nargs = 1,
	metavar = 'e', help = """Example default.""", default = [0])

# Add flags
parser.add_argument('-u',
	action = 'store_const', dest = 'u',
	const = True, default = False,
	help = 'Example flag.')

# Parse arguments
args = parser.parse_args()

# Assign to in-script variables
test = args.test[0]
example = args.e[0]
flag = args.u

# FUNCTIONS ====================================================================

# RUN ==========================================================================

# END ==========================================================================

################################################################################
