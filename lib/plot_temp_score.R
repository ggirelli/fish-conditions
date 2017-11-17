#!/usr/bin/Rscript

# ------------------------------------------------------------------------------
# 
# Author: Gabriele Girelli
# Email: gigi.ga90@gmail.com
# Version: 2.0.0
# Date: 20171118
# Description: plot behaviour of temperature condition score.
# 
# ------------------------------------------------------------------------------



# DEPENDENCIES =================================================================

suppressMessages({
	if (!require('argparser')) {
		install.packages("argparser", repos = "http://cran.uk.r-project.org")
		if (!require('argparser')) {
			stop(paste0('The R package argparser was not found. ',
				'Tried to automatically install the package but failed.',
				'Install it manually with install.packages(argparser)'))
		}
	}
	if (!require('ggplot2')) {
		install.packages("ggplot2", repos = "http://cran.uk.r-project.org")
		if (!require('ggplot2')) {
			stop(paste0('The R package ggplot2 was not found. ',
				'Tried to automatically install the package but failed.',
				'Install it manually with install.packages(ggplot2)'))
		}
	}
})

# INPUT ========================================================================

# Create arguent parser
parser = arg_parser('Plot condition score behaviour.',
	name = 'plot_temp_score.R')

# Define mandatory arguments
parser = add_argument(parser, arg = 'temp_score_tsv',
	help = 'Path to temp.score.tsv file.')
parser = add_argument(parser, arg = 'output_pdf',
	help = 'Path to pdf output file.')

# Parse arguments
p = parse_args(parser)

# Attach argument values to variables
attach(p['' != names(p)])

# RUN ==========================================================================

t = read.delim(temp_score_tsv, as.is = T, header = T)

pdf(output_pdf, width = 10, height = 10)

p2 = ggplot(t, aes(x = t, y = norm, col = hybr)) + geom_line()
vline = data.frame(t = c(t$t[which.max(t$score[t$hybr == 'H1'])],
		t$t[which.max(t$score[t$hybr == 'H2'])]),
	hybr = c('H1', 'H2'), stringsAsFactors = F)
p2 = p2 + geom_vline(data = vline, aes(xintercept = t, color = hybr),
	linetype = 2)
p2 = p2 + xlab("Temperature [C]") + ylab("Normalized score")
print(p2)

p1 = ggplot(t, aes(x = t, y = score, col = hybr)) + geom_line()
vline = data.frame(t = c(t$t[which.max(t$score[t$hybr == 'H1'])],
		t$t[which.max(t$score[t$hybr == 'H2'])]),
	hybr = c('H1', 'H2'), stringsAsFactors = F)
p1 = p1 + geom_vline(data = vline, aes(xintercept = t, color = hybr),
	linetype = 2)
p1 = p1 + xlab("Temperature [C]") + ylab("Score")
print(p1)

graphics.off()



# END --------------------------------------------------------------------------

################################################################################
