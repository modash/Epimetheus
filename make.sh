#!/bin/sh
gcc -Wall -Wextra -posix -Werror -std=c99 -D_POSIX_C_SOURCE=2001 -lm  bed_utils_light2.c -o bed_utils_light2
gcc -o genetersect genetersect.c
Rscript packageCheck.R
