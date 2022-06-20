#! /usr/bin/gnuplot

set term qt font "arial,12" enhanced #persist
set encoding utf8
set tmargin 1
set tics out
set grid
set nokey

set samples 10000    # x-axis
set isosamples 10000 # y-axis

set style line 1 lt rgb "#ff0000" lw 1.5 #Scarlet, thick

plot "./wfscratch.dat" using 1:2 smooth csplines linestyle 1

pause mouse close
