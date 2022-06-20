#! /usr/bin/gnuplot

set term qt font "arial,12" enhanced
set encoding utf8
set tmargin 1
set tics out
set grid

set samples 10000    # x-axis
set isosamples 10000 # y-axis

set colorsequence podo

plot "./k01wfData.dat" using 1:2 smooth csplines, \
"./k02wfData.dat" using 1:2 smooth csplines, \
"./k03wfData.dat" using 1:2 smooth csplines, \
"./k04wfData.dat" using 1:2 smooth csplines, \
"./k05wfData.dat" using 1:2 smooth csplines, \
"./k06wfData.dat" using 1:2 smooth csplines, \
"./k07wfData.dat" using 1:2 smooth csplines, \
"./k08wfData.dat" using 1:2 smooth csplines, \
"./k09wfData.dat" using 1:2 smooth csplines, \
"./k10wfData.dat" using 1:2 smooth csplines, \
"./k11wfData.dat" using 1:2 smooth csplines, \
"./k12wfData.dat" using 1:2 smooth csplines, \
"./k13wfData.dat" using 1:2 smooth csplines


pause mouse close
