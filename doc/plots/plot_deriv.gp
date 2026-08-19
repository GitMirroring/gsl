#!/usr/bin/env gnuplot

set logscale x
set format x "10^{%L}"
unset key
set grid back
set xlabel "step size h"

set term pngcairo enh col size 1000,480

## adaptive differentiation diagnostic f(x) = x^(3/2), x=1

set out "../images/deriv_1.png"
file = '../examples/diff_diag_x1.txt'

set multiplot layout 1,2

set ylabel "derivative estimate"

set title "Derivative estimate vs h, x=1"
plot file index 1 us 1:2 w lp lt -1 pt 7

set ylabel "error estimate"
set logscale y
set format y "10^{%L}"

set title "Error estimate vs h, x=1"
plot file index 1 us 1:3 w lp lt -1 pt 7

unset logscale y
set format y "%g"

unset multiplot

## adaptive differentiation diagnostic f(x) = x^(3/2), x=0

set out "../images/deriv_2.png"
file = '../examples/diff_diag_x0.txt'

set multiplot layout 1,2

set ylabel "derivative estimate"

set title "Derivative estimate vs h, x=0"
plot file index 1 us 1:2 w lp lt -1 pt 7

set ylabel "error estimate"
set logscale y
set format y "10^{%L}"

set title "Error estimate vs h, x=0"
plot file index 1 us 1:3 w lp lt -1 pt 7

unset multiplot
