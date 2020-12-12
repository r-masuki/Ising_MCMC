set terminal postscript eps color enhanced "Times" 20 size 4, 2

set output "result/E.eps"

set xlabel "T"
set ylabel "E/L^2"
set grid
plot "result/T_dependence.txt" using 1:2:3 with yerrorbars notitle pointtype 7 linecolor rgb 'red'

set output "result/M.eps"

set xlabel "T"
set ylabel "m"
set grid
plot "result/T_dependence.txt" using 1:4:5 with yerrorbars notitle pointtype 7 linecolor rgb 'red'

set output "result/M2.eps"

set xlabel "T"
set ylabel "m^2"
set grid
plot "result/T_dependence.txt" using 1:6:7 with yerrorbars notitle pointtype 7 linecolor rgb 'red'

set output "result/Cv.eps"

set xlabel "T"
set ylabel "c_{v}"
set grid
plot "result/T_dependence.txt" using 1:8:9 with yerrorbars notitle pointtype 7 linecolor rgb 'red'

set output "result/chi.eps"

set xlabel "T"
set ylabel "{/Symbol c}_{v}"
set grid
plot "result/T_dependence.txt" using 1:10:11 with yerrorbars notitle pointtype 7 linecolor rgb 'red'

