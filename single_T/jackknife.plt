set terminal postscript eps color enhanced "Times" 20 size 4, 2

set output "result/jakknife_E.eps"

set xlabel "blocksize"
set ylabel "E_{err}"
set grid
plot "result/jackknife_result.txt" using 1:3 notitle pointtype 7 linecolor rgb 'red'

set output "result/jakknife_M.eps"

set xlabel "blocksize"
set ylabel "M_{err}"
set grid
plot "result/jackknife_result.txt" using 1:5 notitle pointtype 7 linecolor rgb 'red'

set output "result/jakknife_M2.eps"

set xlabel "blocksize"
set ylabel "M2_{err}"
set grid
plot "result/jackknife_result.txt" using 1:7 notitle pointtype 7 linecolor rgb 'red'

set output "result/jakknife_Cv.eps"

set xlabel "blocksize"
set ylabel "C_{v,err}"
set grid
plot "result/jackknife_result.txt" using 1:9 notitle pointtype 7 linecolor rgb 'red'

set output "result/jakknife_chi.eps"

set xlabel "blocksize"
set ylabel "{/Symbol c}_{v,err}"
set grid
plot "result/jackknife_result.txt" using 1:11 notitle pointtype 7 linecolor rgb 'red'

