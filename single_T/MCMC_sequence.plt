set terminal png #size 4, 2

set output "result/MCMC_sequence_E.png"

set xlabel "MC step({/Symbol \264}10^4)"
set ylabel "E"
set grid
plot "result/MCMC_sequence.txt" using ($1/10000):2 linecolor rgb 'blue'  notitle

set output "result/MCMC_sequence_M.png"

set xlabel "MC step"
set ylabel "M"
set grid
plot "result/MCMC_sequence.txt" using ($1/10000):3 linecolor rgb 'blue' notitle

set output "result/MCMC_sequence_M2.png"

set xlabel "MC step"
set ylabel "M^2"
set grid
plot "result/MCMC_sequence.txt" using ($1/10000):4 linecolor rgb 'blue' notitle

