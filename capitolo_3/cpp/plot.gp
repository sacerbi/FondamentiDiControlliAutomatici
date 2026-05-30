set title 'Analisi del Sistema'
set xlabel 'Tempo (t)'
set ylabel 'Ampiezza'
set grid
plot 'dati.txt' using 1:2 with lines lw 2 title "Impulso", \
     'dati.txt' using 1:3 with lines lw 2 title "Rampa"
