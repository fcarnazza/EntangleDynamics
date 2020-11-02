set term postscript eps color enhanced font "Helvetica, 12"
set output "valencia_mutual_information10.eps"

set xlabel "time"
set ylabel "mutual information"
set multiplot layout 2 rowsfirst
# --- GRAPH a
@TMARGIN; @LMARGIN
@NOXTICS; @YTICS
set label 1 'a' @POS
plot f(x) with lines ls 1
# --- GRAPH b
@TMARGIN; @RMARGIN
@NOXTICS; @NOYTICS
set label 1 'b' @POS
plot g(x) with lines ls 1
unset multiplot
p './statevector_mutual_information10.txt' every  ::1::195  u 1:2 w lines title "statevector", './valencia10.txt' u 1:2:4 w errorlines title "extrapolated from ibmq valencia" 
