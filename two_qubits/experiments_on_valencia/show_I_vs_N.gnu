# python main.py &> out &
# grep "time" out > t; grep "shot" out > s; grep "Mutual" out > m; paste t s m > I_vs_shots; rm s; rm m; rm t

#set term postscript color enhanced font "Helvetica,24"
#set output "I_versus_N.eps"

#set xlabel "1/sqrt(Nshot)"
#set ylabel "I(A,B)"

f(x)=m*x+qqqSOBSTITUTETHIS
fit[:0.04] f(x) './I_vs_shotsSOBSTITUTETHIS' u (1.0/sqrt($1)):2:3 via m,qqqSOBSTITUTETHIS

#p[0:] './I_vs_shots' u (1.0/sqrt($2)):7:9 w errorlines title "QASM", 1.386292 title "STATEVECTOR", f(x) title "a+b/sqrt(N)", q+0.0142 notitle, q-0.0142 notitle
