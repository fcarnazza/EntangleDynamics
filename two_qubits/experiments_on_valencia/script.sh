#!/bin/bash
#script to extrapolate over n1,n2,...N shots.
t=0
until [ $t -gt 195 ]
do
  awk '/'"time $t shots"'/{ print $4,$5,$6; }' val10 > I_vs_shots$t
  sed 's/SOBSTITUTETHIS/'$t'/g' show_I_vs_N.gnu > show_I_vs_N$t.gnu
  gnuplot show_I_vs_N$t.gnu
  rm show_I_vs_N$t.gnu
  rm I_vs_shots$t
  ((t=t+1))
done
grep "+/-" fit.log > t; paste t > mutual1; rm t
sed 's/qqq/q /g' mutual1 > mutual2
awk '/'"q"'/{ print $2/100,$4,$5,$6; }' mutual2 > mutual3
gnuplot show.gnu
