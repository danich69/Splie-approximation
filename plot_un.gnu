#set terminal postfile
#set output "plot.ps"
set title "Interpolation"
set key left box
set xlabel "x"
set ylabel "f(x) approximated"
plot "result.dat" w l, "data.dat" using 1:2 w p
pause -1 "This is your plot"
