execute: prog
	./prog
prog: main.o spline_approximation.o
	gfortran $^ -o $@ -g
main.o: main.f90 spline_approximation.o
	gfortran $^ -c -g
spline_approximation.mod spline_approximation.o: spline_approximation.f90
	gfortran $^ -c -g
data: ./creator
	./creator
./creator: creator.f90
	gfortran creator.f90 -o creator
clean: 
	rm -f *.o *mod
graph: result.dat
	gnuplot plot_un.gnu
