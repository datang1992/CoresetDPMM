make:	*.o
*.o:	*.cpp
	g++ -Wall -I/opt/local/include -c *.cpp -fopenmp
	g++ -L/opt/local/lib *.o -lgsl -lgslcblas -lm -fopenmp -o main
