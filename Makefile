all: conv, zero, pi

run_pi: pi
	mpirun -np 4 --use-hwthread-cpus ./pi

run_zero: zero
	mpirun -np 2 --use-hwthread-cpus ./zero

run_conv: conv
	mpirun -np 4 --use-hwthread-cpus ./conv_par 500 500| python3 visualize.py
conv:
	$(CC) -Wall -o conv_seq conv.c -lm
	mpicc -Wall -o conv_par conv_parallel.c -lm

zero:
	mpicc -Wall -o zero zero.c
pi:
	mpicc -Wall -o pi pi.c

clean:
	rm -f conv_par conv_seq pi zero 

