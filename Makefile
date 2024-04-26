CFLAGS = -O3

a.out:	main.o algorithm.cpp functions.cpp solution.cpp reduce_sum.cpp
		g++ $^ $(CFLAGS) -o $@

main.o: main.cpp inc.h
		g++ -c $(CFLAGS) main.cpp

clean:
		rm -f *.out *.o
		