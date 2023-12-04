all:
	gcc -o genRandomChain genRandomChain.c -lm -Wall
	./genRandomChain 110 77 0.9 0.9 100 300
