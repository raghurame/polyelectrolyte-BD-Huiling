all:
	gcc -o genRandomChain genRandomChain.c -lm -Wall
	./genRandomChain 60 0.9 1 50
