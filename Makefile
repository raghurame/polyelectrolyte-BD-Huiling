all:
	gcc -o genRandomChain genRandomChain.c -lm -Wall
	./genRandomChain 60 1 300 300
