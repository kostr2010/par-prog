CFLAGS = -Wall -pthread -lm

all : main.out clean

main.out : main.o
	gcc -o $@ $^ $(CFLAGS)

clean:
	rm *.o