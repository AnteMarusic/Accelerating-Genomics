CC = gcc
CFLAGS = -Wall -g
OUTPUT = PairHMM

all: $(OUTPUT)

OBJECTS = main.o algo.o list.o


#algo.o: algo.h algo.c #when this changes do
		#$(CC) $(CFLAGS) -c algo.c #this

#list.o: list.h list.c
		#$(CC) $(CFLAGS) -c list.c

$(OBJECTS): %.o: %.c
		$(CC) $(CFLAGS) -c $< -o $@

$(OUTPUT): algo.o list.o main.c
		$(CC) $(CFLAGS) -o $(OUTPUT) main.c list.o algo.o

run:
		./$(OUTPUT)

clean:
		rm *.o $(OUTPUT)
