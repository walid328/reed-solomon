CC= gcc
FLAGS= -Wall -g
OBJFILES= main.o polynomial.o finite_field.o rs_code.o
TARGET= main

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(FLAGS) $^ -o $@

%.o : %.c
	$(CC) $(FLAGS) -c $^ -o $@

.PHONY:clean
clean:
	rm -r $(OBJFILES) $(TARGET)