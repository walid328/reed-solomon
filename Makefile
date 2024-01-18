CC= gcc
FLAGS= -Wall -g
OBJFILES= polynomial.o finite_field.o rs_code.o field_spec.o
TARGET= main

all: $(TARGET)

$(TARGET): $(TARGET).o $(OBJFILES)
	$(CC) $(FLAGS) $^ -o $@

test_polynomial: test_polynomial.o $(OBJFILES)
	$(CC) $(FLAGS) $^ -o $@

%.o : %.c
	$(CC) $(FLAGS) -c $^ -o $@

.PHONY:clean test

clean:
	rm -r *.o $(TARGET)

test:
	./test_polynomial "new_free"
	./test_polynomial "add"
	./test_polynomial "mul"