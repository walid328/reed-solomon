CC= gcc
FLAGS= -Wall -g
OBJ= build/polynomial.o build/finite_field.o build/field_settings.o build/array.o build/rs_code.o
TEST_TARGET= build/test_polynomial build/test_rs_code build/test_compare build/test_performances build/main

all: $(TEST_TARGET)

$(TEST_TARGET): % : %.o $(OBJ)
	$(CC) $(FLAGS) $^ -o $@

build/%.o : %.c
	mkdir -p $(dir $@)
	$(CC) $(FLAGS) -o $@ $< -c

.PHONY:clean test

clean:
	rm -r build

test:
	./build/test_polynomial 193 "polynomial"
	./build/test_rs_code 193 "rs_code"
	./build/test_compare 193 "compare"

performances:
	./build/test_performances