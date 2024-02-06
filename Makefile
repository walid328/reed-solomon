CC= gcc
FLAGS= -Wall -g
OBJ= build/polynomial.o build/finite_field.o build/field_settings.o build/array.o build/rs_code.o
TEST_TARGET= build/test_polynomial build/test_rs_code build/main 

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
	./build/test_polynomial "new_free"
	./build/test_polynomial "add"
	./build/test_polynomial "mul"
	./build/test_polynomial "copy"
	./build/test_polynomial "rev"
	./build/test_polynomial "deriv"
	./build/test_polynomial "mul_scalar"
	./build/test_polynomial "euc_div"
	./build/test_polynomial "xgcd"
	./build/test_polynomial "interpol"
	./build/test_polynomial "dft"
	./build/test_polynomial "fft"
	./build/test_polynomial "fast_mul"
	./build/test_polynomial "fast_euc_div"
	./build/test_polynomial "fast_xgcd"
	./build/test_rs_code "encode"
	./build/test_rs_code "fast_encode"
	./build/test_rs_code "decode"
	./build/test_rs_code "encode_decode"
	./build/test_rs_code "encode_decode_2"
	./build/test_rs_code "fast_encode_decode"