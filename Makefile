CC= gcc
FLAGS= -Wall -g

SRC_DIR= src
TEST_DIR= test
BUILD_DIR= build

SRC= $(notdir $(wildcard $(SRC_DIR)/*.c))
TEST= $(notdir $(wildcard $(TEST_DIR)/*.c))
MAIN= main.c

SRC_OBJ= $(addprefix $(BUILD_DIR)/, $(SRC:.c=.o))
TEST_OBJ= $(addprefix $(BUILD_DIR)/, $(TEST:.c=.o))
MAIN_OBJ= $(addprefix $(BUILD_DIR)/, $(MAIN:.c=.o))

TEST_TARGET= $(addprefix $(BUILD_DIR)/, $(TEST:.c=))
MAIN_TARGET= $(addprefix $(BUILD_DIR)/, $(MAIN:.c=))

all: $(BUILD_DIR) $(MAIN_TARGET) $(TEST_TARGET)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(MAIN_TARGET): $(MAIN_OBJ) $(SRC_OBJ)
	$(CC) $(FLAGS) $^ -o $@

$(TEST_TARGET): %: %.o $(SRC_OBJ)
	$(CC) $(FLAGS) $^ -o $@

$(MAIN_OBJ): $(MAIN)
	$(CC) $(FLAGS) $^ -o $@ -c

$(TEST_OBJ): $(BUILD_DIR)/%.o: $(TEST_DIR)/%.c
	$(CC) $(FLAGS) $< -o $@ -c

$(SRC_OBJ): $(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(FLAGS) $< -o $@ -c

.PHONY:clean test perf
clean:
	rm -rf build results_perf

test:
	./build/test_polynomial 193 "polynomial"
	./build/test_rs_code 193 "rs_code"
	./build/test_compare 193 "compare"

perf:
	mkdir -p results_perf
	./build/test_perf