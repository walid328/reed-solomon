CC= gcc
FLAGS= -Wall
OBJFILES= main.o finite_field.o
TARGET= main

all: $(TARGET)

$(TARGET): $(OBJFILES)
	$(CC) $(FLAGS) $^ -o $@

%.o : %.c
	$(CC) $(FLAGS) -c $^ -o $@

.PHONY:clean
clean:
	rm -r $(OBJFILES) $(TARGET)