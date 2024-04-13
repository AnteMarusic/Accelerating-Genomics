# Define the output executable name
CC = gcc  # Set your compiler here (e.g., clang)
CFLAGS = -Wall -g  # Compiler flags (warnings, debugging)
OUTPUT = program

# Define object files based on source files
OBJECTS = main.o int_matrix.o algo.o # Replace with your source file names

# Main rule: builds the executable
$(OUTPUT): $(OBJECTS)
	$(CC) $(CFLAGS) $(OBJECTS) -o $(OUTPUT)

# Object file rule: compiles a source file to an object file
$(OBJECTS): %.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Rule to clean up object files
clean:
	rm -f $(OBJECTS) $(OUTPUT)

# Phony target for running the program (replace with your actual command)
run: $(OUTPUT)
	./$(OUTPUT)

.PHONY: clean run  # Declare clean and run as phony targets
