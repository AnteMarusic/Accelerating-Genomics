import random

# Define the minimum and maximum length for the random strings
MIN_LEN = 450
MAX_LEN = 500

# Function to generate a random string of A, T, G, C of random length between MIN_LEN and MAX_LEN
def generate_random_string():
    length = random.randint(MIN_LEN, MAX_LEN)  # Use randint for inclusive range
    characters = "ATGC"
    random_string = ""
    for _ in range(length):
        random_char = random.choice(characters)
        random_string += random_char
    return random_string

# Define the number of alignments (replace with your desired number)
NUM_OF_ALIGNMENTS = 500
print("ciao")
# Open the output file in write mode
with open("input.txt", "w") as output_file:
    output_file.write(str(NUM_OF_ALIGNMENTS) + "\n")
    # Generate and write 2 * NUM_OF_ALIGNMENTS strings to the file
    for _ in range(2 * NUM_OF_ALIGNMENTS):
        random_string = generate_random_string()
        output_file.write(random_string + "\n")
