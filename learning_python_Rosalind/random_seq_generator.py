import random

def generate_random_dna_sequence(length):
    return ''.join(random.choice('ATGC') for _ in range(length))

s_length = 800  #change as needed
t_length = 1200  #change as needed

s = generate_random_dna_sequence(s_length)
t = generate_random_dna_sequence(t_length)

print("s:", s)
print("t:", t)
