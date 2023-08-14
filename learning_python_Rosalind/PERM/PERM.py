#Given: A positive integer n≤7

#Return: The total number of permutations of length n, followed by a list of all such permutations (in any order).

import itertools

def generate_permutations(n):
    permutations = list(itertools.permutations(range(1, n + 1)))
    return permutations

def main():
    n = int(input())
    permutations = generate_permutations(n)

    total_permutations = len(permutations)
    print(total_permutations)

    for p in permutations:
        print(" ".join(str(x) for x in p))

if __name__ == "__main__":
    main()

