#Given: A positive integer n≤6.

#Return: The total number of signed permutations of length n, followed by a list of all such permutations (you may list the signed permutations in any order).

from itertools import permutations

def signed_permutations(n):
    numbers = list(range(1, n + 1))
    signed_permuts = []

    for p in permutations(numbers):
        for signs in product([-1, 1], repeat=n):
            signed_permuts.append([str(p[i] * signs[i]) for i in range(n)])

    return signed_permuts

def main():
    n = int(input())
    signed_permuts = signed_permutations(n)
    total_permuts = len(signed_permuts)

    print(total_permuts)
    for permut in signed_permuts:
        print(" ".join(permut))

if __name__ == "__main__":
    main()

