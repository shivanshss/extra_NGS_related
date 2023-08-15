#Given: Positive integers n≤40 and k≤5.

#Return: The total number of rabbit pairs that will be present after n months, if we begin with 1 pair and in each generation, every pair of reproduction-age rabbits produces a litter of k rabbit pairs (instead of only 1 pair).



def rabbit_pairs(n, k):
    if n == 1:
        return 1
    elif n == 2:
        return 1

    prev_prev = 1
    prev = 1
    current = 0

    for i in range(3, n + 1):
        current = prev + k * prev_prev
        prev_prev = prev
        prev = current

    return current

# Example usage
n, k = map(int, input().split())  # Replace with input values
result = rabbit_pairs(n, k)
print(result)

