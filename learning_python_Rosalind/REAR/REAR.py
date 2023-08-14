#Given: A collection of at most 5 pairs of permutations, all of which have length 10.

#Return: The reversal distance between each permutation pair.

def flip(arr, k):
    return arr[:k][::-1] + arr[k:]

def pancake_sort(arr):
    n = len(arr)
    result = []

    for target in range(n, 0, -1):
        target_index = arr.index(target)
        if target_index != target - 1:
            if target_index != 0:
                arr = flip(arr, target_index + 1)
                result.append(target_index + 1)
            arr = flip(arr, target)
            result.append(target)

    return result

def reversal_distance(p1, p2):
    return len(pancake_sort(p1 + 1)) + len(pancake_sort(p2 + 1))

def main():
    permutations = []
    while True:
        try:
            p1 = list(map(int, input().split()))
            p2 = list(map(int, input().split()))
            permutations.append((p1, p2))
        except EOFError:
            break
    
    for p1, p2 in permutations:
        distance = reversal_distance(p1, p2)
        print(distance, end=" ")

if __name__ == "__main__":
    main()

