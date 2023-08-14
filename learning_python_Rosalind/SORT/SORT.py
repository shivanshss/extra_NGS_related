#Given: Two permutations π and γ, each of length 10.

#Return: The reversal distance drev(π,γ), followed by a collection of reversals sorting π into γ. If multiple collections of such reversals exist, you may return any one.

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
    return len(pancake_sort(p1 + 1))

def main():
    p1 = list(map(int, input().split()))
    p2 = list(map(int, input().split()))

    distance = reversal_distance(p1, p2)
    reversals = pancake_sort(p1 + 1)

    print(distance)
    for rev in reversals:
        print(rev, end=" ")

if __name__ == "__main__":
    main()

