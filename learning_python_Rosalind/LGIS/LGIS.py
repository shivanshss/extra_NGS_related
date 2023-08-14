#Given: A positive integer n≤10000 followed by a permutation π of length n

#Return: A longest increasing subsequence of π, followed by a longest decreasing subsequence of π.

def longest_increasing_subsequence(arr):
    n = len(arr)
    lis = [1] * n

    for i in range(1, n):
        for j in range(i):
            if arr[i] > arr[j] and lis[i] < lis[j] + 1:
                lis[i] = lis[j] + 1

    max_lis_length = max(lis)
    max_lis = []
    last_index = lis.index(max_lis_length)

    for i in range(n - 1, -1, -1):
        if lis[i] == max_lis_length and arr[i] <= arr[last_index]:
            max_lis.append(arr[i])
            max_lis_length -= 1
            last_index = i

    max_lis.reverse()
    return max_lis

def longest_decreasing_subsequence(arr):
    reverse_arr = arr[::-1]
    return longest_increasing_subsequence(reverse_arr)

def main():
    n = int(input())
    permutation = list(map(int, input().split()))

    increasing_subsequence = longest_increasing_subsequence(permutation)
    decreasing_subsequence = longest_decreasing_subsequence(permutation)

    print(" ".join(str(x) for x in increasing_subsequence))
    print(" ".join(str(x) for x in decreasing_subsequence))

if __name__ == "__main__":
    main()

