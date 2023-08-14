#Given: Two protein strings s and t in FASTA format (with each string having length at most 1000 aa).

#Return: The edit distance dE(s,t) followed by two augmented strings s′ and t′ representing an optimal alignment of s and t.

def optimal_alignment(s, t):
    m = len(s)
    n = len(t)

    # Create a 2D table to store edit distances
    dp = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize the first row and column
    for i in range(m + 1):
        dp[i][0] = i
    for j in range(n + 1):
        dp[0][j] = j

    # Fill in the table using dynamic programming
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            cost = 0 if s[i - 1] == t[j - 1] else 1
            dp[i][j] = min(dp[i - 1][j] + 1, dp[i][j - 1] + 1, dp[i - 1][j - 1] + cost)

    # Backtrack to find the optimal alignment
    s_aligned = []
    t_aligned = []
    i, j = m, n
    while i > 0 and j > 0:
        if s[i - 1] == t[j - 1]:
            s_aligned.append(s[i - 1])
            t_aligned.append(t[j - 1])
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j - 1] + 1:
            s_aligned.append(s[i - 1])
            t_aligned.append(t[j - 1])
            i -= 1
            j -= 1
        elif dp[i][j] == dp[i - 1][j] + 1:
            s_aligned.append(s[i - 1])
            t_aligned.append("-")
            i -= 1
        else:
            s_aligned.append("-")
            t_aligned.append(t[j - 1])
            j -= 1

    while i > 0:
        s_aligned.append(s[i - 1])
        t_aligned.append("-")
        i -= 1
    while j > 0:
        s_aligned.append("-")
        t_aligned.append(t[j - 1])
        j -= 1

    return dp[m][n], "".join(reversed(s_aligned)), "".join(reversed(t_aligned))

def main():
    # Read input sequences from FASTA format
    sequences = []
    with open("input.fasta", "r") as fasta_file:
        lines = fasta_file.readlines()
        sequences = [line.strip() for line in lines[1:]]

    s = sequences[0]  # First protein string
    t = sequences[1]  # Second protein string

    # Find the optimal alignment
    distance, s_aligned, t_aligned = optimal_alignment(s, t)
    print(distance)
    print(s_aligned)
    print(t_aligned)

if __name__ == "__main__":
    main()

