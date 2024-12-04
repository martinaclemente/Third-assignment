def longestCommonSubsequence(s, t):
    m, n = len(s), len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    
    for i, x in enumerate(s):
        for j, y in enumerate(t):
            if x == y:
                dp[i + 1][j + 1] = 1 + dp[i][j]
            else:
                dp[i + 1][j + 1] = max(dp[i + 1][j], dp[i][j + 1])
    
    longest_common_subsequence = ""
    while m * n != 0:
        if dp[m][n] == dp[m - 1][n]:
            m -= 1
        elif dp[m][n] == dp[m][n - 1]:
            n -= 1
        else:
            longest_common_subsequence += s[m - 1]
            m -= 1
            n -= 1
    
    return longest_common_subsequence[::-1]

def read_fasta(file_path):
    with open(file_path, 'r') as fa:
        seq_name = None
        seq_string = ""
        for line in fa:
            line = line.strip()
            if line.startswith(">"):
                if seq_name is not None:
                    yield seq_name, seq_string
                seq_name = line[1:] 
                seq_string = ""
            else:
                seq_string += line
        if seq_name is not None:
            yield seq_name, seq_string



if __name__ == "__main__":
    fasta_file_path = "rosalind_lcsq.txt"  
    sequences = list(read_fasta(fasta_file_path))
    
    if len(sequences) >= 2:
        s, t = sequences[0][1], sequences[1][1] 
        lcs_result = longestCommonSubsequence(s, t)
        print(lcs_result)