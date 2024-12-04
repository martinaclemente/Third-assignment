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
def edit_distance(s, t):
    m, n = len(s), len(t)
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    for i in range(m + 1):
        dp[i][0] = i 
    for j in range(n + 1):
        dp[0][j] = j
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if s[i - 1] == t[j - 1]:
                dp[i][j] = dp[i - 1][j - 1] 
            else:
                dp[i][j] = min(
                    dp[i - 1][j] + 1,  
                    dp[i][j - 1] + 1,   
                    dp[i - 1][j - 1] + 1 
                )
    
    return dp[m][n]

if __name__ == "__main__":
    fasta_file_path = "rosalind_edit.txt" 
    sequences = list(read_fasta(fasta_file_path))
    
    if len(sequences) >= 2:
        s, t = sequences[0][1], sequences[1][1]  
        edit_distance_result = edit_distance(s, t)
        print(edit_distance_result)