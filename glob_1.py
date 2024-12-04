def get_blosum62():
    return {
        'A': {'A': 4, 'R': -1, 'N': -2, 'D': -2, 'C': 0, 'Q': -1, 'E': -1, 'G': 0, 'H': -2, 'I': -1, 'L': -1, 'K': -1, 'M': -1, 'F': -2, 'P': 1, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': 0},
        'R': {'A': -1, 'R': 5, 'N': 0, 'D': 0, 'C': -2, 'Q': 1, 'E': 0, 'G': -1, 'H': -1, 'I': 0, 'L': -1, 'K': 2, 'M': 1, 'F': -1, 'P': -2, 'S': 0, 'T': -1, 'W': -3, 'Y': -2, 'V': -1},
        'N': {'A': -2, 'R': 0, 'N': 6, 'D': 1, 'C': -3, 'Q': 0, 'E': 0, 'G': 0, 'H': 1, 'I': -1, 'L': -1, 'K': 0, 'M': 0, 'F': -1, 'P': -2, 'S': 1, 'T': 0, 'W': -2, 'Y': -2, 'V': -1},
        'D': {'A': -2, 'R': 0, 'N': 1, 'D': 6, 'C': -3, 'Q': 0, 'E': 2, 'G': 1, 'H': 0, 'I': -1, 'L': -1, 'K': 0, 'M': 0, 'F': -1, 'P': -2, 'S': 0, 'T': 1, 'W': -3, 'Y': -2, 'V': -1},
        'C': {'A': 0, 'R': -2, 'N': -3, 'D': -3, 'C': 9, 'Q': -3, 'E': -3, 'G': -2, 'H': -3, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': -2, 'P': -2, 'S': 0, 'T': -1, 'W': -2, 'Y': -2, 'V': 0},
        'Q': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -3, 'Q': 5, 'E': 2, 'G': -1, 'H': 0, 'I': -1, 'L': -1, 'K': 1, 'M': 1, 'F': -1, 'P': -1, 'S': 0, 'T': 0, 'W': -2, 'Y': -2, 'V': -1},
        'E': {'A': -1, 'R': 0, 'N': 0, 'D': 2, 'C': -3, 'Q': 2, 'E': 5, 'G': 2, 'H': -1, 'I': -1, 'L': -1, 'K': 1, 'M': 1, 'F': -2, 'P': -1, 'S': 0, 'T': 0, 'W': -2, 'Y': -2, 'V': -1},
        'G': {'A': 0, 'R': -1, 'N': 0, 'D': 1, 'C': -2, 'Q': -1, 'E': 2, 'G': 6, 'H': -2, 'I': -2, 'L': -2, 'K': -1, 'M': -2, 'F': -2, 'P': -2, 'S': 0, 'T': -1, 'W': -2, 'Y': -2, 'V': -1},
        'H': {'A': -2, 'R': -1, 'N': 1, 'D': 0, 'C': -3, 'Q': 0, 'E': -1, 'G': -2, 'H': 8, 'I': -2, 'L': -2, 'K': 0, 'M': 0, 'F': -1, 'P': -2, 'S': 0, 'T': 1, 'W': -2, 'Y': -1, 'V': -1},
        'I': {'A': -1, 'R': 0, 'N': -1, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': 4, 'L': 2, 'K': 0, 'M': 1, 'F': 0, 'P': -1, 'S': -1, 'T': 0, 'W': -1, 'Y': 3, 'V': 3},
        'L': {'A': -1, 'R': -1, 'N': -1, 'D': -1, 'C': -1, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': 2, 'L': 4, 'K': 1, 'M': 2, 'F': 0, 'P': -1, 'S': -1, 'T': 0, 'W': -1, 'Y': 1, 'V': 1},
        'K': {'A': -1, 'R': 2, 'N': 0, 'D': 0, 'C': -2, 'Q': 1, 'E': 1, 'G': -1, 'H': 0, 'I': 0, 'L': 1, 'K': 5, 'M': 1, 'F': -2, 'P': -2, 'S': 0, 'T': 0, 'W': -2, 'Y': -2, 'V': -1},
        'M': {'A': -1, 'R': 1, 'N': 0, 'D': 0, 'C': -1, 'Q': 1, 'E': 1, 'G': -2, 'H': 0, 'I': 1, 'L': 2, 'K': 1, 'M': 5, 'F': 0, 'P': -1, 'S': -1, 'T': -1, 'W': -1, 'Y': -1, 'V': 1},
        'F': {'A': -2, 'R': -1, 'N': -1, 'D': -1, 'C': -2, 'Q': -1, 'E': -2, 'G': -2, 'H': -1, 'I': 0, 'L': 0, 'K': -2, 'M': 0, 'F': 6, 'P': -2, 'S': -2, 'T': -2, 'W': 2, 'Y': 3, 'V': 0},
        'P': {'A': 1, 'R': -2, 'N': -2, 'D': -2, 'C': -2, 'Q': -1, 'E': -1, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': -2, 'P': 7, 'S': -1, 'T': -1, 'W': -2, 'Y': -1, 'V': -1},
        'S': {'A': 0, 'R': 0, 'N': 1, 'D': 0, 'C': 0, 'Q': 0, 'E': 0, 'G': 0, 'H': 0, 'I': -1, 'L': -1, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 4, 'T': 1, 'W': -2, 'Y': -2, 'V': -1},
        'T': {'A': -1, 'R': -1, 'N': 0, 'D': 1, 'C': -1, 'Q': 0, 'E': 0, 'G': -1, 'H': 1, 'I': 0, 'L': 0, 'K': 0, 'M': -1, 'F': -2, 'P': -1, 'S': 1, 'T': 5, 'W': -2, 'Y': -2, 'V': 0},
        'W': {'A': -3, 'R': -3, 'N': -2, 'D': -3, 'C': -2, 'Q': -2, 'E': -2, 'G': -2, 'H': -2, 'I': -1, 'L': -1, 'K': -2, 'M': -1, 'F': 2, 'P': -2, 'S': -2, 'T': -2, 'W': 11, 'Y': 2, 'V': -3},
        'Y': {'A': -2, 'R': -2, 'N': -2, 'D': -2, 'C': -2, 'Q': -2, 'E': -2, 'G': -2, 'H': -1, 'I': 3, 'L': 1, 'K': -2, 'M': -1, 'F': 3, 'P': -1, 'S': -2, 'T': -2, 'W': 2, 'Y': 7, 'V': -1},
        'V': {'A': 0, 'R': -1, 'N': -1, 'D': -1, 'C': 0, 'Q': -1, 'E': -1, 'G': -1, 'H': -1, 'I': 3, 'L': 1, 'K': -1, 'M': 1, 'F': 0, 'P': -1, 'S': -1, 'T': 0, 'W': -3, 'Y': -1, 'V': 4}
    }
 

def global_alignment_score(v, w, scoring_matrix, sigma):
    '''Return the global alignment score of v and w subject to the given scoring matrix and indel penalty sigma.'''

    S = [[0 for j in range(len(w) + 1)] for i in range(len(v) + 1)]

    for i in range(1, len(v) + 1):
        S[i][0] = -i * sigma
    for j in range(1, len(w) + 1):
        S[0][j] = -j * sigma

    for i in range(1, len(v) + 1):
        for j in range(1, len(w) + 1):
            scores = [
                S[i - 1][j] - sigma,
                S[i][j - 1] - sigma,
                S[i - 1][j - 1] + scoring_matrix[v[i - 1]][w[j - 1]]
            ]
            S[i][j] = max(scores)

    return S[len(v)][len(w)]

def read_fasta(file_path):
    sequences = []
    current_sequence = ''

    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line.startswith('>'): 
                if current_sequence:
                    sequences.append(current_sequence)
                    current_sequence = ''
            else:
                current_sequence += line
        if current_sequence: 
            sequences.append(current_sequence)

    return sequences

if __name__ == '__main__':
    fasta_file_path = "rosalind_glob.txt" 
    sequences = read_fasta(fasta_file_path)
    if len(sequences) < 2:
        print("Error: The FASTA file must contain at least two sequences.")
        exit(1)

    s = sequences[0]
    t = sequences[1]

    scoring_matrix = get_blosum62()
    sigma = 5 
    score = global_alignment_score(s, t, scoring_matrix, sigma)

    print(score)