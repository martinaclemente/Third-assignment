def find_spliced_motif(s, t):
    position = [0]
    for i in t:
        position.append(s[position[-1]:].index(i) + 1 + position[-1])
    for p in position[1:]:
        print(p, end=" ")

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
    fasta_file_path = "rosalind_sseq.txt"
    sequences = list(read_fasta(fasta_file_path))
    
    if len(sequences) >= 2:
        s, t = sequences[0][1], sequences[1][1] 
        find_spliced_motif(s, t)