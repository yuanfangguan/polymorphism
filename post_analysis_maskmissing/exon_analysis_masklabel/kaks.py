from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo

def calculate_ka_ks(seq1, seq2):
    # Define a function to count synonymous and nonsynonymous mutations
    def count_mutations(codon1, codon2):
        # Initialize counters
        synonymous = 0
        nonsynonymous = 0
        # Ensure codons are of the same length
        min_len = min(len(codon1), len(codon2))
        # Loop through each position in the codons
        for i in range(min_len):
            if codon1[i] != codon2[i]:
                # If the codons differ
                if Seq(codon1[i:i+3]).translate() == Seq(codon2[i:i+3]).translate():
                    # If the mutations are synonymous
                    synonymous += 1
                else:
                    # If the mutations are nonsynonymous
                    nonsynonymous += 1
        return synonymous, nonsynonymous

    # Initialize counts
    total_synonymous = 0
    total_nonsynonymous = 0

    # Iterate over codons in the aligned sequences
    for i in range(0, len(seq1), 3):
        codon1 = seq1[i:i+3]
        codon2 = seq2[i:i+3]
        synonymous, nonsynonymous = count_mutations(codon1, codon2)
        total_synonymous += synonymous
        total_nonsynonymous += nonsynonymous

    # Calculate Ka and Ks
    try:
        ka = total_nonsynonymous / (len(seq1) // 3)
        ks = total_synonymous / (len(seq1) // 3)
        ka_ks_ratio = ka / ks
        return ka_ks_ratio,ka,ks
    except ZeroDivisionError:
        return None, ka,ks

# Example sequences
seq1 = "ATGCGTACGATGAGACCTAG"
seq2 = "ATGCGTACGATGAGACCGAG"

# Calculate Ka/Ks ratio
ka_ks_ratio,ka,ks = calculate_ka_ks(seq1, seq2)
print("Ka/Ks ratio:", ka_ks_ratio,ka,ks)

