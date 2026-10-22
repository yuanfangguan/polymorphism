from Bio import Phylo
from Bio.Phylo.PAML import codeml
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

def calculate_dn_ds(seq1, seq2):
    # Convert sequences to Biopython Seq objects
    seq1 = Seq(seq1)
    seq2 = Seq(seq2)

    # Create a MultipleSeqAlignment object with the two sequences
    seq_record1 = SeqRecord(seq1, id="seq1")
    seq_record2 = SeqRecord(seq2, id="seq2")
    alignment = MultipleSeqAlignment([seq_record1, seq_record2])

    # Calculate the distance matrix
    calculator = DistanceCalculator('identity')
    distance_matrix = calculator.get_distance(alignment)

    # Construct a phylogenetic tree using UPGMA method
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(distance_matrix)

    # Save the tree in Newick format (required by codeml)
    Phylo.write(tree, "tree.nwk", "newick")

    # Prepare the input control file for codeml
    with open("codeml.ctl", "w") as ctl_file:
        ctl_file.write(f"seqfile = seqfile.txt\n")
        ctl_file.write(f"treefile = tree.nwk\n")
        ctl_file.write(f"outfile = yn00.out\n")

    # Prepare the sequence file for codeml
    with open("seqfile.txt", "w") as seq_file:
        seq_file.write(f">{seq_record1.id}\n{seq_record1.seq}\n")
        seq_file.write(f">{seq_record2.id}\n{seq_record2.seq}\n")

    # Run codeml and calculate the dn/ds ratio
#    codeml.run(command="codeml", control_file="codeml.ctl")

    # Parse the codeml results to obtain dn and ds
    with open("yn00.out") as result_file:
        for line in result_file:
            if "dN/dS=" in line:
                dn_ds_ratio = float(line.split("dN/dS=")[-1])
                break

    return dn_ds_ratio

# Example usage
sequence1 = "ATGCTAGCGTAACTG"
sequence2 = "ATGCTAGAGTAACAG"
result = calculate_dn_ds(sequence1, sequence2)
print(f"dn/ds ratio: {result}")

