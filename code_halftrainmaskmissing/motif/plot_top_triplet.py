import pandas as pd
import matplotlib.pyplot as plt

# Define file paths
whole_genome_file = "top_100000_windows.txt.triplet"
mask_function_file = "../motif_maskfunctional/top_100000_windows.txt.triplet"
mask_label_file = "../motif_masklabel/top_100000_windows.txt.triplet"

# Function to read triplet data from a file
def read_triplet_file(file_path):
    triplets = []
    frequencies = []
    with open(file_path, "r") as file:
        for line in file:
            if line.strip():  # Ignore empty lines
                triplet, freq = line.strip().split(":")
                triplets.append(triplet.strip())
                frequencies.append(float(freq.strip()))
    return pd.DataFrame({"Triplet": triplets, "Frequency": frequencies})

# Read the data from the files
whole_genome_df = read_triplet_file(whole_genome_file)
mask_function_df = read_triplet_file(mask_function_file)
mask_label_df = read_triplet_file(mask_label_file)

# Find the top 10 triplets from the whole genome
top_10_triplets = whole_genome_df.nlargest(10, "Frequency")

# Extract corresponding frequencies from the other files
top_10_triplets["Mask Function"] = top_10_triplets["Triplet"].map(
    dict(zip(mask_function_df["Triplet"], mask_function_df["Frequency"]))
)
top_10_triplets["Mask Label"] = top_10_triplets["Triplet"].map(
    dict(zip(mask_label_df["Triplet"], mask_label_df["Frequency"]))
)

# Plot the top 10 triplets
def plot_triplet_frequencies(df, title, save_path):
    triplets = df["Triplet"]
    x = range(len(triplets))  # x positions for triplets

    plt.figure(figsize=(6, 4))
    plt.plot(x, df["Frequency"], marker="o", label="Whole genome model")
    plt.plot(x, df["Mask Function"], marker="o", label="Mask functional regions")
    plt.plot(x, df["Mask Label"], marker="o", label="Mask all labels")

    # Formatting
    plt.xticks(x, triplets, rotation=45)
    plt.ylabel("Frequency")
    plt.title(title)
    plt.legend()
    plt.tight_layout()

    # Save the figure
    plt.savefig(save_path)
    print(f"Figure saved to {save_path}")
    plt.close()

# Plot the top 10 triplets and save the figure
output_file = "top_10_triplet_frequencies.png"
plot_triplet_frequencies(top_10_triplets, "Top 10 Triplets and Their Frequency Changes", output_file)

