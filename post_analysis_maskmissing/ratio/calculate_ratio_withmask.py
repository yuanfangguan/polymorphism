import os
import numpy as np
import pickle

def load_predictions(pred_file):
    """
    Load predictions from a .npy file.
    """
    return np.load(pred_file)

def load_labels(label_file):
    """
    Load true labels from a .pkl file.
    """
    with open(label_file, 'rb') as f:
        labels = pickle.load(f)
    return labels

def load_genome_matrix(genome_file):
    """
    Load the genome encoding matrix from a .pkl file.
    """
    with open(genome_file, 'rb') as f:
        genome_matrix = pickle.load(f)
    return genome_matrix

def calculate_ratio(predictions, labels, genome_matrix):
    """
    Calculate the ratio of predictions to labels. If the genome row is all-zero,
    set the ratio to -1.
    """
    ratios = np.zeros_like(predictions, dtype=np.float32)

    for i in range(len(predictions)):
        if np.all(genome_matrix[i] == 0):  # Check if the genome row is all-zero
            ratios[i] = -1
        else:
            if labels[i] != 0:
                ratios[i] = labels[i] / predictions[i]
            else:
                ratios[i] = np.inf  # Avoid division by zero; set to infinity

    return ratios

def process_files(pred_dir, label_dir, genome_dir, output_dir):
    """
    Process all files and generate the ratio files.
    """
    os.makedirs(output_dir, exist_ok=True)

    for pred_file in sorted(os.listdir(pred_dir)):
        if pred_file.endswith(".pred.npy"):
            # Construct corresponding label and genome file paths
            base_name = pred_file.replace(".pred.npy", "")
            label_file = os.path.join(label_dir, f"{base_name}_label.pkl")
            genome_file = os.path.join(genome_dir, f"{base_name}.pkl")

            # Load data
            predictions = load_predictions(os.path.join(pred_dir, pred_file))
            labels = load_labels(label_file)
            genome_matrix = load_genome_matrix(genome_file)

            # Calculate ratio
            ratios = calculate_ratio(predictions, labels, genome_matrix)

            # Save the result
            output_file = os.path.join(output_dir, f"{base_name}_ratio.npy")
            np.save(output_file, ratios)
            print(f"Saved: {output_file}")

# Example usage
pred_dir = "../../calibrated_predictions/"
label_dir = "../../data/"
genome_dir = "../../data/"
output_dir = "./"

process_files(pred_dir, label_dir, genome_dir, output_dir)

