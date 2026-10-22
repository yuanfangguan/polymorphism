import numpy as np
import pickle
import os
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import mean_squared_error

# Step 1: Load the training predictions, labels, and genome data
def load_data(pred_file, label_file, genome_file):
    # Load predictions
    predictions = np.load(pred_file)

    # Load labels
    with open(label_file, 'rb') as f:
        labels = pickle.load(f)

    # Load genome (nX4 matrix)
    with open(genome_file, 'rb') as f:
        genome_data = pickle.load(f)  # nX4 matrix of binarized genome

    # Ensure predictions, labels, and genome_data have the same length by trimming
    print(len(predictions), len(labels), len(genome_data))
    min_length = min(len(predictions), len(labels), len(genome_data))
    predictions = predictions[:min_length]
    labels = np.array(labels[:min_length])
    genome_data = genome_data[:min_length]

    print(f"Loaded data: {len(predictions)} samples")
    return predictions, labels, genome_data

# Step 2: Exclude regions based on consecutive zeros in labels and all-zero genome positions
def filter_data(predictions, labels, genome_data, window_size=512, step_size=256):
    valid_indices = set()  # Use a set to prevent duplicate indices

    # Sliding window exclusion for consecutive zero labels
    for i in range(0, len(labels) - window_size + 1, step_size):
        window_labels = labels[i:i + window_size]
        window_genome = genome_data[i:i + window_size]

        # Exclude if all labels in the window are zero or if the genome is all-zero (nX4)
        if not (np.all(window_labels == 0) or np.any(np.all(window_genome == 0, axis=1))):
            valid_indices.update(range(i, i + window_size))  # Update set with unique indices

    # Convert set to sorted list to preserve order
    valid_indices = sorted(valid_indices)

    # Filter predictions and labels based on valid indices
    predictions = predictions[valid_indices]
    labels = labels[valid_indices]

    print(f"Filtered down to {len(predictions)} valid samples after exclusions.")
    return predictions, labels


# Step 3: Train an Isotonic Regression model for calibration
def train_calibration_model(predictions, labels):
    iso_reg = IsotonicRegression(out_of_bounds='clip')
    iso_reg.fit(predictions, labels)
    calibrated_preds = iso_reg.predict(predictions)

    # Print training RMSE to check calibration quality
    rmse = np.sqrt(mean_squared_error(labels, calibrated_preds))
    print(f"Calibration RMSE: {rmse:.5f}")
    return iso_reg

# Step 4: Apply calibration to other prediction files
def calibrate_predictions(pred_file, calibration_model):
    predictions = np.load(pred_file)
    calibrated_predictions = calibration_model.predict(predictions)
    return calibrated_predictions

# Step 5: Save calibrated predictions
def save_calibrated_predictions(predictions, output_file):
    np.save(output_file, predictions)
    print(f"Saved calibrated predictions to {output_file}")

# Main function to execute calibration and apply it to other predictions
def main():
    # Paths to training data for calibration
    #NC_000002.12.pred.npy
    train_pred_file = "train_chr/NC_000001.11.pred.npy"
    train_label_file = "../../data/NC_000001.11_label.pkl"
    train_genome_file = "../../data/NC_000001.11.pkl"

    # Load the training data and genome data
    train_predictions, train_labels, train_genome_data = load_data(train_pred_file, train_label_file, train_genome_file)

    # Filter data to exclude regions with all-zero labels and genome rows with all-zero values
    filtered_predictions, filtered_labels = filter_data(train_predictions, train_labels, train_genome_data)

    # Train the calibration model
    calibration_model = train_calibration_model(filtered_predictions, filtered_labels)

    # Directory containing other predictions to be calibrated
    prediction_dir = "./"
    output_dir = "/local/disk4/gyuanfan/rate/calibrated_predictions"

    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)

    # Calibrate and save predictions for each file
    for filename in os.listdir(prediction_dir):
        if filename.endswith(".npy"):
            pred_file = os.path.join(prediction_dir, filename)
            calibrated_preds = calibrate_predictions(pred_file, calibration_model)

            # Save the calibrated predictions
            output_file = os.path.join(output_dir, filename)
            save_calibrated_predictions(calibrated_preds, output_file)

if __name__ == "__main__":
    main()

