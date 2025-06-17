# Knockoff-augmented neural networks for identifying risk variants in family-based association studies

**KNOT** is a **k**nockoff-augmented **n**eural network **o**n **t**rio data for stabilized variable selection with false discovery rate control in family-based genome-wide association studies.

The repository includes scripts for simulation data generation, model training, and feature importance computation.

![Pipeline](figure/framework.jpg)


## Repository Structure

- `function/model_combine.py`: Defines the DNN architecture with `PositionalEncoding` for sequence position encoding, `LocallyConnected1D` for non-shared weight convolutions, and `DNN` for the full model, including a Siamese encoder and classifier.
- `function/utils.py`: Contains the `Args` class for managing hyperparameters (e.g., learning rate, epochs, latent dimension) and configuration settings.
- `function/run_multiple_server.py`: Main script to execute experiments, handling data loading, model training, and feature importance computation using SHAP and gradients.
- `function/callback_prediction_quan_combine.py`: Trainer class for quantitative tasks, implementing MSE loss, distance loss, L1 regularization, and R² validation.
- `function/callback_prediction_combine.py`: Trainer class for classification tasks, using BCE loss, distance loss, L1 regularization, and ROC-AUC validation.


## Installation

1. **Clone the Repository**:

   ```bash
   git clone https://github.com/EddieFua/KNOT.git
   cd KNOT
   ```

2. **Set Up Environment**:

   ```bash
   python -m venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -r requirements.txt
   ```

   **Required Dependencies**:

   - `torch`
   - `numpy`
   - `pandas`
   - `pyreadr`
   - `shap`
   - `scikit-learn`
   - `tqdm`

## Usage

### Running an Experiment

Use `run_multiple_server.py` to run experiments:

```bash
python run_multiple_server.py --sample_size 3000 --quan False --data_path /path/to/data
```

#### Arguments:

- `--sample_size`: Number of samples (default: 3000).
- `--quan`: True for quantitative, False for classification (default: False).
- `--data_path`: Directory with data files (e.g., `child_array.RData`, `dad_array.RData`, `mom_array.RData`, `weight.csv`, `y.RData`).

#### Data Format:

- **Genomic Data**: `.RData` files with shape `(num_samples, num_features, num_knockoffs)`.
- **Weights**: `weight.csv` with shape `(num_knockoffs + 1, num_features)`.
- **Labels**: Quantitative tasks use `y.RData` (continuous values); classification tasks generate binary tensors.

#### Outputs:

- `FI_nn_final_gradient.csv`: Gradient-based feature importance.
- `FI_nn_final_shap.csv`: SHAP-based feature importance.

### Example Workflow

1. **Prepare Data**: Use the R code `R/10_replicate.R` to generate simulation data.

2. **Run Experiment**:

   ```bash
   python run_multiple_server.py --sample_size 5000 --quan True --data_path ./data
   ```

3. **Analyze Results**: Feature importance scores are saved as CSV files for downstream analysis. Use the R script `R/control_fdr_function.R` to analyze these results under knockoff FDR control.

## Model Details

### Training

- **Early Stopping**: Patience of 5 epochs.
- **Loss**:
  - Quantitative: MSE + contrastive loss + L1 regularization.
  - Classification: BCE + contrastive loss + L1 regularization.
- **Hyperparameters**: Configurable via `Args` (e.g., learning rate: 0.0001, epochs: 50, latent_dim: 32).

## Applications

- **GWAS**: Identifies significant features using p-values.
- **Pathway Enrichment**: Maps features to biological pathways.
- **PRS**: Computes risk scores based on feature effect sizes.

## Contact

For questions, open an issue on GitHub or [email](yinghao.fu@my.cityu.edu.hk), or visit my [personal homepage](https://eddiefua.github.io/).
