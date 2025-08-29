# Knockoff-Augmented Neural Networks for Identifying Risk Variants in Family-Based Association Studies

**KNOT** (**K**nockoff-augmented **N**eural network **o**n **T**rio data) is designed for stabilized variable selection with false discovery rate (FDR) control in family-based genome-wide association studies (GWAS).

This repository provides scripts for simulation data generation, model training, and feature importance computation.

![Pipeline](figure/framework.jpg)

---

## Repository Structure

- `KNOT/model_combine.py`: Defines the DNN architecture with:
  - `PositionalEncoding` for sequence position encoding  
  - `LocallyConnected1D` for non-shared weight convolutions  
  - `DNN` for the full model, including a Siamese encoder and classifier  
- `KNOT/utils.py`: Contains the `Args` class for managing hyperparameters (e.g., learning rate, epochs, latent dimension) and configuration settings.  
- `KNOT/run_multiple_server.py`: Main script to execute experiments, handling data loading, model training, and feature importance computation using SHAP and gradients.  
- `KNOT/callback_prediction_quan_combine.py`: Trainer class for quantitative tasks, implementing MSE loss, contrastive (distance) loss, L1 regularization, and R² validation.  
- `KNOT/callback_prediction_combine.py`: Trainer class for classification tasks, using BCE loss, contrastive (distance) loss, L1 regularization, and ROC-AUC validation.  
- `KNOT/generate_knockoffs.R`: Functions for generating knockoffs.  
- `KNOT/permutation_test.R`: Function for permutation tests based on SHAP interaction values.  

---

## Installation

1. **Clone the repository**:

   ```bash
   git clone https://github.com/EddieFua/KNOT.git
   cd KNOT

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

### Generating knockoffs:
Use `generate_knockoffs.R` to create knockoffs.
The file `original.RData` contains simulated genotype data with samples ordered as: dad → mom → offspring in each trio.
```R
load("./example_data/Binary/original.RData")
dat1 = knockofftrio_create_knockoff(
  dat = sim$dat,
  pos = sim$pos,
  M = 10,
  hap = TRUE,
  dat.hap = sim$dat.hap,
  xchr = FALSE,
  sex = sim$sex,
  phasing.dad = sim$phasing.dad,
  phasing.mom = sim$phasing.mom
)
```


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
- **Interaction Identification**: Permutation test based on SHAP interaction value.
- **PRS**: Computes risk scores based on feature effect sizes.

```R
### Example: Interaction Identification
# genotype_dat: matrix with rows as samples and columns as selected variants
# y: phenotype labels
# gene_closer: mapping of each variant to its gene
# N: number of permutations

res = compute_shap_interaction_pvalues(
  genotype_dat,
  y,
  gene_closer,
  N = 100,
  cores = 10,
  seed = 10
)
```

## Contact

For questions, open an issue on GitHub or [email](yinghao.fu@my.cityu.edu.hk), or visit my [personal homepage](https://eddiefua.github.io/).
