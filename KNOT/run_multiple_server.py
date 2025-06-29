#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import torch
import numpy as np
import pandas as pd
import pyreadr
from torch.utils.data import DataLoader, TensorDataset, random_split
from model_combine import DNN
from utils import Args          
import torch.optim as optim
from tqdm.auto import tqdm  
import shap
from torch.utils.data import random_split

def run_experiment(sample_size, quan, data_path, seed=42):
    """Run a single experiment for genomic data analysis.
    
    Args:
        sample_size: Number of samples to use
        quan: Boolean indicating quantitative (True) or classification (False) task
        data_path: Path to data directory
        seed: Random seed for reproducibility
    """
    torch.manual_seed(seed)
    np.random.seed(seed)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    
    if quan:
        from callback_prediction_quan_combine import TrainerWithCallback
        alpha = 50
    else:
        from callback_prediction_combine import TrainerWithCallback
        alpha = 8

    # Load data
    child_array = pyreadr.read_r(f'{data_path}/child_array.RData')['child_array'].values
    dad_array = pyreadr.read_r(f'{data_path}/dad_array.RData')['dad_array'].values
    mom_array = pyreadr.read_r(f'{data_path}/mom_array.RData')['mom_array'].values
    weight = pd.read_csv(f'{data_path}/weight.csv').to_numpy()
    weight_tensor = torch.tensor(weight, dtype=torch.float32)

    num_samples, num_features, num_knockoffs = child_array.shape[:3]
    merged_array = np.stack((dad_array, mom_array, child_array), axis=-1)
    X_tensor = torch.tensor(merged_array[:sample_size], dtype=torch.float32)

    if quan:
        Y_data = pyreadr.read_r(f'{data_path}/y.RData')['y'].values
        Y_data = (Y_data - Y_data.mean()) / Y_data.std()
        Y_tensor = torch.tensor(Y_data[:sample_size], dtype=torch.float32)
    else:
        Y_tensor = torch.stack([torch.zeros(sample_size), torch.zeros(sample_size), torch.ones(sample_size)], dim=1)

    # Normalize weights
    args = Args(num_features, num_knockoffs - 1, lr=0.0001, num_epochs=50, device=device, weight=weight_tensor, quan=quan)
    for i in range(args.Num_knock + 1):
        col_min, col_max = args.weight[i].min().item(), args.weight[i].max().item()
        args.weight[i] = (args.weight[i] - col_min) / (col_max - col_min)

    # Prepare data loaders
    dataset = TensorDataset(X_tensor, Y_tensor)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_dataset, test_dataset = random_split(dataset, [train_size, test_size])

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False)

    model = DNN(args).to(device)
    optimizer = optim.Adam(model.parameters(), lr=args.lr, weight_decay=args.weight_decay)
    trainer = TrainerWithCallback(args, model, optimizer, device)
    trained_model, feature_importances = trainer.train(train_loader, test_loader, num_epochs=args.num_epochs, alpha=alpha)
    df = pd.DataFrame(feature_importances.T)
    df.to_csv(f'{data_path}/FI_nn_final_gradient.csv', index=False, header=False)


    test_data = torch.cat([X_b for X_b, Y_b in test_loader], dim=0)
    dad_tensor = test_data[:, :, :, 0].float()
    mom_tensor = test_data[:, :, :, 1].float()
    child_tensor_test = test_data[:, :, :, 2].float()
    mean_tensor_test = ((dad_tensor + mom_tensor) / 2)

    train_data = torch.cat([X_b for X_b, Y_b in train_loader], dim=0)
    dad_tensor = train_data[:, :, :, 0].float()
    mom_tensor = train_data[:, :, :, 1].float()
    child_tensor_train = train_data[:, :, :, 2].float()
    mean_tensor_train = ((dad_tensor + mom_tensor) / 2)

    all_mean_tensor = torch.cat((mean_tensor_train, mean_tensor_test), dim=0)
    all_child_tensor = torch.cat((child_tensor_train, child_tensor_test), dim=0)
    trained_model.return_y_only = True
    trained_model.eval()
    torch.cuda.empty_cache()
    all_mean_tensor = torch.cat((mean_tensor_train, mean_tensor_test), dim=0)
    all_child_tensor = torch.cat((child_tensor_train, child_tensor_test), dim=0)
    shap_values_list = []
    batch_size = all_child_tensor.shape[0]
    print("开始按一一配对计算 SHAP 值...")
    for i in tqdm(range(0, all_child_tensor.shape[0], batch_size), desc="Batch SHAP"):
        batch_start = i
        batch_end = min(i + batch_size, all_child_tensor.shape[0])
        background = all_mean_tensor[batch_start:batch_end].to(device)  # [batch_size, num_features, num_knockoffs]
        explanation = all_child_tensor[batch_start:batch_end].to(device)  # [batch_size, num_features, num_knockoffs]
        explainer = shap.GradientExplainer(trained_model, background)
        shap_val = explainer.shap_values(explanation)
        if isinstance(shap_val, list):
            shap_val = shap_val[0]
        shap_values_list.append(shap_val)
    shap_values_all = np.concatenate(shap_values_list, axis=0)
    mean_abs_shap = np.mean(np.abs(shap_values_all), axis=0)
    # mean_abs_shap = np.mean(shap_values_all, axis=0)
    if mean_abs_shap.ndim > 2:
        mean_abs_shap = np.squeeze(mean_abs_shap, axis=-1)
    df = pd.DataFrame(mean_abs_shap.T)
    df.to_csv(f'{data_path}/FI_nn_final_shap.csv', index=False, header=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run genomic data analysis experiment.")
    parser.add_argument('--sample_size', type=int, default=3000, help="Number of samples")
    parser.add_argument('--quan', type=lambda x: x.lower() == 'true', default=False, help="Quantitative task flag")
    parser.add_argument('--data_path', type=str, default='/your_path', help="Path to data directory")
    args = parser.parse_args()
    run_experiment(args.sample_size, args.quan, args.data_path)
