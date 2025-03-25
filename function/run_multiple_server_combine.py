import pyreadr
import numpy as np
import pandas as pd
import torch
from model_combine import DNN
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm
import shap
import argparse
import os
import random
import xarray as xr
import json
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, random_split
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score, classification_report
from sklearn.metrics import mean_squared_error, r2_score

# 命令行参数解析
parser = argparse.ArgumentParser()
parser.add_argument('--replicate', type=int, default=0)
parser.add_argument('--sample_size', type=int, default=3000)
parser.add_argument('--linear', type=lambda x: x.lower() == 'true', default=True)
parser.add_argument('--quan', type=lambda x: x.lower() == 'true', default=False)
args = parser.parse_args()

sample_size = args.sample_size
linear = args.linear
quan = args.quan
replicate = args.replicate

device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

if quan:
    from callback_prediction_quan_combine import TrainerWithCallback
else:
    from callback_prediction_combine import TrainerWithCallback
def init_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    
class Args:
    def __init__(self, num_features, num_knockoffs, lr, num_epochs, device, weight, quan, weight_cnn, gate_init):
        self.device = device
        self.pVal = num_features
        self.FILTER = 12   ##nonlinear quan 8
        self.dropout = 0.3
        self.gate_init = gate_init
        self.Num_knock = num_knockoffs
        self.KERNEL = 3
        self.STRIDE = 2
        self.batch_size = 128
        self.num_epochs = num_epochs
        self.l1_regularization = 0.00002
        self.weight_decay = 0.01
        self.lr = lr
        self.lr_milestones = [50]
        self.latent_dim = 32
        self.weight = weight
        self.weight_cnn = weight_cnn
        self.bias = True
        self.quan = quan

def run_one_experiment(sample_size, linear, replicate, k, quan):
    if linear and quan:
        linear_or_not = 'linear'
        quan_or_not = 'quan_'
        alpha = 12
        epoch = 80
    elif not linear and quan:
        linear_or_not = 'nonlinear'
        quan_or_not = 'quan_'
        alpha = 60  #10
        epoch = 80
    elif linear and not quan:
        linear_or_not = 'linear'
        quan_or_not = ''
        alpha = 3
        epoch = 50
    else:
        linear_or_not = 'nonlinear'
        quan_or_not = ''
        alpha = 8
        epoch = 40 #80
    if quan:
        base_path = f'/gpfs1/home/yinghaofu2/scratch/simulation/sim_one_class_{linear_or_not}_no_control_{sample_size}_500_{quan_or_not}{replicate + 1}y'
    else:
        base_path = f'/gpfs1/home/yinghaofu2/scratch/simulation/sim_one_class_{linear_or_not}_no_control_{sample_size}_500_{quan_or_not}{replicate + 1}'
    child_array = pyreadr.read_r(f'{base_path}/child_array.RData')['child_array'].values
    dad_array = pyreadr.read_r(f'{base_path}/dad_array.RData')['dad_array'].values
    mom_array = pyreadr.read_r(f'{base_path}/mom_array.RData')['mom_array'].values

    num_samples, num_features, num_knockoffs = child_array.shape
    weight = pd.read_csv(f'{base_path}/weight.csv', header=0)
    weight = weight.apply(pd.to_numeric, errors='coerce').to_numpy()
    weight_tensor = torch.tensor(weight, dtype=torch.float32)

    reduced_dad_array = dad_array
    reduced_mom_array = mom_array
    reduced_child_array = child_array
    lr = 0.0001
    num_epochs = 50

    merged_array = np.stack((reduced_dad_array, reduced_mom_array, reduced_child_array), axis=-1)
    X_tensor = torch.tensor(merged_array, dtype=torch.float32)

    if quan:
        result = pyreadr.read_r(f'{base_path}/y.RData')
        Y_data = None
        for name, data in result.items():
            Y_data = data.values
        Y_tensor = torch.tensor(Y_data, dtype=torch.float32)
        # Standardize the targets
        Y_data = (Y_data - Y_data.mean()) / Y_data.std()
        Y_tensor = torch.tensor(Y_data, dtype=torch.float32)
    else:
        dad_y = torch.zeros(num_samples)
        mom_y = torch.zeros(num_samples)
        child_y = torch.ones(num_samples)
        Y_tensor = torch.stack((dad_y, mom_y, child_y), dim=1)

    args = Args(num_features, num_knockoffs - 1, lr, num_epochs, device, weight_tensor, quan, weight_cnn=weight_tensor, gate_init=20)
    args.weight = weight_tensor.clone()
    args.weight_cnn = weight_tensor.clone() 
    for i_kn in range(args.Num_knock + 1):
        col_min = args.weight[i_kn].min().item()
        col_max = args.weight[i_kn].max().item()
        args.weight[i_kn] = (args.weight[i_kn] - col_min) / (col_max - col_min)
    
    if quan:
        if not linear:
            args.weight_decay = 5e-2
            args.dropout = 0.3  ##nonlinear quan 0.1 0.3
            args.FILTER = 28  ##nonlinear quan 12 20
            args.gate_init = 20 ##nonlinear quan 10 20
            # args.l1_regularization = 0.000002
            args.l1_regularization = 0.000002
        else:
            args.FILTER = 8
            args.dropout = 0.1
            args.weight_decay = 0.5
            args.gate_init = 20
            args.l1_regularization = 0.000002
    else:
        if not linear:
            args.FILTER = 28   ####20
            args.dropout = 0.3
            args.weight_decay = 5e-2
            args.gate_init = 20
            args.l1_regularization = 0.000002
        else:
            args.weight_decay = 0.1
            args.gate_init = 20
            args.l1_regularization = 0.00002
    
    dataset = TensorDataset(X_tensor, Y_tensor)
    train_size = int(0.8 * len(dataset))
    test_size = len(dataset) - train_size
    train_dataset, test_dataset = random_split(dataset, [train_size, test_size])

    train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=False)
    test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False)

    model = DNN(args).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=args.weight_decay)
    # optimizer = optim.AdamW(model.parameters(), lr=lr, weight_decay=args.weight_decay)
    trainer = TrainerWithCallback(args, model, optimizer, device)
    trained_model, feature_importances = trainer.train(train_loader, test_loader, num_epochs=epoch, epoch_param_current=[epoch], alpha=alpha)
    df = pd.DataFrame(feature_importances.T)
    df.to_csv(f'{base_path}/FI_nn_final_gradient.csv', index=False, header=False)

    # ------------------ SHAP 解释 ------------------
    # 获取 train/test 的全部数据
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

    # 设置模型为只返回预测值（便于计算梯度）
    trained_model.return_y_only = True
    trained_model.eval()
    torch.cuda.empty_cache()

    # combined_tensor = torch.cat((mean_tensor_train, mean_tensor_test, child_tensor_train, child_tensor_test), dim=0)
    # model.return_y_only = True
    # model.eval()
    # n = combined_tensor.shape[0]
    # background_data = combined_tensor[:n//2].to(device)
    # explanation_data = combined_tensor[n//2:].to(device)    
    # explainer = shap.GradientExplainer(model, background_data.to(device))
    # shap_values = explainer.shap_values(explanation_data.to(device))
    # feature_importances = np.abs(shap_values).mean(axis=0)
    # df = pd.DataFrame(feature_importances.squeeze(2).T)
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
    df.to_csv(f'{base_path}/FI_nn_final_shap.csv', index=False, header=False)
init_seed(replicate)
run_one_experiment(sample_size, linear, replicate, 1, quan)