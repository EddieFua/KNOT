import pyreadr
import numpy as np
import pandas as pd
import torch
from model_combine_AGP import DNN
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
from tqdm import tqdm
import shap
import argparse
import os
import xarray as xr
import json
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset, random_split
from sklearn.model_selection import KFold
from sklearn.metrics import accuracy_score, classification_report
from sklearn.metrics import mean_squared_error, r2_score
import h5py
linear = False
quan = False
if quan:
    from callback_prediction_quan_combine import TrainerWithCallback
else:
    from callback_prediction_combine import TrainerWithCallback
device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
def init_seed(seed):
    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
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
init_seed(123)      
linear_or_not = 'nonlinear'
quan_or_not = ''
alpha = 2  ##1.5
epoch = 30
h5_filepath = '/home/yyang773/KONet-Trio/AGP_data/child_array.h5'
with h5py.File(h5_filepath, 'r') as f:
    child_array = np.array(f['child_array'])
child_array = child_array.transpose(2, 1, 0)
h5_filepath = '/home/yyang773/KONet-Trio/AGP_data/dad_array.h5'
with h5py.File(h5_filepath, 'r') as f:
    dad_array = np.array(f['dad_array'])
dad_array = dad_array.transpose(2, 1, 0)
h5_filepath = '/home/yyang773/KONet-Trio/AGP_data/mom_array.h5'
with h5py.File(h5_filepath, 'r') as f:
    mom_array = np.array(f['mom_array'])
mom_array = mom_array.transpose(2, 1, 0)
################
weight = pd.read_csv('/home/yyang773/KONet-Trio/AGP_data/weight.csv', header=0)
weight = weight.apply(pd.to_numeric, errors='coerce').to_numpy()
weight = torch.tensor(weight, dtype=torch.float32)
weight_tensor = weight

reduced_dad_array = dad_array
reduced_mom_array = mom_array
reduced_child_array = child_array
num_samples, num_features, num_knockoffs = reduced_child_array.shape
print(reduced_dad_array.shape)


lr = 0.0001
num_epochs = 50
merged_array = np.stack((reduced_dad_array, reduced_mom_array, reduced_child_array), axis=-1)
X_tensor = torch.tensor(merged_array, dtype=torch.float32)
dad_y = torch.zeros(num_samples)
mom_y = torch.zeros(num_samples)
child_y = torch.ones(num_samples)
Y_tensor = torch.stack((dad_y, mom_y, child_y), dim=1)
args = Args(num_features, num_knockoffs - 1, lr, num_epochs, device, weight_tensor, quan, weight_cnn=weight_tensor, gate_init=20)
args.weight = weight_tensor.clone()
args.weight_cnn = weight_tensor.clone() 
for k in range(args.Num_knock + 1):
    col_min = args.weight[k].min().item()
    col_max = args.weight[k].max().item()
    args.weight[k] = (args.weight[k] - col_min) / (col_max - col_min)
    
    
args.weight_decay = 0.05
args.dropout = 0.3  ##nonlinear quan 0.1
args.FILTER = 4  ##nonlinear quan 12
args.gate_init = 5 ##nonlinear quan 10
args.l1_regularization = 0.0000002

dataset = TensorDataset(X_tensor, Y_tensor)
train_size = int(0.8 * len(dataset))
test_size = len(dataset) - train_size
train_dataset, test_dataset = random_split(dataset, [train_size, test_size])

train_loader = DataLoader(train_dataset, batch_size=args.batch_size, shuffle=True)
test_loader = DataLoader(test_dataset, batch_size=args.batch_size, shuffle=False)

model = DNN(args).to(device)
optimizer = optim.Adam(model.parameters(), lr=lr, weight_decay=args.weight_decay)
trainer = TrainerWithCallback(args, model, optimizer, device)
trained_model, feature_importances = trainer.train(train_loader, test_loader, num_epochs=epoch, epoch_param_current=[epoch], alpha=alpha)
df = pd.DataFrame(feature_importances.T)
df.to_csv('/home/yyang773/KONet-Trio/AGP_data/FI_gradient.csv', index=False, header=False)


# # ------------------ SHAP 解释 ------------------
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
# model_cpu = trained_model.cpu()
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
if mean_abs_shap.ndim > 2:
    mean_abs_shap = np.squeeze(mean_abs_shap, axis=-1)
df = pd.DataFrame(mean_abs_shap.T)
df.to_csv('/home/yyang773/KONet-Trio/AGP_data/FI_shap.csv', index=False, header=False)