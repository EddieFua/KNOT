import torch
import torch.nn as nn
import torch.nn.functional as F
import numpy as np
from sklearn.metrics import r2_score
import copy

class TrainerWithCallback:
    """Trainer for quantitative genomic data analysis."""
    def __init__(self, args, model, optimizer, device):
        self.args = args
        self.model = model.to(device)
        self.optimizer = optimizer
        self.criterion_mse = nn.MSELoss()
        self.device = device
        self.scheduler = torch.optim.lr_scheduler.ExponentialLR(self.optimizer, 0.99)
        self.history = []

    def distance_loss(self, mean_embedding, child_embedding):
        """Compute cosine embedding loss for dissimilarity."""
        mean_embedding = F.normalize(mean_embedding, dim=1, eps=1e-15)
        child_embedding = F.normalize(child_embedding, dim=1, eps=1e-15)
        loss_fn = nn.CosineEmbeddingLoss(margin=0.0)
        labels = -torch.ones(mean_embedding.size(0)).to(self.device)
        return loss_fn(mean_embedding, child_embedding, labels)

    def l1_regularizer(self, model):
        """Apply L1 regularization to model weights."""
        l1_loss = sum(torch.sum(torch.abs(param)) for name, param in model.named_parameters()
                      if "weight" in name and param.requires_grad)
        return self.args.l1_regularization * l1_loss

    def validate_and_get_fi(self, train_loader, val_loader):
        """Validate model and compute feature importance."""
        self.model.eval()
        preds_list, labels_list = [], []
        with torch.no_grad():
            for X_batch, Y_batch in val_loader:
                X_batch, Y_batch = X_batch.to(self.device), Y_batch.to(self.device)
                mean_tensor = (X_batch[:, :, :, 0] + X_batch[:, :, :, 1]) / 2
                child_tensor = X_batch[:, :, :, 2]
                combined_tensor = torch.cat((mean_tensor, child_tensor), dim=0)
                self.model.return_y_only = True
                output = self.model(combined_tensor).cpu().numpy()
                preds_list.append(output)
                target_mean_y = ((Y_batch[:, 0] + Y_batch[:, 1]) / 2).reshape(-1, 1).cpu().numpy()
                target_child_y = Y_batch[:, 2].reshape(-1, 1).cpu().numpy()
                labels_list.append(np.concatenate((target_mean_y, target_child_y), axis=0))
        preds_array = np.concatenate(preds_list)
        labels_array = np.concatenate(labels_list)
        r2 = r2_score(labels_array, preds_array)

        all_data = torch.cat([X_b for X_b, _ in train_loader] + [X_b for X_b, _ in val_loader], dim=0)
        mean_big = ((all_data[:, :, :, 0] + all_data[:, :, :, 1]) / 2).to(self.device).requires_grad_(True)
        child_big = all_data[:, :, :, 2].to(self.device).requires_grad_(True)
        combined_tensor = torch.cat((mean_big, child_big), dim=0)
        self.model.return_y_only = True
        out_batch = self.model(combined_tensor)
        grad_batch = torch.autograd.grad(out_batch, combined_tensor, torch.ones_like(out_batch))[0]
        importance = torch.abs(grad_batch).mean(dim=0)
        return r2, importance.cpu().numpy()

    def compute_top5_weighted_fi(self):
        """Calculate weighted average feature importance from top 5 models."""
        self.history.sort(key=lambda x: x[0], reverse=True)
        top5 = self.history[:5]
        scores = np.array([item[0] for item in top5])
        fi_list = [item[1] for item in top5]
        weights = (scores - scores.min()) / (scores.max() - scores.min() + 1e-10)
        return np.average(fi_list, weights=weights, axis=0)

    def train(self, train_loader, test_loader, num_epochs, alpha):
        """Train the model with early stopping."""
        best_val_score, best_model_state, epochs_no_improve = -float('inf'), None, 0
        patience = 5

        for epoch in range(num_epochs):
            self.model.train()
            total_loss = 0.0
            y_mean_true, y_mean_pred, y_child_true, y_child_pred = [], [], [], []

            for X_batch, Y_batch in train_loader:
                X_batch, Y_batch = X_batch.to(self.device), Y_batch.to(self.device)
                mean_tensor = (X_batch[:, :, :, 0] + X_batch[:, :, :, 1]) / 2
                child_tensor = X_batch[:, :, :, 2]
                combined_tensor = torch.cat((mean_tensor, child_tensor), dim=0)
                target_mean_y = ((Y_batch[:, 0] + Y_batch[:, 1]) / 2).reshape(-1, 1)
                target_child_y = Y_batch[:, 2].reshape(-1, 1)

                self.optimizer.zero_grad()
                self.model.return_y_only = False
                z_combined, y_combined = self.model(combined_tensor)
                z_mean, y_mean = z_combined[:mean_tensor.size(0)], y_combined[:mean_tensor.size(0)]
                z_child, y_child = z_combined[mean_tensor.size(0):], y_combined[mean_tensor.size(0):]

                classification_loss = self.criterion_mse(y_mean, target_mean_y) + self.criterion_mse(y_child, target_child_y)
                distance_loss = alpha * self.distance_loss(z_mean, z_child)
                regu_loss = self.l1_regularizer(self.model)
                loss = classification_loss + distance_loss + regu_loss

                loss.backward()
                self.optimizer.step()

                total_loss += loss.item()
                y_mean_true.append(target_mean_y.cpu().numpy())
                y_mean_pred.append(y_mean.cpu().numpy())
                y_child_true.append(target_child_y.cpu().numpy())
                y_child_pred.append(y_child.cpu().numpy())

            val_score, epoch_fi = self.validate_and_get_fi(train_loader, test_loader)
            self.history.append((val_score, epoch_fi, epoch))
            print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {total_loss/len(train_loader):.4f}, Val R²: {val_score:.4f}")
            self.scheduler.step()

            if val_score > best_val_score + 1e-4:
                best_val_score = val_score
                best_model_state = copy.deepcopy(self.model.state_dict())
                epochs_no_improve = 0
            elif epochs_no_improve >= patience:
                print(f"Early stopping at epoch {epoch+1}")
                break
            else:
                epochs_no_improve += 1

        if best_model_state:
            self.model.load_state_dict(best_model_state)
        return self.model, self.compute_top5_weighted_fi()
