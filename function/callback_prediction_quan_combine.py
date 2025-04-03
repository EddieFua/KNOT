import torch
import numpy as np
import torch.nn as nn
import torch.nn.functional as F
import copy
from sklearn.metrics import r2_score

class TrainerWithCallback:
    def __init__(self, args, model, optimizer, device):
        """Initialize the trainer with arguments, model, optimizer, and device."""
        self.history = []
        self.args = args
        self.model = model.to(device)
        self.optimizer = optimizer
        self.criterion_mse = nn.MSELoss()
        self.device = device
        self.scheduler = torch.optim.lr_scheduler.ExponentialLR(self.optimizer, 0.99)

    def distance_loss(self, mean_embedding, child_embedding):
        """Compute cosine embedding loss to encourage dissimilarity between embeddings."""
        mean_embedding = F.normalize(mean_embedding, dim=1, eps=1e-15)
        child_embedding = F.normalize(child_embedding, dim=1, eps=1e-15)
        mean_embedding = mean_embedding.to(self.device)
        child_embedding = child_embedding.to(mean_embedding.device)
        loss_fn = nn.CosineEmbeddingLoss(margin=0.0)
        labels_inter = -torch.ones(mean_embedding.size(0)).to(mean_embedding.device)
        return loss_fn(mean_embedding, child_embedding, labels_inter)

    def l1_regularizer(self, model):
        """Apply L1 regularization to model weights that require gradients."""
        l1_loss = 0.0
        for name, param in model.named_parameters():
            if "weight" in name and param.requires_grad:
                l1_loss += torch.sum(torch.abs(param))
        return self.args.l1_regularization * l1_loss

    def validate_and_get_fi(self, train_loader, val_loader):
        """Validate model and compute feature importance using gradients."""
        self.model.eval()
        # Validation predictions and R² calculation
        preds_list, labels_list = [], []
        with torch.no_grad():
            for X_batch, Y_batch in val_loader:
                X_batch, Y_batch = X_batch.to(self.device), Y_batch.to(self.device)
                dad_tensor = X_batch[:, :, :, 0].float()
                mom_tensor = X_batch[:, :, :, 1].float()
                child_tensor = X_batch[:, :, :, 2].float()
                mean_tensor = (dad_tensor + mom_tensor) / 2
                combined_tensor = torch.cat((mean_tensor, child_tensor), dim=0)

                self.model.return_y_only = True
                output = self.model(combined_tensor).cpu().numpy()
                preds_list.append(output)

                target_dad_y = Y_batch[:, 0].reshape(-1, 1).cpu().numpy()
                target_mom_y = Y_batch[:, 1].reshape(-1, 1).cpu().numpy()
                target_mean_y = (target_dad_y + target_mom_y) / 2
                target_child_y = Y_batch[:, 2].reshape(-1, 1).cpu().numpy()
                cat_y = np.concatenate((target_mean_y, target_child_y), axis=0)
                labels_list.append(cat_y)

        preds_array = np.concatenate(preds_list)
        labels_array = np.concatenate(labels_list)
        r2 = r2_score(labels_array, preds_array)

        # Feature importance computation using gradients
        train_data = torch.cat([X_b for X_b, _ in train_loader], dim=0)
        val_data = torch.cat([X_b for X_b, _ in val_loader], dim=0)
        all_data = torch.cat([train_data, val_data], dim=0)
        dad_big = all_data[:, :, :, 0].float().to(self.device)
        mom_big = all_data[:, :, :, 1].float().to(self.device)
        child_big = all_data[:, :, :, 2].float().to(self.device).requires_grad_(True)

        mean_big = ((dad_big + mom_big) / 2).requires_grad_(True)
        combined_tensor = torch.cat((mean_big, child_big), dim=0)
        
        self.model.return_y_only = True
        out_batch = self.model(combined_tensor)
        grad_batch = torch.autograd.grad(
            outputs=out_batch,
            inputs=combined_tensor,
            grad_outputs=torch.ones_like(out_batch),
            retain_graph=False,
            create_graph=False
        )[0]
        
        importance = torch.abs(grad_batch).mean(dim=0)
        return r2, importance.cpu().detach().numpy()

    def compute_top5_weighted_fi(self):
        """Calculate weighted average feature importance from top 5 validation scores."""
        self.history.sort(key=lambda x: x[0], reverse=True)
        top5 = self.history[:5]
        scores = np.array([item[0] for item in top5])
        fi_list = [item[1] for item in top5]
        score_min, score_max = scores.min(), scores.max()
        weights = (np.ones_like(scores) / len(scores) if score_max == score_min 
                  else (scores - score_min) / (score_max - score_min))
        final_fi = np.zeros_like(fi_list[0])
        for w, fi in zip(weights, fi_list):
            final_fi += w * fi
        return final_fi

    def train(self, dataloader, test_loader, num_epochs, alpha):
        """Train the model with early stopping and return model and feature importance."""
        best_val_score = -float('inf')
        best_model_state = None
        epochs_no_improve = 0
        patience, min_delta = 80, 1e-4

        for epoch in range(num_epochs):
            self.model.train()
            total_loss = total_regu_loss = total_distance_loss = total_classification_loss = 0.0
            y_mean_true, y_mean_pred, y_child_true, y_child_pred = [], [], [], []

            for X_batch, Y_batch in dataloader:
                X_batch, Y_batch = X_batch.to(self.device), Y_batch.to(self.device)
                target_dad_y = Y_batch[:, 0].reshape(-1, 1).float()
                target_mom_y = Y_batch[:, 1].reshape(-1, 1).float()
                target_child_y = Y_batch[:, 2].reshape(-1, 1).float()
                dad_tensor = X_batch[:, :, :, 0].float()
                mom_tensor = X_batch[:, :, :, 1].float()
                child_tensor = X_batch[:, :, :, 2].float()

                self.optimizer.zero_grad()
                mean_tensor = (dad_tensor + mom_tensor) / 2
                target_mean_y = (target_dad_y + target_mom_y) / 2
                combined_tensor = torch.cat((mean_tensor, child_tensor), dim=0)
                self.model.return_y_only = False
                z_combined, y_combined = self.model(combined_tensor)
                z_mean, y_mean = z_combined[:mean_tensor.size(0)], y_combined[:mean_tensor.size(0)]
                z_child, y_child = z_combined[mean_tensor.size(0):], y_combined[mean_tensor.size(0):]

                ce_mean = self.criterion_mse(y_mean, target_mean_y)
                ce_child = self.criterion_mse(y_child, target_child_y)
                classification_loss = ce_mean + ce_child
                distance_loss = alpha * self.distance_loss(z_mean, z_child)
                regu_loss = self.l1_regularizer(self.model)
                loss = classification_loss + distance_loss + regu_loss

                loss.backward()
                self.optimizer.step()

                total_loss += loss.item()
                total_classification_loss += classification_loss.item()
                total_distance_loss += distance_loss.item()
                total_regu_loss += regu_loss.item()
                y_mean_true.append(target_mean_y.cpu().detach().numpy())
                y_mean_pred.append(y_mean.cpu().detach().numpy())
                y_child_true.append(target_child_y.cpu().detach().numpy())
                y_child_pred.append(y_child.cpu().detach().numpy())

            y_mean_true = np.concatenate(y_mean_true)
            y_mean_pred = np.concatenate(y_mean_pred)
            y_child_true = np.concatenate(y_child_true)
            y_child_pred = np.concatenate(y_child_pred)

            avg_total_loss = total_loss / len(dataloader)
            avg_classification_loss = total_classification_loss / len(dataloader)
            avg_distance_loss = total_distance_loss / len(dataloader)
            avg_regu_loss = total_regu_loss / len(dataloader)

            val_score, epoch_fi = self.validate_and_get_fi(dataloader, test_loader)
            self.history.append((val_score, epoch_fi, epoch))
            r2_mean = r2_score(y_mean_true, y_mean_pred)
            r2_child = r2_score(y_child_true, y_child_pred)
            print(f"Epoch [{epoch+1}/{num_epochs}], "
                  f"Total Loss: {avg_total_loss:.4f}, "
                  f"Classification Loss: {avg_classification_loss:.4f}, "
                  f"Distance Loss: {avg_distance_loss:.4f}, "
                  f"Regu Loss: {avg_regu_loss:.4f}, "
                  f"R² Mean: {r2_mean:.4f}, R² Child: {r2_child:.4f}, "
                  f"Val R²: {val_score:.4f}")

            self.scheduler.step()

            if val_score - best_val_score > min_delta:
                best_val_score = val_score
                best_model_state = copy.deepcopy(self.model.state_dict())
                epochs_no_improve = 0
            else:
                epochs_no_improve += 1
                if epochs_no_improve >= patience:
                    print(f"Early stopping at epoch {epoch+1}")
                    break

        if best_model_state is not None:
            self.model.load_state_dict(best_model_state)

        return self.model, self.compute_top5_weighted_fi()