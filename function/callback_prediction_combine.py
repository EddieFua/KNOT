import torch
import numpy as np
import torch.nn as nn
import torch.optim as optim
import torch.nn.functional as F
import copy
from torch.optim.lr_scheduler import StepLR, ReduceLROnPlateau, CosineAnnealingLR
from sklearn.metrics import roc_auc_score
class TrainerWithCallback:
    """Handles training, validation, and feature importance computation for the DNN model."""
    def __init__(self, args, model, optimizer, device):
        self.args = args
        self.model = model.to(device)
        self.optimizer = optimizer
        self.criterion_ce = nn.BCELoss()
        self.device = device
        self.scheduler = optim.lr_scheduler.ExponentialLR(self.optimizer, 0.99)
        self.history = []

    def distance_loss(self, mean_embedding, child_embedding):
        """Computes cosine embedding loss to encourage dissimilarity between embeddings.
        Args:
            mean_embedding: Tensor of shape [batch, latent_dim]
            child_embedding: Tensor of shape [batch, latent_dim]
        Returns:
            Scalar loss value
        """
        mean_embedding = nn.functional.normalize(mean_embedding, dim=1, eps=1e-15)
        child_embedding = nn.functional.normalize(child_embedding, dim=1, eps=1e-15)
        loss_fn = nn.CosineEmbeddingLoss(margin=0.0)
        labels = -torch.ones(mean_embedding.size(0), device=self.device)
        return loss_fn(mean_embedding, child_embedding, labels)

    def l1_regularizer(self, model):
        """Applies L1 regularization to model weights.
        Args:
            model: DNN model instance
        Returns:
            Scalar regularization loss
        """
        l1_loss = 0.0
        for name, param in model.named_parameters():
            if "weight" in name and param.requires_grad:
                l1_loss += torch.sum(torch.abs(param))
        return self.args.l1_regularization * l1_loss

    def validate_and_get_fi(self, train_loader, val_loader):
        """Validates the model and computes gradient-based feature importance.
        Args:
            train_loader: DataLoader for training data
            val_loader: DataLoader for validation data
        Returns:
            val_score: ROC AUC score on validation set
            epoch_feature_importance: Feature importance array
        """
        self.model.eval()
        preds_list, labels_list = [], []
        with torch.no_grad():
            for X_batch, Y_batch in val_loader:
                Y_batch = Y_batch.to(self.device)
                mean_tensor = (X_batch[:, :, :, 0] + X_batch[:, :, :, 1]) / 2
                child_tensor = X_batch[:, :, :, 2]
                combined_tensor = torch.cat((mean_tensor, child_tensor), dim=0).to(self.device)
                self.model.return_y_only = True
                output = self.model(combined_tensor)
                pred = (output >= 0.5).float().cpu().numpy()
                preds_list.append(pred)
                target_y = torch.cat((Y_batch[:, 0:1], Y_batch[:, 2:3]), dim=0).cpu().numpy()
                labels_list.append(target_y)
        preds_array = np.concatenate(preds_list)
        labels_array = np.concatenate(labels_list)
        val_score = roc_auc_score(labels_array, preds_array)

        combined_grad_accum, total_samples = None, 0
        for loader in (train_loader, val_loader):
            for X_batch, _ in loader:
                mean_tensor = (X_batch[:, :, :, 0] + X_batch[:, :, :, 1]) / 2
                child_tensor = X_batch[:, :, :, 2]
                mean_tensor = mean_tensor.to(self.device).requires_grad_(True)
                child_tensor = child_tensor.to(self.device).requires_grad_(True)
                combined_tensor = torch.cat((mean_tensor, child_tensor), dim=0)
                out_batch = self.model(combined_tensor)
                grad_batch = torch.autograd.grad(
                    outputs=out_batch,
                    inputs=combined_tensor,
                    grad_outputs=torch.ones_like(out_batch)
                )[0]
                abs_grad = torch.abs(grad_batch)
                batch_grad_sum = torch.sum(abs_grad, dim=0)
                combined_grad_accum = batch_grad_sum if combined_grad_accum is None else combined_grad_accum + batch_grad_sum
                total_samples += combined_tensor.shape[0]
        importance = combined_grad_accum / total_samples
        return val_score, importance.cpu().detach().numpy()

    def compute_top5_weighted_fi(self):
        """Computes weighted average of feature importances from top 5 models based on validation score.
        Returns:
            Weighted feature importance array
        """
        self.history.sort(key=lambda x: x[0], reverse=True)
        top5 = self.history[:5]
        scores = np.array([item[0] for item in top5])
        fi_list = [item[1] for item in top5]
        weights = (scores - scores.min()) / (scores.max() - scores.min() + 1e-10)
        return np.average(fi_list, weights=weights, axis=0)

    def train(self, dataloader, test_loader, num_epochs, epoch_param_current, alpha):
        """Trains the model with a combination of classification, distance, and regularization losses.
        Args:
            dataloader: Training DataLoader
            test_loader: Validation DataLoader
            num_epochs: Number of training epochs
            epoch_param_current: List with current epoch parameter
            alpha: Weight for distance loss
        Returns:
            model: Trained model
            overall_feature_importances: Final feature importance array
        """
        best_val_score, best_model_state, epochs_no_improve = -float('inf'), None, 0
        patience = 80
        for epoch in range(num_epochs):
            self.model.train()
            total_loss, total_cls_loss, total_dist_loss, total_regu_loss = 0.0, 0.0, 0.0, 0.0
            total_correct_mean, total_samples_mean = 0, 0
            total_correct_child, total_samples_child = 0, 0

            for X_batch, Y_batch in dataloader:
                X_batch, Y_batch = X_batch.to(self.device), Y_batch.to(self.device)
                target_mean_y, target_child_y = Y_batch[:, 0:1], Y_batch[:, 2:3]
                mean_tensor = (X_batch[:, :, :, 0] + X_batch[:, :, :, 1]) / 2
                child_tensor = X_batch[:, :, :, 2].requires_grad_(True)
                combined_tensor = torch.cat((mean_tensor, child_tensor), dim=0)

                self.model.return_y_only = False
                z_combined, y_combined = self.model(combined_tensor)
                z_mean, y_mean = z_combined[:mean_tensor.size(0)], y_combined[:mean_tensor.size(0)]
                z_child, y_child = z_combined[mean_tensor.size(0):], y_combined[mean_tensor.size(0):]

                cls_loss = self.criterion_ce(y_mean, target_mean_y) + self.criterion_ce(y_child, target_child_y)
                dist_loss = alpha * self.distance_loss(z_mean, z_child)
                regu_loss = self.l1_regularizer(self.model)
                loss = cls_loss + dist_loss + regu_loss

                self.optimizer.zero_grad()
                loss.backward()
                self.optimizer.step()

                total_loss += loss.item()
                total_cls_loss += cls_loss.item()
                total_dist_loss += dist_loss.item()
                total_regu_loss += regu_loss

                # Accuracy for mean predictions
                pred_mean = (y_mean >= 0.5).float()
                total_correct_mean += (pred_mean == target_mean_y).sum().item()
                total_samples_mean += target_mean_y.size(0)
                # Accuracy for child predictions
                pred_child = (y_child >= 0.5).float()
                total_correct_child += (pred_child == target_child_y).sum().item()
                total_samples_child += target_child_y.size(0)

            avg_loss = total_loss / len(dataloader)
            avg_cls_loss = total_cls_loss / len(dataloader)
            avg_dist_loss = total_dist_loss / len(dataloader)
            avg_regu_loss = total_regu_loss / len(dataloader)
            acc_mean = total_correct_mean / total_samples_mean
            acc_child = total_correct_child / total_samples_child

            val_score, epoch_fi = self.validate_and_get_fi(dataloader, test_loader)
            self.history.append((val_score, epoch_fi, epoch))
            print(f"Epoch [{epoch+1}/{num_epochs}], Loss: {avg_loss:.4f}, "
                  f"Cls: {avg_cls_loss:.4f}, Dist: {avg_dist_loss:.4f}, Regu: {avg_regu_loss:.4f}, "
                  f"Mean Acc: {acc_mean:.4f}, Child Acc: {acc_child:.4f}, Val Score: {val_score:.4f}")

            self.scheduler.step()
            if val_score > best_val_score + 1e-4:
                best_val_score = val_score
                best_model_state = copy.deepcopy(self.model.state_dict())
                epochs_no_improve = 0
            else:
                epochs_no_improve += 1
                if epochs_no_improve >= patience:
                    print(f"Early stopping at epoch {epoch+1}")
                    break

        if best_model_state:
            self.model.load_state_dict(best_model_state)
        return self.model, self.compute_top5_weighted_fi()
