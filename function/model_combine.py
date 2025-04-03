import torch
import torch.nn as nn
import torch.optim as optim
import numpy as np
import pandas as pd
import math
from torch.utils.data import DataLoader, TensorDataset, random_split
from sklearn.metrics import roc_auc_score
import h5py
from tqdm import tqdm
import shap
import copy

# ------------------------------------------------------------------
# Positional Encoding
# ------------------------------------------------------------------
class PositionalEncoding(nn.Module):
    """Adds fixed positional encodings to input embeddings to capture sequence position."""
    def __init__(self, seq_len, embed_dim):
        super(PositionalEncoding, self).__init__()
        position = torch.arange(0, seq_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, embed_dim, 2).float() * (-math.log(10000.0) / embed_dim)
        )
        pe = torch.zeros(seq_len, embed_dim)
        pe[:, 0::2] = torch.sin(position * div_term)  # Sine for even indices
        pe[:, 1::2] = torch.cos(position * div_term)  # Cosine for odd indices
        self.register_buffer('pe', pe.unsqueeze(0))  # Shape: [1, seq_len, embed_dim]

    def forward(self, x):
        """Adds positional encoding to the input tensor.
        Args:
            x: Tensor of shape [batch, seq_len, embed_dim]
        Returns:
            Tensor with positional encoding added, same shape as input.
        """
        return x + self.pe

# ------------------------------------------------------------------
# Locally Connected 1D Layer
# ------------------------------------------------------------------
class LocallyConnected1D(nn.Module):
    """Implements a 1D locally connected layer without weight sharing, where each position has unique filters."""
    def __init__(self, in_channels, out_channels, input_length, kernel_size, stride=1, bias=True):
        super(LocallyConnected1D, self).__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.kernel_size = kernel_size
        self.stride = stride
        self.input_length = input_length
        self.output_length = (input_length - kernel_size) // stride + 1
        self.weight = nn.Parameter(
            torch.randn(self.output_length, out_channels, in_channels, kernel_size)
        )
        if bias:
            self.bias = nn.Parameter(torch.randn(out_channels, self.output_length))
        else:
            self.register_parameter('bias', None)

    def forward(self, x):
        """Performs the forward pass of the locally connected layer.
        Args:
            x: Tensor of shape [batch, in_channels, input_length]
        Returns:
            Tensor of shape [batch, out_channels, output_length]
        """
        x_unfolded = x.unfold(dimension=2, size=self.kernel_size, step=self.stride)
        x_unfolded = x_unfolded.permute(0, 2, 1, 3)  # [batch, output_length, in_channels, kernel_size]
        outputs = torch.einsum('blik,loik->bol', x_unfolded, self.weight)
        if self.bias is not None:
            outputs += self.bias.unsqueeze(0)
        return outputs

# ------------------------------------------------------------------
# Deep Neural Network (DNN) Model
# ------------------------------------------------------------------
class DNN(nn.Module):
    """Deep Neural Network for genomic data analysis, combining locally connected layers, attention, and classification."""
    def __init__(self, args):
        super(DNN, self).__init__()
        self.args = args
        self.return_y_only = False
        self.dropout = nn.Dropout(args.dropout)

        # Locally connected layer for feature extraction
        self.cnn1 = LocallyConnected1D(
            in_channels=args.Num_knock + 1,
            out_channels=args.FILTER,
            input_length=args.pVal,
            kernel_size=1
        )
        self.embed_dim = args.FILTER
        self.proj_bn = nn.BatchNorm1d(self.embed_dim)
        self.pos_encoding = PositionalEncoding(seq_len=args.pVal, embed_dim=self.embed_dim)
        self.mha = nn.MultiheadAttention(
            embed_dim=self.embed_dim,
            num_heads=4,
            dropout=args.dropout,
            batch_first=True
        )
        self.attn_norm = nn.LayerNorm(self.embed_dim)
        self.flatten = nn.Flatten()
        output_size = args.pVal * self.embed_dim

        # Latent space projection
        self.local2 = nn.Linear(output_size, args.latent_dim)
        self.elu = nn.ELU()
        self.batchnorm = nn.BatchNorm1d(args.latent_dim, eps=1e-15)

        # Classifier layers
        self.dense2 = nn.Linear(args.latent_dim, args.latent_dim // 2)
        self.dense3 = nn.Linear(args.latent_dim // 2, 1)
        self.sigmoid = nn.Sigmoid()

        # Weight initialization
        nn.init.constant_(self.cnn1.weight, 0.1)
        nn.init.kaiming_normal_(self.local2.weight)
        nn.init.kaiming_normal_(self.dense2.weight)
        nn.init.kaiming_normal_(self.dense3.weight)

    def encoder(self, x):
        """Encodes input into a latent representation.
        Args:
            x: Tensor of shape [batch, pVal, Num_knock+1]
        Returns:
            z: Latent tensor of shape [batch, latent_dim]
        """
        x = x.permute(0, 2, 1)  # [batch, Num_knock+1, pVal]
        x = self.cnn1(x)         # [batch, FILTER, pVal]
        x = self.proj_bn(x)      # Batch normalization
        x = self.elu(x)          # Activation
        x = x.permute(0, 2, 1)   # [batch, pVal, FILTER]
        x = self.pos_encoding(x) # Add positional encoding
        attn_out, _ = self.mha(x, x, x)  # Multi-head attention
        x = x + attn_out         # Residual connection
        x = self.attn_norm(x)    # Layer normalization
        x = self.flatten(x)      # [batch, pVal * FILTER]
        x = self.elu(self.local2(x))  # [batch, latent_dim]
        z = self.batchnorm(x)    # Batch normalization
        return z

    def classifier(self, z):
        """Classifies the latent representation.
        Args:
            z: Tensor of shape [batch, latent_dim]
        Returns:
            y: Output tensor of shape [batch, 1]
        """
        z = self.elu(self.dense2(z))
        z = self.dropout(z)
        y = self.dense3(z)
        if not self.args.quan:
            y = self.sigmoid(y)  # Apply sigmoid for classification
        return y

    def forward(self, x):
        """Forward pass of the DNN.
        Args:
            x: Input tensor of shape [batch, pVal, Num_knock+1]
        Returns:
            If return_y_only is True: y [batch, 1]
            Else: (z, y) where z is [batch, latent_dim], y is [batch, 1]
        """
        z = self.encoder(x)
        y = self.classifier(z)
        return y if self.return_y_only else (z, y)