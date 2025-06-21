import torch
import torch.nn as nn
import math

class PositionalEncoding(nn.Module):
    """Adds fixed positional encodings to input embeddings to capture sequence position."""
    def __init__(self, seq_len, embed_dim):
        super().__init__()
        position = torch.arange(seq_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(torch.arange(0, embed_dim, 2).float() * (-math.log(10000.0) / embed_dim))
        pe = torch.zeros(seq_len, embed_dim)
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer('pe', pe.unsqueeze(0))

    def forward(self, x):
        """Adds positional encoding to the input tensor.
        
        Args:
            x: Tensor of shape [batch, seq_len, embed_dim]
        Returns:
            Tensor with positional encoding added, same shape as input.
        """
        return x + self.pe

class LocallyConnected1D(nn.Module):
    """Implements a 1D locally connected layer without weight sharing."""
    def __init__(self, in_channels, out_channels, input_length, kernel_size, stride=1, bias=True):
        super().__init__()
        self.output_length = (input_length - kernel_size) // stride + 1
        self.weight = nn.Parameter(torch.randn(self.output_length, out_channels, in_channels, kernel_size))
        self.bias = nn.Parameter(torch.randn(out_channels, self.output_length)) if bias else None

    def forward(self, x):
        """Performs the forward pass of the locally connected layer.
        
        Args:
            x: Tensor of shape [batch, in_channels, input_length]
        Returns:
            Tensor of shape [batch, out_channels, output_length]
        """
        x_unfolded = x.unfold(dimension=2, size=self.kernel_size, step=self.stride).permute(0, 2, 1, 3)
        outputs = torch.einsum('blik,loik->bol', x_unfolded, self.weight)
        if self.bias is not None:
            outputs += self.bias.unsqueeze(0)
        return outputs

class DNN(nn.Module):
    """Deep Neural Network for genomic data analysis with locally connected layers and attention."""
    def __init__(self, args):
        super().__init__()
        self.args = args
        self.return_y_only = False
        self.dropout = nn.Dropout(args.dropout)

        self.cnn1 = LocallyConnected1D(
            in_channels=args.Num_knock + 1,
            out_channels=args.FILTER,
            input_length=args.pVal,
            kernel_size=1
        )
        self.embed_dim = args.FILTER
        self.proj_bn = nn.BatchNorm1d(self.embed_dim)
        self.pos_encoding = PositionalEncoding(seq_len=args.pVal, embed_dim=self.embed_dim)
        self.mha = nn.MultiheadAttention(embed_dim=self.embed_dim, num_heads=4, dropout=args.dropout, batch_first=True)
        self.attn_norm = nn.LayerNorm(self.embed_dim)
        self.flatten = nn.Flatten()
        output_size = args.pVal * self.embed_dim

        self.local2 = nn.Linear(output_size, args.latent_dim)
        self.elu = nn.ELU()
        self.batchnorm = nn.BatchNorm1d(args.latent_dim, eps=1e-15)
        self.dense2 = nn.Linear(args.latent_dim, args.latent_dim // 2)
        self.dense3 = nn.Linear(args.latent_dim // 2, 1)
        self.sigmoid = nn.Sigmoid()

        nn.init.constant_(self.cnn1.weight, 0.1)
        nn.init.kaiming_normal_(self.local2.weight)
        nn.init.kaiming_normal_(self.dense2.weight)
        nn.init.kaiming_normal_(self.dense3.weight)

    def encoder(self, x):
        """Encodes input into a latent representation."""
        x = x.permute(0, 2, 1)
        x = self.cnn1(x)
        x = self.proj_bn(x)
        x = self.elu(x)
        x = x.permute(0, 2, 1)
        x = self.pos_encoding(x)
        attn_out, _ = self.mha(x, x, x)
        x = x + attn_out
        x = self.attn_norm(x)
        x = self.flatten(x)
        x = self.elu(self.local2(x))
        return self.batchnorm(x)

    def classifier(self, z):
        """Classifies the latent representation."""
        z = self.elu(self.dense2(z))
        z = self.dropout(z)
        y = self.dense3(z)
        if not self.args.quan:
            y = self.sigmoid(y)
        return y

    def forward(self, x):
        """Forward pass of the DNN."""
        z = self.encoder(x)
        y = self.classifier(z)
        return y if self.return_y_only else (z, y)
