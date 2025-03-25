import torch
import torch.nn as nn
import math

# ------------------------------------------------------------------
# Positional Encoding
# ------------------------------------------------------------------
class PositionalEncoding(nn.Module):
    def __init__(self, seq_len, embed_dim):
        """
        Adds a fixed positional encoding to each token in the sequence.
        """
        super(PositionalEncoding, self).__init__()
        position = torch.arange(0, seq_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, embed_dim, 2).float() * (-math.log(10000.0) / embed_dim)
        )
        pe = torch.zeros(seq_len, embed_dim)
        # Even indices -> sin, odd indices -> cos
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        # Register as buffer so it's not considered a learnable parameter
        self.register_buffer('pe', pe.unsqueeze(0))  # shape: [1, seq_len, embed_dim]

    def forward(self, x):
        """
        x shape: [batch, seq_len, embed_dim]
        Returns x plus the positional encoding of the same shape.
        """
        return x + self.pe

# ------------------------------------------------------------------
# Locally Connected 1D
# ------------------------------------------------------------------
class LocallyConnected1D(nn.Module):
    def __init__(self, in_channels, out_channels, input_length, kernel_size, stride=1, bias=True):
        """
        A 1D layer similar to convolution, but without weight sharing.
        Each position has its own filters.
        """
        super(LocallyConnected1D, self).__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.kernel_size = kernel_size
        self.stride = stride
        self.input_length = input_length
        
        # Compute the output length (similar to conv output formula)
        self.output_length = (input_length - kernel_size) // stride + 1
        
        # Weights: [output_length, out_channels, in_channels, kernel_size]
        self.weight = nn.Parameter(
            torch.randn(self.output_length, out_channels, in_channels, kernel_size)
        )
        if bias:
            self.bias = nn.Parameter(torch.randn(out_channels, self.output_length))
        else:
            self.register_parameter('bias', None)

    def forward(self, x):
        """
        x shape: [batch_size, in_channels, input_length]
        Returns shape: [batch_size, out_channels, output_length]
        """
        # Unfold to extract local patches
        x_unfolded = x.unfold(dimension=2, size=self.kernel_size, step=self.stride)
        # x_unfolded shape: [batch_size, in_channels, output_length, kernel_size]
        
        # Rearrange dims for einsum
        x_unfolded = x_unfolded.permute(0, 2, 1, 3)
        # shape: [batch_size, output_length, in_channels, kernel_size]
        
        # Multiply patches by position-specific weights
        outputs = torch.einsum('blik,loik->bol', x_unfolded, self.weight)
        # shape: [batch_size, out_channels, output_length]
        
        # Add bias if applicable
        if self.bias is not None:
            outputs += self.bias.unsqueeze(0)
        
        return outputs

# ------------------------------------------------------------------
# Weighted Input
# ------------------------------------------------------------------
class WeightedInput(nn.Module):
    def __init__(self, weight, gate_init = 20, temperature_init = 1):
        """
        weight: [pVal, Num_knock+1] -> stored transposed as a parameter.
        Multiplying the input x by this weight implements a learnable
        per-feature scaling.
        """
        super(WeightedInput, self).__init__()
        self.register_buffer('weight', weight)
        self.temperature_param = nn.Parameter(torch.tensor(temperature_init, dtype=torch.float32))
        self.gate = nn.Parameter(
            torch.full_like(weight, gate_init),  # 与 weight 同形的门控参数
            requires_grad=True
        )

    def forward(self, x):
        temperature = torch.sigmoid(self.temperature_param) * 10  # Scale temperature to be between 0 and 10
        adaptive_gate = torch.sigmoid(self.gate / temperature)  # Apply temperature to gate
        
        dynamic_weight = (adaptive_gate * self.weight + 
                          (1 - adaptive_gate) * torch.ones_like(self.weight))

        return x * dynamic_weight.T.unsqueeze(0)

# ------------------------------------------------------------------
# Main Model: DNN
# ------------------------------------------------------------------
class DNN(nn.Module):
    def __init__(self, Args):
        """
        Args should include:
          - weight:        initial weights for WeightedInput
          - Num_knock:     number of knockoffs
          - pVal:          sequence length (number of SNPs)
          - FILTER, KERNEL, STRIDE: parameters for cnn2 (if used)
          - latent_dim:    dimension of final hidden layer
          - quan:          if True => regression; else => classification
        """
        super(DNN, self).__init__()
        self.args = Args
        self.return_y_only = False
        self.dropout = nn.Dropout(0)
        # 1) Weighted input layer
        self.weighted_input = WeightedInput(Args.weight, Args.gate_init)

        # 2) First LocallyConnected1D
        self.cnn1 = LocallyConnected1D(
            in_channels=self.args.Num_knock + 1,
            out_channels=self.args.FILTER,
            input_length=self.args.pVal,
            kernel_size=1
        )

        self.embed_dim = self.args.FILTER
        self.proj_bn = nn.BatchNorm1d(self.embed_dim)
        self.pos_encoding = PositionalEncoding(
            seq_len=self.args.pVal, 
            embed_dim=self.embed_dim
        )
        self.mha = nn.MultiheadAttention(
            embed_dim=self.embed_dim,
            num_heads=4,   ###4
            dropout=self.args.dropout,   ###0.1
            batch_first=True
        )
        self.attn_norm = nn.LayerNorm(self.embed_dim)
        self.flatten = nn.Flatten()
        output_size = self.args.pVal * self.embed_dim
        
        # 8) Linear layer from 'output_size' -> latent_dim
        self.local2 = nn.Linear(output_size, Args.latent_dim)
        self.elu = nn.ELU()
        # 1D BatchNorm on shape [batch, latent_dim]
        self.batchnorm = nn.BatchNorm1d(Args.latent_dim, eps=1e-15)

        # 9) Final classifier layers
        self.dense2 = nn.Linear(Args.latent_dim, Args.latent_dim // 2)
        self.dense3 = nn.Linear(Args.latent_dim // 2, 1)
        self.sigmoid = nn.Sigmoid()

        # ----------------------------------------------------------------------
        # Initialization
        # ----------------------------------------------------------------------
        nn.init.constant_(self.cnn1.weight, 0.1)
        nn.init.kaiming_normal_(self.local2.weight)
        nn.init.kaiming_normal_(self.dense2.weight)
        nn.init.kaiming_normal_(self.dense3.weight)
        

    def encoder(self, x):
        """
        Returns a latent vector z of shape: [batch, latent_dim].
        """
        # x shape: [batch, pVal, Num_knock+1]

        # 1) Weighted input
        # x = self.weighted_input(x)                # [batch, pVal, Num_knock+1]
        x = x.permute(0, 2, 1)                    # [batch, (Num_knock+1), pVal]

        # 2) LocallyConnected1D => [batch, 1, pVal]
        # 在DNN类的encoder部分：
        x = self.cnn1(x)          # 先卷积
        x = self.proj_bn(x)  # BatchNorm before activation
        x = self.elu(x)      # Activation
        
        # x = self.flatten(x)                       # [batch, pVal]
        # x = self.elu(self.proj(x))                          # [batch, pVal, embed_dim]

        # If needed, you can use cnn2:
        # x = self.elu(self.cnn2(x))                # [batch, FILTER, output_length]
        # x = self.flatten(x)                       # [batch, FILTER * output_length]
    
    
        # 3) Project to embed_dim
        x = x.permute(0, 2, 1)                    # [batch, pVal, 1]
        x = self.pos_encoding(x)                  # [batch, pVal, embed_dim]
        attn_out, _ = self.mha(x, x, x)
        x = x + attn_out
        x = self.attn_norm(x)          # [batch, pVal, embed_dim]
        x = self.flatten(x)
        
        # 7) Map to latent_dim
        x = self.elu(self.local2(x))             # [batch, latent_dim]                     # [batch, latent_dim]
        z = self.batchnorm(x)
        return z

    def classifier(self, z):
        """
        Takes z of shape [batch, latent_dim] and outputs y of shape [batch, 1].
        """
        z = self.elu(self.dense2(z))
        z = self.dropout(z)
        if self.args.quan:
            y = self.dense3(z)
        else:
            y = self.sigmoid(self.dense3(z))
        return y

    def forward(self, x):
        """
        Forward pass: returns (z, y) unless self.return_y_only = True,
        in which case it just returns y.
        """
        z = self.encoder(x)
        y = self.classifier(z)
        if self.return_y_only:
            return y
        else:
            return z, y