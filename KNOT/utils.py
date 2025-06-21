class Args:
    """Configuration class for DNN model and training."""
    def __init__(self, num_features, num_knockoffs, lr, num_epochs, device, weight, quan):
        self.device = device
        self.pVal = num_features
        self.Num_knock = num_knockoffs
        self.FILTER = 28
        self.dropout = 0.3
        self.batch_size = 128
        self.num_epochs = num_epochs
        self.l1_regularization = 0.000002
        self.weight_decay = 0.05
        self.lr = lr
        self.latent_dim = 32
        self.weight = weight
        self.quan = quan
