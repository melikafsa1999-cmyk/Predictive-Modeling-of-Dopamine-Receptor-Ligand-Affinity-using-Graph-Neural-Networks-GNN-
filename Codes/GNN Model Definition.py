

"""
GNN Model Definition
"""

import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import NNConv, global_mean_pool

class GNNModel(nn.Module):
    def __init__(self, node_in_feats, edge_in_feats, hidden_dim=64, num_layers=3, output_dim=1):
        super(GNNModel, self).__init__()
        self.num_layers = num_layers
        
        # Edge-conditioned convolution layers
        self.convs = nn.ModuleList()
        self.bns = nn.ModuleList()
        
        for i in range(num_layers):
            in_channels = node_in_feats if i == 0 else hidden_dim
            edge_nn = nn.Sequential(
                nn.Linear(edge_in_feats, hidden_dim * in_channels),
                nn.ReLU(),
                nn.Linear(hidden_dim * in_channels, hidden_dim * in_channels)
            )
            conv = NNConv(
                in_channels, hidden_dim,
                nn=edge_nn,
                aggr='mean'
            )
            self.convs.append(conv)
            self.bns.append(nn.BatchNorm1d(hidden_dim))
        
        # Fully connected layers after pooling
        self.fc1 = nn.Linear(hidden_dim, hidden_dim)
        self.fc2 = nn.Linear(hidden_dim, output_dim)

    def forward(self, x, edge_index, edge_attr, batch):
        for i in range(self.num_layers):
            x = self.convs[i](x, edge_index, edge_attr)
            x = self.bns[i](x)
            x = F.relu(x)
        
        # Global mean pooling
        x = global_mean_pool(x, batch)
        
        # Fully connected layers
        x = F.relu(self.fc1(x))
        x = self.fc2(x)
        return x
