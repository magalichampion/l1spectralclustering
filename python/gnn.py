import torch
import numpy as np
from torch_geometric.datasets import Planetoid
from torch_geometric.nn import VGAE, GCNConv
from torch_geometric.utils import from_scipy_sparse_matrix
from scipy.sparse import csr_matrix
from sklearn.cluster import KMeans

class Encoder(torch.nn.Module):
    def __init__(self, in_channels, out_channels):
        super().__init__()
        self.conv1 = GCNConv(in_channels, 2 * out_channels)
        self.conv_mu = GCNConv(2 * out_channels, out_channels)
        self.conv_logstd = GCNConv(2 * out_channels, out_channels)

    def forward(self, x, edge_index):
        x = self.conv1(x, edge_index).relu()
        return self.conv_mu(x, edge_index), self.conv_logstd(x, edge_index)

def run_vgae(adj_matrix, num_clusters):
    # Convert R matrix to Torch Geometric format
    adj_sparse = csr_matrix(adj_matrix)
    edge_index, _ = from_scipy_sparse_matrix(adj_sparse)
    
    # Use Identity matrix as initial features if none exist
    num_nodes = adj_matrix.shape[0]
    x = torch.eye(num_nodes) 
    
    model = VGAE(Encoder(num_nodes, 16))
    optimizer = torch.optim.Adam(model.parameters(), lr=0.01)

    # Simple training loop
    model.train()
    for epoch in range(100):
        optimizer.zero_grad()
        z = model.encode(x, edge_index)
        loss = model.recon_loss(z, edge_index) + (1 / num_nodes) * model.kl_loss()
        loss.backward()
        optimizer.step()

    # Cluster the latent embeddings (z)
    model.eval()
    with torch.no_grad():
        z = model.encode(x, edge_index)
        kmeans = KMeans(n_clusters=num_clusters, n_init=10).fit(z.cpu().numpy())
    
    return kmeans.labels_ + 1 # +1 to match R's 1-based indexing
