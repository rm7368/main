import pandas as pd
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader

class EnterocyteDataset(Dataset):
  """Dataset for enterocyte domain classification."""
  def __init__(self, expr_file, labels_file):
    """
    Args:
    expr_file: Path to CSV with expression (cellxgene)
    labels_file: Path to CSV with domain labels
    """
    self.expression = pd.read_csv(expr_file).values.astype(np.float32)
    labels_df = pd.read_csv(labels_file)
  
    domain_map = {
    'Domain_A': 0,
    'Domain_B': 1,
    'Domain_C': 2,
    'Domain_D': 3,
    'Domain_E': 4
    }
    self.domain_labels = np.array([domain_map[d] for d in labels_df['domain']])
  
    self.n_cells, self.n_genes = self.expression.shape
    self.n_domains = 5
  
    print(f"Loaded {self.n_cells} cells x {self.n_genes} genes")
    print(f"Domain distribution: {np.bincount(self.domain_labels)}")
    
  def __len__(self):
    return self.n_cells
  
  def __getitem(self,idx):
    """
    Returns:
      dict with 'expression' and 'domain' tensors
    """
    return {
      'expression': torch.FloatTensor(self.expression[idx]),
      'domain': torch.LongTensor([self.domain_labels[idx]])[0]
    }
    
  def get_dataloaders(train_expr, train_labels, test_expr, test_labels, 
                      batch_size=128, num_workers=4):
    
      
