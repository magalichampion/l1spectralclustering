from rsc.clustering import RSC
import numpy as np

if not hasattr(np, "int"):
    np.int = int
if not hasattr(np, "float"):
    np.float = float
    
def run_rsc(A, k, theta=10):
    # theta is the percentage of corrupted edges to remove
    A_np = np.array(A, dtype=np.float64)
    
    # Initialize the model
    model = RSC(k=int(k), theta=float(theta))
    
    # fit_predict often fails on the C-level if A is not contiguous
    labels = model.fit_predict(np.ascontiguousarray(A_np))
    
    return labels
