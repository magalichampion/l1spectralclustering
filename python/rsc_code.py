from rsc.rsc import RSC
import numpy as np

def run_rsc(A, k, theta=10):
    # theta is the percentage of corrupted edges to remove
    model = RSC(k=int(k), theta=theta)
    labels = model.fit_predict(A)
    return labels
