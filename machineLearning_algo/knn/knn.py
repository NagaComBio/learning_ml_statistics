import numpy as np
from collections import Counter

def euc_dist(x, y):
    d = np.sqrt(np.sum((x - y)**2))
    return(d)

class KNN:
    
    def __init__(self, k):
        self.k = k
    
    def fit(self, X, y):
        self.train_x = X
        self.train_y = y
    
    def predict(self, X):
        pred_labels = [self._predict(x) for x in X]
        pred_labels = np.array(pred_labels)
        return(pred_labels)
    
    def _predict(self, x):
        distances = [euc_dist(x, x_train) for x_train in self.train_x]
        distances = np.array(distances)

        indices = np.argsort(distances)[:self.k]
        nearest_labels = self.train_y[indices]

        pred_label = Counter(nearest_labels).most_common(1)
        pred_label = pred_label[0][0]

        return(pred_label)
