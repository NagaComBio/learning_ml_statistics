import numpy as np

class LogisticRegression():

    def __init__(self, learning_rate=0.001, iteration=50):
        self.weights = None
        self.bias = None
        self.lr = learning_rate
        self.iter = iteration
    
    def fit(self, X, y):
        n_samples, n_features = X.shape

        self.weights = np.random.sample(n_features)
        self.bias = np.random.sample(1)

        for _ in range(self.iter):
            linear_reg = np.dot(X, self.weights) + self.bias
            y_pred = self._sigmoid(linear_reg)
            dw = (1/ n_samples) * np.dot(X.T, (y_pred - y))
            db = (1/ n_samples) * np.sum(y_pred - y)
        
            self.weights -= self.lr * dw
            self.bias  -= self.lr * db
    
    def predict(self, X):
        linear_reg = np.dot(X, self.weights) + self.bias
        y_pred = self._sigmoid(linear_reg)
        y_class = [1 if x > 0.5 else 0 for x in y_pred]
        return np.array(y_class)
    
    def _sigmoid(self, x):
        return 1 / ( 1+ np.exp(-x))
    
