import numpy as np

class LinearRegression():

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
            y_pred = np.dot(X, self.weights) + self.bias
            dw = (1/ n_samples) * np.dot(X.T, (y_pred - y))
            db = (1/ n_samples) * np.sum(y_pred - y)
        
            self.weights -= self.lr * dw
            self.bias  -= self.lr * db
    
    def predict(self, X):
        y_pred = np.dot(X, self.weights) + self.bias
        return(y_pred)
    
