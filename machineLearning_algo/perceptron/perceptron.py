import numpy as np


class Perceptron():
    def __init__(self, learning_rate=0.01, iteration=1000):
        self.lr = learning_rate
        self.iter = iteration
        self.weights = None
        self.bias = None
        self.activation_func = self._update_rule
    
    def fit(self, X, y):

        n_samples, n_feature = X.shape
        self.weights = np.zeros(n_feature)
        self.bias = 0

        for idx, x in enumerate(X):
            y_linear = np.dot(x, self.weights) + self.bias
            y_pred = self.activation_func(y_linear)
            update = self.lr * (y[idx] - y_pred)

            self.weights += update * x
            self.bias += update
    
    def predict(self, X):
        y_linear = np.dot(X, self.weights) + self.bias
        y_pred = self.activation_func(y_linear)
        return(y_pred)

    def _update_rule(self, x):
        return np.where(x > 0, 1, 0)