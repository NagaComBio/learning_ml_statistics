import numpy as np

class BaseRegression():

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
            y_pred = self._approximation(X, self.weights, self.bias)
            dw = (1/ n_samples) * np.dot(X.T, (y_pred - y))
            db = (1/ n_samples) * np.sum(y_pred - y)
        
            self.weights -= self.lr * dw
            self.bias  -= self.lr * db
    
    def predict(self, X):
        return self._predict(X, self.weights, self.bias)

    def _approximation(self, X, w, b):
        raise NotImplementedError()
    
    def _predict(self, X, w, b):
        raise NotImplementedError()

class LinearRegression(BaseRegression):
    
    def _approximation(self, X, w, b):
        return np.dot(X, w) + b
    
    def _predict(self, X, w, b):
        return np.dot(X, w) + b

class LogisticRegression(BaseRegression):

    def _approximation(self, X, w, b):
        linear_reg = np.dot(X, w) + b
        y_pred = self._sigmoid(linear_reg)
        return(y_pred)

    def _predict(self, X, w, b):
        linear_reg = np.dot(X, w) + b
        y_pred = self._sigmoid(linear_reg)
        y_class = [1 if x > 0.5 else 0 for x in y_pred]
        return np.array(y_class)
    
    def _sigmoid(self, x):
        return 1 / ( 1+ np.exp(-x))
    
class Utility():
    def __init__(self):
        pass

    def r_squared(self, y, y_pred):
        y_mean = np.mean(y)
        r1 = np.sum((y - y_pred)**2)
        r2 = np.sum((y - y_mean)**2)
        r_sqr = 1 - (r1/r2)
    
        return r_sqr