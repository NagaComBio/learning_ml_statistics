import numpy as np

class Naive_Bayes():

    def __init__(self):
        pass

    def fit(self, X, y):
        self.classes = np.unique(y)
        self.class_len = len(self.classes)
        self.mean = np.zeros((self.class_len, X.shape[1]))
        self.var = np.zeros((self.class_len, X.shape[1]))
        self.prior_y = np.zeros(self.class_len)

        for idx, c in enumerate(np.unique(y)):
            X_c = X[y == c,:]
            self.mean[idx,:] = np.mean(X_c)
            self.var[idx,:] = np.var(X_c)
            self.prior_y[idx] = X_c.shape[0]/self.class_len

    def predict(self, X):
        y_pred = [self._predict(x) for x in X]
        return np.array(y_pred)

    def _predict(self, x):
        posteriors = []
        for i in range(self.mean.shape[0]):
            mean = self.mean[i,:]
            var = self.var[i,:]
            posterior = np.sum(np.log(self._pdf(x, mean, var)))
            posterior = posterior + np.log(self.prior_y[i])
            posteriors.append(posterior)
        
        return self.classes[np.argmax(posteriors)]

    def _pdf(self, x, mean, var):
        numerator = np.exp(-(x-mean)**2 / (2 * var))
        denomintor = np.sqrt(2 * np.pi * var)

        return numerator / denomintor