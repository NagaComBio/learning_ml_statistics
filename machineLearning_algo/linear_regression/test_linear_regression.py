import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.model_selection import train_test_split

X, y = datasets.make_regression(n_samples = 100, n_features = 3, noise = 20, random_state = 42)
print(X.shape)
print(y.shape)
print(X[0].shape)
plt.figure()
plt.scatter(X[:,0], y)
plt.show()

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size= 0.2, random_state = 42)

from linear_regression import LinearRegression
reg = LinearRegression(iteration=7000, learning_rate=0.01)
reg.fit(train_X, train_y)
y_train_pred = reg.predict(train_X)

train_error = np.sum((train_y - y_train_pred)**2)
print(f"Train error: {train_error}")

y_test_pred = reg.predict(test_X)
test_error = np.sum((test_y - y_test_pred)**2)
print(f"Test error: {test_error}")

plt.figure()
plt.scatter(train_X[:,0], train_y)
plt.axline(xy1 = (0, reg.bias), slope = reg.weights[0])
plt.show()



