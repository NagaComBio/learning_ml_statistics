import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.model_selection import train_test_split


print("#"*20)
print("Testing Logistic regression")

bc = datasets.load_breast_cancer()
X, y = bc.data, bc.target
print(X.shape)
print(y.shape)

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size= 0.2, random_state = 42)

def accuracy(y, y_pred):
    return np.sum(y == y_pred) / len(y)

from regression import LogisticRegression
reg = LogisticRegression(iteration=6000, learning_rate=0.01)
reg.fit(train_X, train_y)
y_train_pred = reg.predict(train_X)

train_acc = accuracy(train_y, y_train_pred)
print(f"Train acc: {train_acc}")

y_test_pred = reg.predict(test_X)
test_acc = accuracy(test_y, y_test_pred)
print(f"Test acc: {test_acc}")

############################
print("#"*20)
print("Testing linear regression")

X, y = datasets.make_regression(n_samples = 100, n_features = 3, noise = 20, random_state = 42)
train_X, test_X, train_y, test_y = train_test_split(X, y, test_size= 0.2, random_state = 42)

from regression import LinearRegression
reg = LinearRegression(iteration=7000, learning_rate=0.01)
reg.fit(train_X, train_y)
y_train_pred = reg.predict(train_X)

from regression import Utility
local_metrics = Utility()
train_rsqr = local_metrics.r_squared(train_y, y_train_pred)
print(f"Train R2: {train_rsqr}")

from sklearn import metrics
print(f"SKlearn: {metrics.r2_score(train_y, y_train_pred)}")

y_test_pred = reg.predict(test_X)
test_rsqr = local_metrics.r_squared(test_y, y_test_pred)
print(f"Test R2: {test_rsqr}")
print(f"Sklearn: {metrics.r2_score(test_y, y_test_pred)}")