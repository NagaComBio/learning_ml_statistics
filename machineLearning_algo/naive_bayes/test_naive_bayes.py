import numpy as np
import matplotlib.pyplot as plt
from sklearn import datasets
from sklearn.model_selection import train_test_split

bc = datasets.load_breast_cancer()
X, y = bc.data, bc.target
print(X.shape)
print(y.shape)

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size= 0.2, random_state = 42)

def accuracy(y, y_pred):
    return np.sum(y == y_pred) / len(y)

from naive_bayes import Naive_Bayes
reg = Naive_Bayes()
reg.fit(train_X, train_y)
y_train_pred = reg.predict(train_X)

train_acc = accuracy(train_y, y_train_pred)
print(f"Train acc: {train_acc}")

y_test_pred = reg.predict(test_X)
test_acc = accuracy(test_y, y_test_pred)
print(f"Test acc: {test_acc}")
