import numpy as np
from sklearn import datasets
from sklearn.model_selection import train_test_split

iris = datasets.load_iris()
X = iris.data
y = iris.target

train_X, test_X, train_y, test_y = train_test_split(X, y, test_size = 0.2)

from knn import KNN
for k in range(3, 10):
    print(f"Accuracy for K: {k}")
    cls = KNN(k)
    cls.fit(train_X, train_y)

    pred_train_labels = cls.predict(train_X)
    acc = np.sum(train_y == pred_train_labels) / len(pred_train_labels)
    print(f"\ttrain acc: {acc}")

    pred_test_labels = cls.predict(test_X)
    acc = np.sum(test_y == pred_test_labels) / len(pred_test_labels)
    print(f"\ttrain acc: {acc}")

