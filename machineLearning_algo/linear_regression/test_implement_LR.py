## Test the implementation
## by comparing it with sklearn
from sklearn.linear_model import LinearRegression
import pandas as pd
import numpy as np

## Sample housing data
dat = pd.read_csv("data/kc_house_data.csv")
print(dat.head())

## Regression from sklearn
lr = LinearRegression()
x = dat['sqft_living15'].to_numpy()[1:1000].reshape(-1,1)
y = dat['price'].to_numpy()[1:1000].reshape(-1,1)
x = x / np.mean(x)
y = y / np.mean(y)
print(x.shape)
print(y.shape)
lr.fit(x, y)

print(lr.coef_)
print(lr.intercept_)

####################
# Now from the implementation
from implement_linear_regression import get_differentiate, gradient_descent
#x = list(x)
#y = list(y)
x = x / np.mean(x)
y = y/ np.mean(y)
gradient_descent(x, y, learning_rate=0.0005, n_iter=15000, mini_batch=1)