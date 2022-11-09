import numpy as np
import random
import matplotlib.pyplot as plt

def generate_dataset(n, b, e):
    '''
    Generate a random data set between 0, 1
    '''
    x = np.random.rand(n)
    y = x * b + e

    return x, y

def cost_function(x, true_y, pred_b, pred_e=0):
    '''
    Sum of squared errors
    '''
    pred_y = x * pred_b + pred_e
    diff_y = np.square(pred_y - true_y)
    m = x.shape[0]
    cost = 1/(2*m) * np.sum(diff_y)

    return(cost)

def brute_force_statergy(x, y):
    '''
    Brute force staregy to find the parameters with minimum 
    '''
    b, e, z = list(), list(), list()
    for i in np.arange(0, 1, 0.01):
        for k in np.arange(0, 1, 0.01):
            cst = cost_function(x, y, i, k)
            z.append(cst)
            b.append(i)
            e.append(k)
            #print(f"{i}\t{k}\t{cst}")

    # Z = np.array(z).reshape((len(b), len(e)))
    # plt.figure()
    # cm = plt.cm.get_cmap('viridis')
    # plt.contour(b, e, Z)
    # plt.show()

    ## Find the parameters with minimum cost
    min_cost_idx = np.argmin(z)
    print("With brute force")
    print(f"Min cost: {z[min_cost_idx]}\nb: {b[min_cost_idx]}\ne: {e[min_cost_idx]}")

def get_differentiate(x, y, b, e, mini_batch = 1):
    '''
    Compute diffentiate for stochastic & mini batch gradient descent
    '''
    diff_b = 0
    diff_e = 0
    
    if mini_batch >= 1:
        # stochastic & mini batch
        m = mini_batch
        batch_loop = random.sample(range(m), mini_batch)
    else:
        # default gradient descent
        m = x.shape[0]
        batch_loop = range(m)
   
    for i in batch_loop:
        f_x = b * x[i] + e
        diff_b_i = (f_x - y[i]) * x[i]
        diff_e_i = (f_x - y[i])
        diff_b += diff_b_i
        diff_e += diff_e_i
    
    diff_b = diff_b / m
    diff_e = diff_e / m 

    return diff_b, diff_e

def gradient_descent(x, y, n_iter=10000, learning_rate=0.001, mini_batch=1):
    '''
    Gradient descent
    '''
    b = np.random.rand()
    e = np.random.rand()
    hist_cost, hist_b, hist_e = list(), list(), list()
    for i in range(n_iter):
        diff_b, diff_e = get_differentiate(x, y, b, e, mini_batch)
        #print(diff_b)
        #print(diff_e)
        b = b - learning_rate * diff_b
        e = e - learning_rate * diff_e
        cst = cost_function(x, y, b, e)
        hist_cost.append(cst)
        hist_b.append(b)
        hist_e.append(e)
    
        if i % 200 == 0:
            print(f"Iter: {i}\tCost: {cst}\tb: {b}\te: {e}")
            #print(f"Iter: {i}\tCost: {cst}")

        if i > 10:
            if abs(cst - hist_cost[i-2]) < 1e-10:
                print(f"Final: Iter: {i}\tCost: {cst}\tb: {b}\te: {e}")
                break
        
        if np.isnan(cst):
            break


if __name__ == "__main__":
    # Generate the datasets
    x, y = generate_dataset(100, 0.4, 0.3)
    #brute_force_statergy(x, y)

    gradient_descent(x, y, learning_rate=0.005, n_iter=55000, mini_batch=10)