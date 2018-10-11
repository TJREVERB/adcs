import numpy as np
def crossm(u):
    np = [0 -u(3) u(2);  u(3) 0 -u(1); -u(2) u(1) 0]
    return np

matrix = np.array([1,2];[2,3])
print(crossm(matrix))
