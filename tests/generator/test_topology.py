import numpy as np
from topology import tetrahedron

si = np.array([0.0, 0.0, 0.0])

o = tetrahedron(si, 1.62)

for i, p in enumerate(o):
    print(i + 1, p)
