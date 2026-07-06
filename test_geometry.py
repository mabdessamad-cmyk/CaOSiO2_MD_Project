import numpy as np

from geometry import *

box = 10.0

rng = np.random.default_rng(12345)

print("Random position")
print(random_position(box, rng))

print()

p1 = np.array([9.8, 5.0, 5.0])
p2 = np.array([0.2, 5.0, 5.0])

print("Normal difference")
print(p1 - p2)

print()

print("Minimum image")
print(minimum_image(p1 - p2, box))

print()

print("Distance")
print(distance(p1, p2, box))
