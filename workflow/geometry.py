"""
workflow/geometry.py

Geometry utilities for periodic simulation boxes.
"""

import numpy as np


def minimum_image(dx, L):
    """
    Apply the minimum image convention.

    Parameters
    ----------
    dx : float
        Coordinate difference.
    L : float
        Box length.

    Returns
    -------
    float
    """

    return dx - L * np.rint(dx / L)


def vector(r1, r2, box_lengths):
    """
    Shortest vector from r1 to r2 under periodic boundary conditions.
    """

    dr = np.asarray(r2) - np.asarray(r1)

    for i in range(3):
        dr[i] = minimum_image(dr[i], box_lengths[i])

    return dr


def distance_squared(r1, r2, box_lengths):
    """
    Squared distance under periodic boundary conditions.
    """

    dr = vector(r1, r2, box_lengths)

    return np.dot(dr, dr)


def distance(r1, r2, box_lengths):
    """
    Distance under periodic boundary conditions.
    """

    return np.sqrt(distance_squared(r1, r2, box_lengths))


def angle(v1, v2):
    """
    Angle (degrees) between two vectors.
    """

    v1 = np.asarray(v1)
    v2 = np.asarray(v2)

    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)

    if norm1 == 0 or norm2 == 0:
        return 0.0

    cosine = np.dot(v1, v2) / (norm1 * norm2)

    cosine = np.clip(cosine, -1.0, 1.0)

    return np.degrees(np.arccos(cosine))
