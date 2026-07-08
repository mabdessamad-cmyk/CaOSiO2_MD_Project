"""
geometry.py

Geometric utilities for periodic simulation boxes.
"""

import numpy as np


def random_position(box_length, rng):
    """
    Generate a random position inside a cubic box.


    Parameters
    ----------
    box_length : float
        Box length in Å.

    rng : numpy.random.Generator
        Random number generator.

    Returns
    -------
    numpy.ndarray
        Position (x,y,z)
    """

    return rng.random(3) * box_length


def wrap(position, box_length):
    """
    Wrap a point back into the simulation box.
    """

    return position % box_length


def minimum_image(dr, box_length):
    """
    Apply the minimum-image convention.

     Parameters
    ----------
    dr : ndarray
        position_i - position_j
    """

    return dr - box_length * np.round(dr / box_length)


def distance_squared(pos1, pos2, box_length):
    """
    Squared distance using periodic boundaries.
    """

    dr = minimum_image(pos1 - pos2, box_length)

    return np.dot(dr, dr)


def distance(pos1, pos2, box_length):
    """
    Distance using periodic boundaries.
    """

    return np.sqrt(distance_squared(pos1, pos2, box_length))
