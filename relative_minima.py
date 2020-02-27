import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.modeling import models
from scipy.signal import argrelmin


def find_minima(z: np.ndarray, mode: str = 'wrap'):
    """
    Finds local minima in a N-dimensional array.

    Parameters
    ----------
    z : np.ndarray
        Array in which to find the local minima.
    mode : str
        Sets the behaviour for pixels at the borders. This parameter
        is passed directly to scipy.signal.argrelmin, and can be
        either *wrap* or *clip*.

    Returns
    -------
    minima : list
        A list of tuples, each specifying the coordinates of a local
        mininum in the array.

    Notes
    -----
    This function is not capable of identifying flat minima! It
    strictly requires that the values at both sides of a minimum
    be higher than the minimum.

    See also
    --------
    scipy.signal.argrelmin, scipy.signal.extrema, scipy.signal.find_peaks
    """

    def flat_warning():
        warnings.warn('Array has neighbouring duplicate values, possibly related to a flat minimum.', RuntimeWarning,
                      stacklevel=2)

    if z.ndim > 1:
        axis = []
        minima = []
        for dimension in range(z.ndim):
            if np.any(np.diff(z, axis=dimension) == 0):
                flat_warning()
            c = argrelmin(z, axis=dimension, mode=mode)
            d = [tuple([c[j][i] for j in range(z.ndim)]) for i in range(len(c[0]))]
            axis.append(d)

        for c in axis[0]:
            if np.all([c in axis[i] for i in range(1, z.ndim)]):
                minima.append(c)
    else:
        if np.any(np.diff(z) == 0):
            flat_warning()
        minima = argrelmin(z, mode=mode)[0]

    return minima


def test_3d():
    k = np.linspace(-10, 10, 20)
    x = np.array(np.meshgrid(k, k, k))

    def gauss3d(radius: np.ndarray, sigma: float):
        return np.exp(-radius ** 2 / (2 * (sigma ** 2)))

    image = np.zeros_like(x[0])
    for i in range(3):
        x0 = np.random.choice(k, size=3)
        print('x0 =', x0)
        r = np.sqrt(np.sum([np.square(x[_] - x0[_]) for _ in range(3)], axis=0))
        image += -gauss3d(r, 2.0)

    minima_coordinates = find_minima(image)

    print('minima_coordinates =', minima_coordinates)

    return image, minima_coordinates


def test_2d(plot: bool = False):
    k = np.linspace(-10, 10, 20)
    y, x = np.meshgrid(k, k)

    image = np.zeros_like(x)
    for i in range(3):
        x0 = np.random.choice(k, size=2)
        print(x0)
        image += models.Gaussian2D(amplitude=-1, x_mean=x0[0], y_mean=x0[1])(x, y)

    minima_coordinates = find_minima(image)

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(image, origin='lower')
        for coordinate in minima_coordinates:
            ax.scatter(*coordinate)
        ax.plot(minima_coordinates)
        plt.show()

    print(minima_coordinates)

    return image, minima_coordinates


def test_1d(plot: bool = False):
    x = np.linspace(0, 10, 20)
    image = np.zeros_like(x)

    for i in range(3):
        image += models.Gaussian1D(-1, np.random.choice(x))(x)

    minima_coordinates = find_minima(image)

    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(x, image)
        for coordinate in minima_coordinates:
            ax.axvline(x[coordinate])
        plt.show()

    print(minima_coordinates)

    return image, minima_coordinates


if __name__ == '__main__':
    image, min = test_2d(plot=True)
