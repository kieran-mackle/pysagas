from __future__ import annotations
import numpy as np
from typing import Union


class Vector:
    """An N-dimensional vector."""

    PRECISION = 3

    def __init__(self, x: float = None, y: float = None, z: float = None):
        """Define a new vector.

        Parameters
        ----------
        x : float, optional
            The x-component of the vector.

        y : float, optional
            The y-component of the vector.

        z : float, optional
            The z-component of the vector.
        """
        self._x = x
        self._y = y
        self._z = z

        self._non_none = [str(i) for i in [x, y, z] if i is not None]
        self._round_non_none = [
            str(round(i, self.PRECISION)) for i in [x, y, z] if i is not None
        ]
        self._dimensions = len(self._non_none)

    def __str__(self) -> str:
        s = f"{self._dimensions}-dimensional vector: ({', '.join(self._round_non_none)})"
        return s

    def __repr__(self) -> str:
        return f"Vector({', '.join(self._round_non_none)})"

    def __neg__(self):
        """Returns the vector pointing in the opposite direction."""
        return Vector(x=-1 * self.x, y=-1 * self.y, z=-1 * self.z)

    def __add__(self, other):
        """Element-wise vector addition.

        Parameters
        ----------
        other : Vector
            Another Vector object to be added. This Vector must be of the same
            dimension as the one it is being added to.
        """
        if not isinstance(other, Vector):
            raise Exception(f"Cannot add a {type(other)} to a vector.")
        return Vector(x=self.x + other.x, y=self.y + other.y, z=self.z + other.z)

    def __sub__(self, other):
        """Element-wise vector subtraction.

        Parameters
        ----------
        other : Vector
            Another Vector object to be added. This Vector must be of the same
            dimension as the one it is being added to.
        """
        if not isinstance(other, Vector):
            raise Exception(f"Cannot add a {type(other)} to a vector.")
        return Vector(x=self.x - other.x, y=self.y - other.y, z=self.z - other.z)

    def __truediv__(self, denominator: Union[float, int]):
        """Element-wise vector division.

        Parameters
        ----------
        denominator : float | int
            The denominator to use in the division.
        """
        return Vector(
            x=self.x / denominator, y=self.y / denominator, z=self.z / denominator
        )

    def __mul__(self, multiple: Union[float, int]):
        """Element-wise vector multiplication.

        Parameters
        ----------
        multiple : float | int
            The multiple to use in the multiplication.
        """
        return Vector(x=self.x * multiple, y=self.y * multiple, z=self.z * multiple)

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Vector):
            raise Exception(f"Cannot compare a {type(other)} to a vector.")
        return (self.x == other.x) & (self.y == other.y) & (self.z == other.z)

    @property
    def x(self) -> float:
        return self._x

    @property
    def y(self) -> float:
        return self._y

    @property
    def z(self) -> float:
        return self._z

    @property
    def vec(self) -> np.array:
        """The vector represented as a Numpy array."""
        return np.array([float(i) for i in self._non_none])

    @property
    def unit(self) -> Vector:
        """The unit vector associated with the Vector."""
        return self / self.norm

    @property
    def norm(self) -> Vector:
        """The norm associated with the Vector."""
        return np.linalg.norm(self.vec)

    @classmethod
    def from_coordinates(cls, coordinates: np.array) -> Vector:
        """Constructs a Vector object from an array of coordinates.

        Parameters
        ----------
        coordinates : np.array
            The coordinates of the vector.

        Returns
        -------
        Vector

        Examples
        --------
        >>> Vector.from_coordinates([1,2,3])
        Vector(1, 2, 3)
        """
        return cls(*coordinates)
