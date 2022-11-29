import numpy as np
from typing import Union


class Vector:
    """An N-dimensional vector."""

    PRECISION = 3

    def __init__(self, x: float = None, y: float = None, z: float = None):
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
        """Vector division.

        Parameters
        ----------
        denominator : float | int
            The denominator to use in the division.
        """
        return Vector(
            x=self.x / denominator, y=self.y / denominator, z=self.z / denominator
        )

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
        """The vector as a Numpy array."""
        return np.array([float(i) for i in self._non_none])


class Cell:
    """A cell."""

    def __init__(self, p0: Vector, p1: Vector, p2: Vector) -> None:
        """Constructs a cell, defined by three points.

        Parameters
        ----------
        p0 : Vector
            The first point defining the cell.
        p1 : Vector
            The second point defining the cell.
        p2 : Vector
            The third point defining the cell.
        """
        # Save points
        self.p0 = p0
        self.p1 = p1
        self.p2 = p2

        # Calculate normal vector
        self.n = self.calc_normal(p0, p1, p2)
        self.A = self.calc_area(p0, p1, p2)

        # Calculate geometric sensitivities
        # TODO - below not validated
        _, self.dn = create_sensitivity_matrix(self.p0, self.p1, self.p2)
        _, self.dA = create_area_sensitivity_matrix(self.p0, self.p1, self.p2)

    def __repr__(self) -> str:
        return f"Cell({self.p0}, {self.p1}, {self.p2})"

    def __str__(self) -> str:
        return "A Cell"

    @staticmethod
    def calc_normal(p0: Vector, p1: Vector, p2: Vector) -> Vector:
        """Calculates the normal vector of a cell defined by three points.

        Parameters
        ----------
        p0 : Vector
            The first point defining the cell.
        p1 : Vector
            The second point defining the cell.
        p2 : Vector
            The third point defining the cell.

        Returns
        --------
        normal : Vector
            The unit normal vector of the cell defined by the points.

        References
        -----------
        https://www.khronos.org/opengl/wiki/Calculating_a_Surface_Normal
        """
        # Define cell edges
        A = p1 - p0
        B = p2 - p0

        # Calculate normal components
        Nx = A.y * B.z - A.z * B.y
        Ny = A.z * B.x - A.x * B.z
        Nz = A.x * B.y - A.y * B.x

        # Construct normal vector
        normal = Vector(x=Nx, y=Ny, z=Nz)

        # Convert to unit vector
        unit_normal = normal / np.linalg.norm(normal.vec)

        return unit_normal

    @staticmethod
    def calc_area(p0: Vector, p1: Vector, p2: Vector) -> float:
        """Calculates the area of a cell defined by three points.

        Parameters
        ----------
        p0 : Vector
            The first point defining the cell.
        p1 : Vector
            The second point defining the cell.
        p2 : Vector
            The third point defining the cell.

        Returns
        --------
        area : float
            The area of the cell defined by the points.

        References
        -----------
        https://en.wikipedia.org/wiki/Cross_product
        """
        # Define cell edges
        A = p1 - p0
        B = p2 - p0

        # Calculate area
        S = 0.5 * np.linalg.norm(np.cross(A.vec, B.vec))
        return S


def calculate_3d_normal(p0, p1, p2):
    """Calculates the normal vector of a plane defined by 3 points."""
    n = np.cross(
        np.array([p1.x - p0.x, p1.y - p0.y, p1.z - p0.z]),
        np.array([p2.x - p0.x, p2.y - p0.y, p2.z - p0.z]),
    )
    n = n / np.sqrt(np.sum(n**2))
    return n


def create_sensitivity_matrix(A: Vector, B: Vector, C: Vector):
    """Calculates the sensitivity of face normals to three points, A, B and C."""
    g = (
        ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) ** 2
        + ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) ** 2
        + ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) ** 2
    )

    n = np.empty(3)
    n[0] = ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) / g**0.5
    n[1] = ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) / g**0.5
    n[2] = ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) / g**0.5

    M_sense = np.empty((3, 9))
    for row in range(3):
        if row == 0:
            f = (B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)  # f_x
        elif row == 1:
            f = (B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)  # f_y
        elif row == 2:
            f = (B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)  # f_z
        else:
            raise Exception("Not sure how I got here.")

        for col in range(9):
            if col == 0:  # d/dA.x
                g_dash = 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                    -(B.z - A.z) + (C.z - A.z)
                ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                    -(C.y - A.y) + (B.y - A.y)
                )
            elif col == 1:  # d/dA.y
                g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                    -(C.z - A.z) + (B.z - A.z)
                ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                    -(B.x - A.x) + (C.x - A.x)
                )
            elif col == 2:  # d/dA.z
                g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                    -(B.y - A.y) + (C.y - A.y)
                ) + 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                    -(C.x - A.x) + (B.x - A.x)
                )
            elif col == 3:  # d/dB.x
                g_dash = 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                    -(C.z - A.z)
                ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                    (C.y - A.y)
                )
            elif col == 4:  # d/dB.y
                g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                    (C.z - A.z)
                ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                    -(C.x - A.x)
                )
            elif col == 5:  # d/dB.z
                g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                    -(C.y - A.y)
                ) + 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                    (C.x - A.x)
                )
            elif col == 6:  # d/dC.x
                g_dash = 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                    (B.z - A.z)
                ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                    -(B.y - A.y)
                )
            elif col == 7:  # d/dC.y
                g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                    -(B.z - A.z)
                ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                    (B.x - A.x)
                )
            elif col == 8:  # d/dC.z
                g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                    (B.y - A.y)
                ) + 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                    -(B.x - A.x)
                )
            else:
                raise Exception("???")

            if row == 0:  # for n..x
                # ((B.y-A.y)*(C.z-A.z) - (B.z-A.z)*(C.y-A.y))  # f_x
                if col == 0:  # d/dA.x
                    f_dash = 0
                elif col == 1:  # d/dA.y
                    f_dash = -(C.z - A.z) + (B.z - A.z)
                elif col == 2:  # d/dA.z
                    f_dash = -(B.y - A.y) + (C.y - A.y)
                elif col == 3:  # d/dB.x
                    f_dash = 0
                elif col == 4:  # d/dB.y
                    f_dash = C.z - A.z
                elif col == 5:  # d/dB.z
                    f_dash = -(C.y - A.y)
                elif col == 6:  # d/dC.x
                    f_dash = 0
                elif col == 7:  # d/dC.y
                    f_dash = -(B.z - A.z)
                elif col == 8:  # d/dC.z
                    f_dash = B.y - A.y
            if row == 1:  # for n.y  # CHECK
                # f = ((B.z-A.z)*(C.x-A.x) - (B.x-A.x)*(C.z-A.z))  # f_y
                if col == 0:  # d/dA.x
                    f_dash = -(B.z - A.z) + (C.z - A.z)
                elif col == 1:  # d/dA.y
                    f_dash = 0
                elif col == 2:  # d/dA.z
                    f_dash = -(C.x - A.x) + (B.x - A.x)  # CHECK
                elif col == 3:  # d/dB.x
                    f_dash = -(C.z - A.z)
                elif col == 4:  # d/dB.y
                    f_dash = 0
                elif col == 5:  # d/dB.z
                    f_dash = C.x - A.x
                elif col == 6:  # d/dC.x
                    f_dash = B.z - A.z
                elif col == 7:  # d/dC.y
                    f_dash = 0
                elif col == 8:  # d/dC.z
                    f_dash = -(B.x - A.x)
            if row == 2:  # for n.z  # CHECK
                # f = ((B.x-A.x)*(C.y-A.y) - (B.y-A.y)*(C.x-A.x))  # f_z
                if col == 0:  # d/dA.x
                    f_dash = -(C.y - A.y) + (B.y - A.y)
                elif col == 1:  # d/dA.y
                    f_dash = -(B.x - A.x) + (C.x - A.x)
                elif col == 2:  # d/dA.z
                    f_dash = 0
                elif col == 3:  # d/dB.x
                    f_dash = C.y - A.y
                elif col == 4:  # d/dB.y
                    f_dash = -(C.x - A.x)
                elif col == 5:  # d/dB.z
                    f_dash = 0
                elif col == 6:  # d/dC.x
                    f_dash = -(B.y - A.y)
                elif col == 7:  # d/dC.y
                    f_dash = B.x - A.x
                elif col == 8:  # d/dC.z
                    f_dash = 0

            M_sense[row, col] = (
                f_dash * (g ** (-0.5)) + f * (-1 / 2) * g ** (-3 / 2) * g_dash
            )
    return n, M_sense


def create_area_sensitivity_matrix(A: Vector, B: Vector, C: Vector):
    """Calculates the sensitivity of a cell's area to its corner points."""
    g = (
        ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) ** 2
        + ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) ** 2
        + ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) ** 2
    )  # g

    A_sense = np.empty(9)

    for col in range(9):
        if col == 0:  # d/dA.x
            g_dash = 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                -(B.z - A.z) + (C.z - A.z)
            ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                -(C.y - A.y) + (B.y - A.y)
            )
        elif col == 1:  # d/dA.y
            g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                -(C.z - A.z) + (B.z - A.z)
            ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                -(B.x - A.x) + (C.x - A.x)
            )
        elif col == 2:  # d/dA.z
            g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                -(B.y - A.y) + (C.y - A.y)
            ) + 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                -(C.x - A.x) + (B.x - A.x)
            )
        elif col == 3:  # d/dB.x
            g_dash = 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                -(C.z - A.z)
            ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                (C.y - A.y)
            )
        elif col == 4:  # d/dB.y
            g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                (C.z - A.z)
            ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                -(C.x - A.x)
            )
        elif col == 5:  # d/dB.z
            g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                -(C.y - A.y)
            ) + 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                (C.x - A.x)
            )
        elif col == 6:  # d/dC.x
            g_dash = 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                (B.z - A.z)
            ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                -(B.y - A.y)
            )
        elif col == 7:  # d/dC.y
            g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                -(B.z - A.z)
            ) + 2 * ((B.x - A.x) * (C.y - A.y) - (B.y - A.y) * (C.x - A.x)) * (
                (B.x - A.x)
            )
        elif col == 8:  # d/dC.z
            g_dash = 2 * ((B.y - A.y) * (C.z - A.z) - (B.z - A.z) * (C.y - A.y)) * (
                (B.y - A.y)
            ) + 2 * ((B.z - A.z) * (C.x - A.x) - (B.x - A.x) * (C.z - A.z)) * (
                -(B.x - A.x)
            )

        A_sense[col] = 0.5 * (1 / 2) * g ** (-1 / 2) * g_dash

    Area = 0.5 * g**0.5

    return Area, A_sense
