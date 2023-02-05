from __future__ import annotations
import numpy as np
import pysagas.flow
from typing import Union, List
from numpy.typing import ArrayLike


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


class Cell:
    """A triangular cell object.

    Attributes
    -----------
    p0 : Vector
        The first vertex of the cell.
    p1 : Vector
        The second vertex of the cell.
    p2 : Vector
        The third vertex of the cell.
    A : float
        The cell face area.
    n : Vector
        The cell normal.
    dndv : np.array
        The sensitivity of the cell's normal vector to each vertex.
    dAdv : np.array
        The sensitivity of the cell's area to each vertex.
    dvdp : np.array
        The sensitivity of the cell's vertices to each geometric parameter.
    dndp : np.array
        The sensitivity of the cell's normal vector to each geometric parameter.
    dAdp : np.array
        The sensitivity of the cell's area to each geometric parameter.
    dcdp : np.array
        The sensitivity of the centroid to each geometric parameter.
    flowstate : FlowState
        The flowstate associated with the cell.
    sensitivities : np.array
        An array containing the [x,y,z] force sensitivities of the cell.
    """

    def __init__(self, p0: Vector, p1: Vector, p2: Vector):
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

        # Calculate cell properites
        self.n = -self.calc_normal(p0, p1, p2)
        self.A = self.calc_area(p0, p1, p2)
        self.c = self.calc_centroid(p0, p1, p2)

        # Calculate geometric sensitivities
        self.dndv = self.n_sensitivity(self.p0, self.p1, self.p2)
        self.dAdv = self.A_sensitivity(self.p0, self.p1, self.p2)
        self.dcdv = self.c_sensitivity(self.p0, self.p1, self.p2)

        # Parameter sensitivities
        self.dvdp: ArrayLike = None  # vertex-parameter sensitivities
        self.dndp: ArrayLike = None  # normal-parameter sensitivities
        self.dAdp: ArrayLike = None  # area-parameter sensitivities
        self.dcdp: ArrayLike = None  # centroid-parameter sensitivities

        # FlowState
        self.flowstate: pysagas.flow.FlowState = None

        # Sensitivities
        self.sensitivities = None
        self.moment_sensitivities = None

    def to_dict(self):
        """Returns the Cell as a dictionary."""
        # TODO - need to think about representing the cells
        # vs. the points
        # The values at the points will need to be averaged
        # across cells
        pass

    @classmethod
    def from_points(cls, points: Union[List[Vector], np.array[Vector]]) -> Cell:
        """Constructs a Vector object from an array of coordinates.

        Parameters
        ----------
        points : Union[List[Vector], np.array[Vector]]
            The points defining the cell.

        Returns
        -------
        Cell
        """
        return cls(*points)

    def __repr__(self) -> str:
        return f"Cell({self.p0}, {self.p1}, {self.p2})"

    def __str__(self) -> str:
        return "A Cell"

    def _add_sensitivities(self, dvdp: np.array) -> None:
        """Adds the cell sensitivities to the cell."""
        self.dvdp = dvdp
        self.dndp = np.dot(self.dndv, self.dvdp)
        self.dAdp = np.dot(self.dAdv, self.dvdp)
        self.dcdp = np.dot(self.dcdv, self.dvdp)

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

        if normal.norm == 0:
            # Degenerate cell
            raise DegenerateCell()

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

    @staticmethod
    def calc_centroid(p0: Vector, p1: Vector, p2: Vector) -> Vector:
        """Calculates the centroid of a cell defined by three points.

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
        c : Vector
            The centroid of the cell defined by the points.

        References
        -----------
        https://en.wikipedia.org/wiki/Centroid
        """
        cx = (p0.x + p1.x + p2.x) / 3
        cy = (p0.y + p1.y + p2.y) / 3
        cz = (p0.z + p1.z + p2.z) / 3
        return Vector(cx, cy, cz)

    @staticmethod
    def n_sensitivity(p0: Vector, p1: Vector, p2: Vector) -> np.array:
        """Calculates the sensitivity of a cell's normal vector
        to the points defining the cell analytically.

        Parameters
        ----------
        p0 : Vector
            The first point defining the cell.
        p1 : Vector
            The second point defining the cell.
        p2 : Vector
            The third point defining the cell.

        Returns
        -------
        sensitivity : np.array
            The sensitivity matrix with size m x n x p. Rows m refer to
            the vertices, columns n refer to the vertex coordinates, and
            slices p refer to the components of the normal vector.
        """
        g = (
            ((p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)) ** 2
            + ((p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)) ** 2
            + ((p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)) ** 2
        )

        # n = np.empty(3)
        # A = p0
        # B = p1
        # C = p2
        # n[0] = ((B.y-A.y)*(C.z-A.z) - (B.z-A.z)*(C.y-A.y)) / g**0.5
        # n[1] = ((B.z-A.z)*(C.x-A.x) - (B.x-A.x)*(C.z-A.z)) / g**0.5
        # n[2] = ((B.x-A.x)*(C.y-A.y) - (B.y-A.y)*(C.x-A.x)) / g**0.5

        M_sense = np.empty((3, 9))
        for row in range(3):
            if row == 0:
                f = (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)  # f_x
            elif row == 1:
                f = (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)  # f_y
            elif row == 2:
                f = (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)  # f_z

            for col in range(9):
                if col == 0:  # d/dp0.x
                    g_dash = 2 * (
                        (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                    ) * (-(p1.z - p0.z) + (p2.z - p0.z)) + 2 * (
                        (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                    ) * (
                        -(p2.y - p0.y) + (p1.y - p0.y)
                    )
                elif col == 1:  # d/dp0.y
                    g_dash = 2 * (
                        (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                    ) * (-(p2.z - p0.z) + (p1.z - p0.z)) + 2 * (
                        (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                    ) * (
                        -(p1.x - p0.x) + (p2.x - p0.x)
                    )
                elif col == 2:  # d/dp0.z
                    g_dash = 2 * (
                        (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                    ) * (-(p1.y - p0.y) + (p2.y - p0.y)) + 2 * (
                        (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                    ) * (
                        -(p2.x - p0.x) + (p1.x - p0.x)
                    )
                elif col == 3:  # d/dp1.x
                    g_dash = 2 * (
                        (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                    ) * (-(p2.z - p0.z)) + 2 * (
                        (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                    ) * (
                        (p2.y - p0.y)
                    )
                elif col == 4:  # d/dp1.y
                    g_dash = 2 * (
                        (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                    ) * ((p2.z - p0.z)) + 2 * (
                        (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                    ) * (
                        -(p2.x - p0.x)
                    )
                elif col == 5:  # d/dp1.z
                    g_dash = 2 * (
                        (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                    ) * (-(p2.y - p0.y)) + 2 * (
                        (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                    ) * (
                        (p2.x - p0.x)
                    )
                elif col == 6:  # d/dp2.x
                    g_dash = 2 * (
                        (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                    ) * ((p1.z - p0.z)) + 2 * (
                        (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                    ) * (
                        -(p1.y - p0.y)
                    )
                elif col == 7:  # d/dp2.y
                    g_dash = 2 * (
                        (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                    ) * (-(p1.z - p0.z)) + 2 * (
                        (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                    ) * (
                        (p1.x - p0.x)
                    )
                elif col == 8:  # d/dp2.z
                    g_dash = 2 * (
                        (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                    ) * ((p1.y - p0.y)) + 2 * (
                        (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                    ) * (
                        -(p1.x - p0.x)
                    )

                if row == 0:  # for n..x
                    # ((p1.y-p0.y)*(p2.z-p0.z) - (p1.z-p0.z)*(p2.y-p0.y))  # f_x
                    if col == 1:  # d/dp0.y
                        f_dash = -(p2.z - p0.z) + (p1.z - p0.z)
                    elif col == 2:  # d/dp0.z
                        f_dash = -(p1.y - p0.y) + (p2.y - p0.y)
                    elif col == 4:  # d/dp1.y
                        f_dash = p2.z - p0.z
                    elif col == 5:  # d/dp1.z
                        f_dash = -(p2.y - p0.y)
                    elif col == 7:  # d/dp2.y
                        f_dash = -(p1.z - p0.z)
                    elif col == 8:  # d/dp2.z
                        f_dash = p1.y - p0.y
                    else:
                        f_dash = 0
                if row == 1:  # for n.y  # CHECK
                    # f = ((p1.z-p0.z)*(p2.x-p0.x) - (p1.x-p0.x)*(p2.z-p0.z))  # f_y
                    if col == 0:  # d/dp0.x
                        f_dash = -(p1.z - p0.z) + (p2.z - p0.z)
                    elif col == 2:  # d/dp0.z
                        f_dash = -(p2.x - p0.x) + (p1.x - p0.x)  # CHECK
                    elif col == 3:  # d/dp1.x
                        f_dash = -(p2.z - p0.z)
                    elif col == 5:  # d/dp1.z
                        f_dash = p2.x - p0.x
                    elif col == 6:  # d/dp2.x
                        f_dash = p1.z - p0.z
                    elif col == 8:  # d/dp2.z
                        f_dash = -(p1.x - p0.x)
                    else:
                        f_dash = 0
                if row == 2:  # for n.z  # CHECK
                    # f = ((p1.x-p0.x)*(p2.y-p0.y) - (p1.y-p0.y)*(p2.x-p0.x))  # f_z
                    if col == 0:  # d/dp0.x
                        f_dash = -(p2.y - p0.y) + (p1.y - p0.y)
                    elif col == 1:  # d/dp0.y
                        f_dash = -(p1.x - p0.x) + (p2.x - p0.x)
                    elif col == 3:  # d/dp1.x
                        f_dash = p2.y - p0.y
                    elif col == 4:  # d/dp1.y
                        f_dash = -(p2.x - p0.x)
                    elif col == 6:  # d/dp2.x
                        f_dash = -(p1.y - p0.y)
                    elif col == 7:  # d/dp2.y
                        f_dash = p1.x - p0.x
                    else:
                        f_dash = 0

                # Assign sensitivity
                M_sense[row, col] = (
                    f_dash * (g ** (-0.5)) + f * (-1 / 2) * g ** (-3 / 2) * g_dash
                )
        return M_sense

    @staticmethod
    def A_sensitivity(p0: Vector, p1: Vector, p2: Vector) -> np.array:
        """Calculates the sensitivity of a cell's area to the points
        defining the cell analytically.

        Parameters
        ----------
        p0 : Vector
            The first point defining the cell.
        p1 : Vector
            The second point defining the cell.
        p2 : Vector
            The third point defining the cell.

        Returns
        -------
        sensitivity : np.array
            The sensitivity matrix with size m x n. Rows m refer to
            the vertices, columns n refer to the vertex coordinates.
        """
        g = (
            ((p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)) ** 2
            + ((p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)) ** 2
            + ((p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)) ** 2
        )

        A_sense = np.empty(9)
        for col in range(9):
            if col == 0:  # d/dp0.x
                g_dash = 2 * (
                    (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                ) * (-(p1.z - p0.z) + (p2.z - p0.z)) + 2 * (
                    (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                ) * (
                    -(p2.y - p0.y) + (p1.y - p0.y)
                )
            elif col == 1:  # d/dp0.y
                g_dash = 2 * (
                    (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                ) * (-(p2.z - p0.z) + (p1.z - p0.z)) + 2 * (
                    (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                ) * (
                    -(p1.x - p0.x) + (p2.x - p0.x)
                )
            elif col == 2:  # d/dp0.z
                g_dash = 2 * (
                    (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                ) * (-(p1.y - p0.y) + (p2.y - p0.y)) + 2 * (
                    (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                ) * (
                    -(p2.x - p0.x) + (p1.x - p0.x)
                )
            elif col == 3:  # d/dp1.x
                g_dash = 2 * (
                    (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                ) * (-(p2.z - p0.z)) + 2 * (
                    (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                ) * (
                    (p2.y - p0.y)
                )
            elif col == 4:  # d/dp1.y
                g_dash = 2 * (
                    (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                ) * ((p2.z - p0.z)) + 2 * (
                    (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                ) * (
                    -(p2.x - p0.x)
                )
            elif col == 5:  # d/dp1.z
                g_dash = 2 * (
                    (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                ) * (-(p2.y - p0.y)) + 2 * (
                    (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                ) * (
                    (p2.x - p0.x)
                )
            elif col == 6:  # d/dp2.x
                g_dash = 2 * (
                    (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                ) * ((p1.z - p0.z)) + 2 * (
                    (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                ) * (
                    -(p1.y - p0.y)
                )
            elif col == 7:  # d/dp2.y
                g_dash = 2 * (
                    (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                ) * (-(p1.z - p0.z)) + 2 * (
                    (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x)
                ) * (
                    (p1.x - p0.x)
                )
            elif col == 8:  # d/dp2.z
                g_dash = 2 * (
                    (p1.y - p0.y) * (p2.z - p0.z) - (p1.z - p0.z) * (p2.y - p0.y)
                ) * ((p1.y - p0.y)) + 2 * (
                    (p1.z - p0.z) * (p2.x - p0.x) - (p1.x - p0.x) * (p2.z - p0.z)
                ) * (
                    -(p1.x - p0.x)
                )

            A_sense[col] = 0.5 * (1 / 2) * g ** (-1 / 2) * g_dash

        return A_sense

    @staticmethod
    def c_sensitivity(p0: Vector, p1: Vector, p2: Vector) -> np.array:
        """Calculates the sensitivity of a cell's centroid to the
        points defining the cell.

        Parameters
        ----------
        p0 : Vector
            The first point defining the cell.
        p1 : Vector
            The second point defining the cell.
        p2 : Vector
            The third point defining the cell.

        Returns
        -------
        sensitivity : np.array
            The sensitivity matrix with size m x n x p. Rows m refer to
            the vertices, columns n refer to the vertex coordinates, and
            slices p refer to the components of the centroid point.
        """
        c_sens = np.array(
            [
                [1 / 3, 0, 0, 1 / 3, 0, 0, 1 / 3, 0, 0],
                [0, 1 / 3, 0, 0, 1 / 3, 0, 0, 1 / 3, 0],
                [0, 0, 1 / 3, 0, 0, 1 / 3, 0, 0, 1 / 3],
            ]
        )
        return c_sens


def calculate_3d_normal(p0: Vector, p1: Vector, p2: Vector) -> np.array:
    """Calculates the normal vector of a plane defined by 3 points.

    Parameters
    ----------
    p0 : Vector
        The first point defining the cell.
    p1 : Vector
        The second point defining the cell.
    p2 : Vector
        The third point defining the cell.

    Returns
    -------
    n : np.array
        The normal vector.
    """
    n = np.cross(
        np.array([p1.x - p0.x, p1.y - p0.y, p1.z - p0.z]),
        np.array([p2.x - p0.x, p2.y - p0.y, p2.z - p0.z]),
    )
    n = n / np.sqrt(np.sum(n**2))
    return n


class DegenerateCell(Exception):
    """Exception raised for degenerate cells."""

    def __init__(self, message="Degenerate cell"):
        self.message = message
        super().__init__(self.message)
