from .geometry.vector import Vector


class GasState:
    """An ideal gas state defined by Mach number, pressure and
    temperature.
    """

    gamma = 1.4
    R = 287  # J/kg·K.3

    def __init__(
        self, mach: float, pressure: float, temperature: float, gamma: float = 1.4
    ) -> None:
        """Define a new gas state.

        Parameters
        -----------
        mach : float
            The flow Mach number.
        pressure : float
            The flow pressure (Pa).
        temperature : float
            The flow temperature (K).
        gamma : float, optional
            The ratio of specific heats. The default is 1.4.
        """
        # Assign properties
        self._T = temperature
        self._P = pressure
        self._M = mach
        self._gamma = gamma

        # Calculate dependents
        self._rho = self.P / (self.R * self.T)
        self._a = (self.gamma * self.R * self.T) ** 0.5
        self._v = self.M * self.a
        self._q = 0.5 * self._rho * self._v**2

    def __str__(self) -> str:
        return f"Mach {self.M} flow condition with P = {self.P}, T = {self.T}."

    def __repr__(self) -> str:
        return f"Flow(M={self.M}, P={self.P}, T={self.T})"

    @property
    def T(self):
        return self._T

    @property
    def P(self):
        return self._P

    @property
    def M(self):
        return self._M

    @property
    def a(self):
        return self._a

    @property
    def rho(self):
        return self._rho

    @property
    def v(self):
        return self._v

    @property
    def q(self):
        return self._q

    @property
    def gamma(self):
        return self._gamma


class FlowState(GasState):
    """An ideal gas state defined by Mach number, pressure and
    temperature, with a flow direction.
    """

    def __init__(
        self,
        mach: float,
        pressure: float,
        temperature: float,
        direction: Vector = None,
        gamma: float = 1.4,
    ) -> None:
        """Define a new flow state.

        Parameters
        -----------
        mach : float
            The flow Mach number.
        pressure : float
            The flow pressure (Pa).
        temperature : float
            The flow temperature (K).
        direction : Vector, optional
            The direction vector of the flow. The default is Vector(1,0,0).
        gamma : float, optional
            The ratio of specific heats. The default is 1.4.
        """
        super().__init__(mach, pressure, temperature, gamma)
        if direction:
            self.direction = direction.unit
        else:
            self.direction = Vector(1, 0, 0)

        # Velocity vector
        self._Vector = self.direction * self.v

    @property
    def vx(self):
        return self._Vector.x

    @property
    def vy(self):
        return self._Vector.y

    @property
    def vz(self):
        return self._Vector.z

    @property
    def vec(self):
        return self._Vector.vec

    @property
    def Vector(self):
        return self._Vector


if __name__ == "__main__":
    flow = FlowState(6, 700, 70)
