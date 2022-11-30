class FlowState:
    """An ideal gas flow state defined by Mach number, pressure and
    temperature.
    """

    gamma = 1.4
    R = 287  # J/kg·K.3

    def __init__(self, mach: float, pressure: float, temperature: float) -> None:
        """Define a new flow state.

        Parameters
        -----------
        mach : float
            The flow Mach number.
        pressure : float
            The flow pressure (Pa).
        temperature : float
            The flow temperature (K).
        """
        # Assign properties
        self._T = temperature
        self._P = pressure
        self._M = mach

        # Calculate dependents
        self._rho = self.P / (self.R * self.T)
        self._a = (self.gamma * self.R * self.T) ** 0.5
        self._v = self.M * self.a

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


if __name__ == "__main__":
    flow = FlowState(6, 700, 70)
