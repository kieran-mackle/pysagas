import pandas as pd
from typing import Optional
from abc import ABC, abstractmethod
from pysagas.cfd.solver import FlowResults, SensitivityResults


class AbstractDeck(ABC):
    @abstractmethod
    def __init__(self) -> None:
        pass

    @abstractmethod
    def __repr__(self) -> None:
        pass

    @abstractmethod
    def __str__(self) -> None:
        pass

    @property
    @abstractmethod
    def deck(self) -> pd.DataFrame:
        pass

    @abstractmethod
    def insert(self):
        pass

    @abstractmethod
    def to_csv(self):
        pass

    # @abstractmethod
    # def interpolate(self, **kwargs):
    #     pass


class Deck(AbstractDeck):
    TYPE = "deck"
    _COLS = []

    def __init__(
        self, inputs: list[str], a_ref: Optional[float] = 1, c_ref: Optional[float] = 1
    ) -> None:
        """Instantiate a new deck.

        Parameters
        ----------
        inputs : list[str]
            A list of the inputs to this aerodeck. For example, ["aoa", "mach"].

        a_ref : float, optional
            The reference area. The default is 1.

        c_ref : float, optional
            The reference length. The default is 1.
        """
        self._deck = pd.DataFrame(columns=inputs + self._COLS)
        self._inputs = inputs
        self._a_ref = a_ref
        self._c_ref = c_ref

    def __repr__(self) -> None:
        return f"{self.TYPE}"

    def __str__(self) -> None:
        return self.__repr__()

    @property
    def deck(self) -> pd.DataFrame:
        return self._deck.drop_duplicates()

    def to_csv(self, file_prefix: Optional[str] = None):
        """Save the deck to CSV file.

        Parameters
        ----------
        file_prefix : str, optional
            The CSV file name prefix. If None is provided, the deck __repr__
            will be used. The default is None.
        """
        file_prefix = file_prefix if file_prefix else self.__repr__()
        self.deck.to_csv(f"{file_prefix}.csv", index=False)


class Aerodeck(Deck):
    """Aerodynamic coefficient deck."""

    TYPE = "aerodeck"
    _COLS = ["CL", "CD", "Cm"]

    def insert(self, result: FlowResults, **kwargs):
        # Check inputs
        inputs_given = [i in kwargs for i in self._inputs]
        if not all(inputs_given):
            raise Exception(
                "Please provide all input values when inserting new result."
            )

        # Get coefficients
        coefficients = result.coefficients(A_ref=self._a_ref, c_ref=self._c_ref)

        # Extract data
        data = {
            "CL": coefficients[0],
            "CD": coefficients[1],
            "Cm": coefficients[2],
        }

        # Add inputs
        data.update(kwargs)

        # Add to deck
        self._deck = pd.concat(
            [self._deck, pd.DataFrame(data, index=[len(self._deck)])]
        )

    @classmethod
    def from_csv(cls, filepath: str, **kwargs):
        """Create an Aerodeck from a CSV file."""
        # Read data from file
        data = pd.read_csv(filepath)

        # Extract inputs
        inputs = list(data.columns)
        [inputs.pop(inputs.index(coef)) for coef in Aerodeck._COLS]

        # Instantiate aerodeck
        aerodeck = cls(inputs=inputs, **kwargs)
        aerodeck._deck = data
        return aerodeck


class Sensdeck(Deck):
    """Aerodynamic coefficient sensitivity deck."""

    TYPE = "sensdeck"
    _COLS = ["dCL", "dCD", "dCm"]

    def __init__(
        self,
        inputs: list[str],
        parameters: list[str],
        a_ref: Optional[float] = 1,
        c_ref: Optional[float] = 1,
    ) -> None:
        """Instantiate a new sensitivity deck.

        Parameters
        ----------
        inputs : list[str]
            A list of the inputs to this deck. For example, ["aoa", "mach"].

        parameters : list[str]
            A list of the parameters to this deck. For example, ["wingspan", "length"].

        a_ref : float, optional
            The reference area. The default is 1.

        c_ref : float, optional
            The reference length. The default is 1.
        """
        super().__init__(inputs, a_ref, c_ref)
        self._parameters = list(parameters)

        # Override self._deck
        base_deck = pd.DataFrame(columns=inputs + Sensdeck._COLS)
        self._deck = {p: base_deck.copy() for p in self._parameters}

    def insert(self, result: SensitivityResults, **kwargs):
        # Check inputs
        inputs_given = [i in kwargs for i in self._inputs]
        if not all(inputs_given):
            raise Exception(
                "Please provide all input values when inserting new result."
            )

        # Get coefficients
        f_sens, m_sens = result.coefficients(A_ref=self._a_ref, c_ref=self._c_ref)

        for param in self._parameters:
            # Extract data for this parameter
            data = {
                "dCL": f_sens.loc[param]["dCL"],
                "dCD": f_sens.loc[param]["dCD"],
                "dCm": m_sens.loc[param]["dMz/dp"],
            }

            # Add inputs
            data.update(kwargs)

            # Add to deck
            self._deck[param] = pd.concat(
                [self._deck[param], pd.DataFrame(data, index=[len(self._deck)])]
            )

    @property
    def deck(self) -> dict[str, pd.DataFrame]:
        decks = {p: df.drop_duplicates() for p, df in self._deck.items()}
        return decks

    def to_csv(self, file_prefix: Optional[str] = None):
        """Save the deck to CSV file.

        Parameters
        ----------
        file_prefix : str, optional
            The CSV file name prefix. If None is provided, the deck __repr__
            will be used. The default is None.
        """
        decks = self.deck
        file_prefix = file_prefix if file_prefix else self.__repr__()

        for p, df in decks.items():
            df.to_csv(f"{p}_{file_prefix}.csv", index=False)

    @classmethod
    def from_csv(cls, param_filepaths: dict[str, str], **kwargs):
        """Create an Aerodeck from a CSV file.

        Parameters
        ----------
        param_filepaths : dict[str, str]
            A dictionary of filepaths, keyed by the associated parameter.
        """
        decks = {}
        for param, filepath in param_filepaths.items():
            # Load data
            data = pd.read_csv(filepath)
            decks[param] = data

        # Extract inputs
        inputs = list(data.columns)
        [inputs.pop(inputs.index(coef)) for coef in Sensdeck._COLS]

        # Instantiate sensdeck
        sensdeck = cls(inputs=inputs, parameters=param_filepaths.keys(), **kwargs)
        sensdeck._deck = decks

        return sensdeck
