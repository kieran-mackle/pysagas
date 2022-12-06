from abc import ABC, abstractmethod


class AbstractWrapper(ABC):
    @abstractmethod
    def __init__(self, **kwargs) -> None:
        pass

    @abstractmethod
    def __repr__(self) -> str:
        return "PySAGAS Solver Wrapper"

    @abstractmethod
    def __str__(self) -> str:
        return "PySAGAS Solver Wrapper"
