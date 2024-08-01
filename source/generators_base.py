"""
Conatins the definition of the ABC Generator. A generator is required for the new ligand generation step. Every Generator implements the generate(...) function which contains more or less the logic of ligand generation.
"""

from abc import ABC, abstractmethod

class Generator(ABC):
    @abstractmethod
    def generate(self, *args, **kwargs):
        pass