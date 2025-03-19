# -*- coding: utf-8 -*-
"""
Created on Mon Oct  9 15:57:03 2023

Qudit depolarization channel

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np
import itertools
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union, TYPE_CHECKING

class depolarizeQudit(cirq.Gate):
    
    def __init__(self,p: float, d: int) -> None:
        self._d = d
        self._p = p
        error_probabilities = {}

        p_depol = p/(d**2) 
        p_identity = 1.0 - p*(d**2-1)/d**2
        array = ["I"]
    
    # Generate the 'Xj', 'Zj', 'Yj', 'Wjk' items for j from 1 to d-1
        for j in range(1, d):
            array.append(f"X{j}")
            array.append(f"Z{j}")
            for k in range(1, d):
                        array.append(f"Y{j}{k}")
        for pauli_tuple in itertools.product(array):
            pauli_string = ''.join(pauli_tuple)
            if pauli_string == 'I':
                error_probabilities[pauli_string] = p_identity
            else:
                error_probabilities[pauli_string] = p_depol
        self._error_probabilities = error_probabilities
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)    

    def create_shift_matrix(self, d, a):
        if d < 2:
            raise ValueError("Dimension 'd' must be at least 2.")
        if a > d:
            raise ValueError("Shift cannot be larger than the Dimension 'd'")
        # Perform the shift operation (x -> x+1 mod d) on the matrix
        shift_matrix = np.roll(np.eye(d), -1*a, axis=0)
        
        return shift_matrix

    def create_roots_of_unity_matrix(self, d, b):
        if d < 1:
            raise ValueError("Dimension 'd' must be at least 1.")
        if b > d:
            raise ValueError("b cannot be larger than the Dimension 'd'")
        # Compute the d squared roots of unity
        roots = [np.exp(2j * np.pi * ((k*b)%d) / d) for k in range(d)]
        
        # Create a diagonal matrix with the roots of unity
        roots_matrix = np.diag(roots)
        
        return roots_matrix
    
    def _mixture_(self) -> Sequence[Tuple[float, np.ndarray]]:
        ps = []
        for pauli in self._error_probabilities:
            Pi = np.identity(1)
            if pauli == 'I':
                    Pi = np.kron(Pi, np.eye(self._d))
            else:
                for i in range(1, self._d):
                    if pauli == f'X{i}':
                        Pi = np.kron(Pi, self.create_shift_matrix(self._d, i))
                        break
                    elif pauli == f'Z{i}':
                        Pi = np.kron(Pi, self.create_roots_of_unity_matrix(self._d, i))
                        break
                    for k in range(1, self._d):
                            if pauli == f'Y{i}{k}':
                                Pi = np.kron(Pi, self.create_shift_matrix(self._d, i) @ self.create_roots_of_unity_matrix(self._d, k))
                                break
            ps.append(Pi)
        return tuple(zip(self._error_probabilities.values(), ps))
    
    def _has_mixture_(self) -> bool:
        return True


    def _circuit_diagram_info_(self, args):
        return f"D({self._p})"

