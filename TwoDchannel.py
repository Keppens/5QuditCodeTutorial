# -*- coding: utf-8 -*-
"""
Created on Thu Aug 29 18:18:44 2024

Two-Qudit depolarization channel

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np
import itertools
from itertools import combinations
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union, TYPE_CHECKING

class depolarizeTwoQudit(cirq.Gate):
    
    def __init__(self,p: float, d: int) -> None:
        self._d = d
        self._p = p
        error_probabilities = {}
        self.elements = self.generate_elements_with_two_subscripts()
        self.combinations = self.generate_combinations()
        p_depol = p/d**4
        p_identity = 1.0 - p*(d**4-1)/d**4
        array = ["I"]
    
    # Generate the 'Xj', 'Zj', 'Yj', 'Wjk' items for j from 1 to d-1
        for combo in self.combinations:
            array.append(combo)
        for pauli_tuple in itertools.product(array):
            pauli_string = ''.join(pauli_tuple)
            if pauli_string == 'I':
                error_probabilities[pauli_string] = p_identity
            else:
                error_probabilities[pauli_string] = p_depol
        self._error_probabilities = error_probabilities
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
    """
    # Function to generate all elements with 'j' ranging from 1 to d-1 and 'i' being either 1 or 2
    def generate_elements_with_two_subscripts(self):
        elements = []
        for j in range(1, self._d):
            for i in range(0, 2):  # 'i' can be 1 or 2
                elements.append(f"X{j}{i}")
                elements.append(f"Z{j}{i}")
                for k in range(1, self._d):
                    elements.append(f"Y{j}{k}{i}")
        return elements
    
    # Function to generate all valid combinations of one and two elements
    def generate_combinations(self):
        all_combinations = []
    
        # Single element combinations as strings
        for elem in self.elements:
            all_combinations.append(elem)
        
        # Two-element combinations as concatenated strings, ensuring different 'i' indices
        for combo in combinations(self.elements, 2):
            if combo[0][-1] != combo[1][-1]:  # Check last character ('i') is different
                combined_string = combo[0] + ' ' + combo[1]
                all_combinations.append(combined_string)
        
        return all_combinations
    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,self._d)    

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

    def _process_part(self, part):
        """
        Helper function to process a single part of the input string and return the corresponding matrix and index.
        """
        letter = part[0]

        if letter == 'X' or letter == 'Z':
            if len(part) != 3:
                raise ValueError(f"Invalid length for part '{part}'. Expected length 3.")
            i = int(part[1])
            matrix_index = int(part[2])

            if letter == 'X':
                matrix = self.create_shift_matrix(self._d, i)
            elif letter == 'Z':
                matrix = self.create_roots_of_unity_matrix(self._d, i)

        elif letter == 'Y':
            if len(part) != 4:
                raise ValueError(f"Invalid length for part '{part}'. Expected length 4.")
            i = int(part[1])
            k = int(part[2])
            matrix_index = int(part[3])
            
            matrix = self.create_shift_matrix(self._d, i) @ self.create_roots_of_unity_matrix(self._d, k)
        else:
            raise ValueError(f"Unexpected letter '{letter}' in part '{part}'.")

        return matrix, matrix_index
    
    def create_matrices_from_string(self, input_str):
        parts = input_str.split()
        
        if len(parts) == 1:
            # Single part case, use the original logic
            part = parts[0]
            matrix, matrix_index = self._process_part(part)
            
            if matrix_index == 0:
                matrix0 = matrix
                matrix1 = np.eye(self._d)
            elif matrix_index == 1:
                matrix0 = np.eye(self._d)
                matrix1 = matrix
            else:
                raise ValueError(f"Unexpected matrix index '{matrix_index}' in part '{part}'.")
        
        elif len(parts) == 2:
            # Two parts case, handle both separately
            part1, part2 = parts
            matrix0, matrix_index0 = self._process_part(part1)
            matrix1, matrix_index1 = self._process_part(part2)

            if matrix_index0 == 0:
                # Part 1 determines matrix0
                matrix0 = matrix0
            elif matrix_index0 == 1:
                # Part 1 determines matrix1
                matrix0 = np.eye(self._d)
                matrix1 = matrix0
            
            if matrix_index1 == 0:
                # Part 2 determines matrix0
                matrix0 = matrix1
            elif matrix_index1 == 1:
                # Part 2 determines matrix1
                matrix1 = matrix1

        else:
            raise ValueError(f"Unexpected number of parts in input string '{input_str}'.")

        return matrix0, matrix1
    
    def _mixture_(self) -> Sequence[Tuple[float, np.ndarray]]:
        ps = []
        for pauli in self._error_probabilities:
            Pi = np.identity(1)
            if pauli == 'I':
                    Pi = np.kron(Pi, np.eye(self._d))
                    Pi = np.kron(Pi, np.eye(self._d))
            else:
                matrix0, matrix1 = self.create_matrices_from_string(pauli)
                Pi = np.kron(Pi, matrix0)
                Pi = np.kron(Pi, matrix1)       
            ps.append(Pi)
        return tuple(zip(self._error_probabilities.values(), ps))
    
    def _has_mixture_(self) -> bool:
        return True


    def _circuit_diagram_info_(self, args):
        return f"D2({self._p})"

