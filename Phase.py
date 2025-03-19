# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 16:41:58 2023

Qudit Phase gate

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np

class Phase(cirq.Gate):
    
    def __init__(self, d: int, b: int) -> None:
        self._d = d
        self._b = b
    
    """A gate that enacts the transformation U|x〉 = w^d|x〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)
    
    def _validate_args(self, qubits):
        return True 

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

    def _unitary_(self):
        # create the unitary matrix
        return np.array(self.create_roots_of_unity_matrix(self._d,self._b),dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return f'[Z({self._b})]'

class Phasedag(cirq.Gate):
    
    def __init__(self, d: int, b: int) -> None:
        self._d = d
        self._b = b
    
    """A gate that enacts the transformation U|x〉 = w^d|x〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)
    
    def _validate_args(self, qubits):
        return True 

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

    def _unitary_(self):
        # create the unitary matrix
        return np.array(np.conjugate(self.create_roots_of_unity_matrix(self._d,self._b)).T,dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return f'[Z*({self._b})]'
    

class Pg(cirq.Gate):
    
    def __init__(self, d: int, g:int) -> None:
        self._d = d
        self._g = g
    
    """A gate that enacts the transformation U|x〉 = w^d|x〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)
    
    def _validate_args(self, qubits):
        return True 

    def create_diagonal_matrix(self, d, g):
        # Calculate the d'th root of unity
        w = np.exp(2j * np.pi / d)
    
        # Create the diagonal matrix
        diagonal_matrix = np.diag([w**(i**2 * g / 2) for i in range(d)])
    
        return diagonal_matrix

    def _unitary_(self):
        # create the unitary matrix
        return np.array(self.create_diagonal_matrix(self._d, self._g),dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return '[Pγ]'

class Pgdag(cirq.Gate):
    
    def __init__(self, d: int, g:int) -> None:
        self._d = d
        self._g = g
    
    """A gate that enacts the transformation U|x〉 = w^d|x〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)
    
    def _validate_args(self, qubits):
        return True 

    def create_diagonal_matrix(self, d, g):
        # Calculate the d'th root of unity
        w = np.exp(2j * np.pi / d)
    
        # Create the diagonal matrix
        diagonal_matrix = np.diag([w**(i**2 * g / 2) for i in range(d)])
    
        return diagonal_matrix

    def _unitary_(self):
        # create the unitary matrix
        return np.array(np.conjugate(self.create_diagonal_matrix(self._d, self._g)).T,dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return '[Pγ-]'
    