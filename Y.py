# -*- coding: utf-8 -*-
"""
Created on Thu 10 apr 16:41:58 2025

Qudit Y gate

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np

class Y(cirq.Gate):
    
    def __init__(self, d: int, a: int, b: int) -> None:
        self._d = d
        self._a = a
        self._b = b
    
    """A gate that enacts the Y_q gate on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)

    def _validate_args(self, qubits):
        return True 

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

    def compute_Y(self):
        X = np.array(self.create_shift_matrix(self._d,self._a),dtype=np.complex64)
        Z = np.array(self.create_roots_of_unity_matrix(self._d,self._b),dtype=np.complex64)
        return X @ Z

    def _unitary_(self):
        # create the unitary matrix
        return np.array(self.compute_Y(),dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return f'[Y({self._a,self._b})]'
    
class Ydag(cirq.Gate):
    
    def __init__(self, d: int, a: int, b: int) -> None:
        self._d = d
        self._a = a
        self._b = b
    
    """A gate that enacts the Y^dagger on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)

    def _validate_args(self, qubits):
        return True 

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

    def compute_Y(self):
        X = np.array(self.create_shift_matrix(self._d,self._a),dtype=np.complex64)
        Z = np.array(self.create_roots_of_unity_matrix(self._d,self._b),dtype=np.complex64)
        return X @ Z

    def _unitary_(self):
        # create the unitary matrix
        return np.array(np.conjugate(self.compute_Y()).T,dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return f'[Y*({self._a,self._b})]'
