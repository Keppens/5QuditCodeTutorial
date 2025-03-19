# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 16:41:58 2023

Qudit Shift gate

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np

class Shift(cirq.Gate):
    
    def __init__(self, d: int, a: int) -> None:
        self._d = d
        self._a = a
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
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

    def _unitary_(self):
        # create the unitary matrix
        return np.array(self.create_shift_matrix(self._d,self._a),dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return f'[X({self._a})]'
    
class Shiftdag(cirq.Gate):
    
    def __init__(self, d: int, a: int) -> None:
        self._d = d
        self._a = a
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
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

    def _unitary_(self):
        # create the unitary matrix
        return np.array(np.conjugate(self.create_shift_matrix(self._d,self._a)).T,dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return f'[X*({self._a})]'
