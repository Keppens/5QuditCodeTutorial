# -*- coding: utf-8 -*-
"""
Created on Fri Dec 15 17:37:30 2023

Qudit X2 gate

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np

class M(cirq.Gate):
    
    
    def __init__(self, d: int, g: int) -> None:
        self._d = d
        self._g = g
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)

    def _validate_args(self, qubits):
        return True 

    def create_qudit_multiplication_gate(self, d, g):
        # Check if g is within the valid range
        if g < 1 or g >= d:
            raise ValueError("Invalid value for 'g'. It must be in the range (1, d-1).")
    
        # Create the qudit multiplication gate matrix
        multiplication_gate_matrix = np.zeros((d, d), dtype=np.complex128)
    
        for i in range(d):
            multiplication_gate_matrix[(i * g) % d, i] = 1
    
        return multiplication_gate_matrix

    def _unitary_(self):
        # create the unitary matrix
        return np.array(self.create_qudit_multiplication_gate(self._d,self._g),dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return f'[x{self._g}]'

class Minv(cirq.Gate):
    
    
    def __init__(self, d: int, g: int) -> None:
        self._d = d
        self._g = g
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)

    def _validate_args(self, qubits):
        return True 

    def create_qudit_multiplication_gate(self, d, g):
        # Check if g is within the valid range
        if g < 1 or g >= d:
            raise ValueError("Invalid value for 'g'. It must be in the range (1, d-1).")
    
        # Create the qudit multiplication gate matrix
        multiplication_gate_matrix = np.zeros((d, d), dtype=np.complex128)
    
        for i in range(d):
            multiplication_gate_matrix[(i * g) % d, i] = 1
    
        return multiplication_gate_matrix

    def _unitary_(self):
        # create the unitary matrix
        return np.array(np.linalg.inv(self.create_qudit_multiplication_gate(self._d,self._g)),dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return f'[x{self._g}-]'

class Mdag(cirq.Gate):
    
    
    def __init__(self, d: int, g: int) -> None:
        self._d = d
        self._g = g
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)

    def _validate_args(self, qubits):
        return True 

    def create_qudit_multiplication_gate(self, d, g):
        # Check if g is within the valid range
        if g < 1 or g >= d:
            raise ValueError("Invalid value for 'g'. It must be in the range (1, d-1).")
    
        # Create the qudit multiplication gate matrix
        multiplication_gate_matrix = np.zeros((d, d), dtype=np.complex128)
    
        for i in range(d):
            multiplication_gate_matrix[(i * g) % d, i] = 1
    
        return multiplication_gate_matrix

    def _unitary_(self):
        # create the unitary matrix
        return np.array(np.conjugate(self.create_qudit_multiplication_gate(self._d,self._g)).T,dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return f'[x{self._g}-]'