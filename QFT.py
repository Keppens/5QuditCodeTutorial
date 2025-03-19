# -*- coding: utf-8 -*-
"""
Created on Wed Dec  6 14:48:31 2023

Qudit QFT gate

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np

class H(cirq.Gate):
    def __init__(self, d: int) -> None:
        self._d = d
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)
    
    def _validate_args(self, qubits):
        return True 

    def create_qft_matrix(self, d):
        # Initialize the QFT matrix
        qft_matrix = np.zeros((d, d), dtype=np.complex128)
    
        # Fill in the QFT matrix
        for a in range(d):
            for b in range(d):
                qft_matrix[a, b] = np.exp(2j * np.pi * a * b / d)
    
        # Normalize the matrix
        qft_matrix /= np.sqrt(d)
    
        return qft_matrix

    def _unitary_(self):
        # create the unitary matrix
        return np.array(self.create_qft_matrix(self._d),dtype=np.complex64)
    def _circuit_diagram_info_(self, args):
        return '[F]'
    
class Hinv(cirq.Gate):
    def __init__(self, d: int) -> None:
        self._d = d
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)
    
    def _validate_args(self, qubits):
        return True 

    def create_qft_matrix(self, d):
        # Initialize the QFT matrix
        qft_matrix = np.zeros((d, d), dtype=np.complex128)
    
        # Fill in the QFT matrix
        for a in range(d):
            for b in range(d):
                qft_matrix[a, b] = np.exp(2j * np.pi * a * b / d)

        # Normalize the matrix
        qft_matrix /= np.sqrt(d)
    
        return qft_matrix

    def _unitary_(self):
        # create the unitary matrix
        return np.array(np.linalg.inv(self.create_qft_matrix(self._d)),dtype=np.complex64)
    def _circuit_diagram_info_(self, args):
        return '[F-]'


class Hdag(cirq.Gate):
    def __init__(self, d: int) -> None:
        self._d = d
    
    """A gate that enacts the transformation U|x〉 = |x + 1 mod d〉 on a qudit.
    """

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        return (self._d,)
    
    def _validate_args(self, qubits):
        return True 

    def create_qft_matrix(self, d):
        # Initialize the QFT matrix
        qft_matrix = np.zeros((d, d), dtype=np.complex128)
    
        # Fill in the QFT matrix
        for a in range(d):
            for b in range(d):
                qft_matrix[a, b] = np.exp(2j * np.pi * a * b / d)
    
        # Normalize the matrix
        qft_matrix /= np.sqrt(d)
    
        return qft_matrix

    def _unitary_(self):
        # create the unitary matrix
        return np.array(np.conjugate(self.create_qft_matrix(self._d)).T,dtype=np.complex64)

    def _circuit_diagram_info_(self, args):
        return '[F*]'