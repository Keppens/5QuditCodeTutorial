# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 17:04:51 2023

Qudit Sum gate

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np

class SUM(cirq.Gate):
    
    """A conditional gate that enacts the transformation SUM|m〉|n〉 = |m〉|n + m mod d〉.
    """
    
    def __init__(self, m, n: int) -> None:
        self._m = m
        self._n = n

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (m,n)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a qubit and a ququart.
        return (self._m,self._n)

    def _validate_args(self, qubits):
        return True 

    def create_block_matrix(self, m, n):
        if m < 1 or n < 1:
            raise ValueError("Both 'm' and 'n' must be at least 1.")
        
        unitary_matrix = np.eye(n,dtype=np.complex128)
        block_matrices = [np.roll(unitary_matrix.astype(np.complex128), i, axis=0) for i in range(m)]
        
        # Initialize the block matrix with zeros
        block_matrix = np.zeros((m*n, m*n),dtype=np.complex128)
        
        # Fill the diagonal with block matrices
        for i in range(m):
            block_matrix[i*n:(i+1)*n, i*n:(i+1)*n] = block_matrices[i]
        
        return block_matrix

    def _unitary_(self):
        return np.array(self.create_block_matrix(self._m,self._n),dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return 'o','[+]'
    
class SUMinv(cirq.Gate):
    
    """A conditional gate that enacts the transformation SUM|m〉|n〉 = |m〉|n + m mod d〉.
    """
    
    def __init__(self, m, n: int) -> None:
        self._m = m
        self._n = n

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (m,n)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a qubit and a ququart.
        return (self._m,self._n)

    def _validate_args(self, qubits):
        return True 

    def create_block_matrix(self, m, n):
        if m < 1 or n < 1:
            raise ValueError("Both 'm' and 'n' must be at least 1.")
        
        unitary_matrix = np.eye(n,dtype=np.complex128)
        block_matrices = [np.roll(unitary_matrix.astype(np.complex128), i, axis=0) for i in range(m)]
        
        # Initialize the block matrix with zeros
        block_matrix = np.zeros((m*n, m*n),dtype=np.complex128)
        
        # Fill the diagonal with block matrices
        for i in range(m):
            block_matrix[i*n:(i+1)*n, i*n:(i+1)*n] = block_matrices[i]
        
        return block_matrix

    def _unitary_(self):
        return np.array(np.linalg.inv(self.create_block_matrix(self._m,self._n)),dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return 'o','[-]'
    
class SUMdag(cirq.Gate):
    
    """A conditional gate that enacts the transformation SUM|m〉|n〉 = |m〉|n + m mod d〉.
    """
    
    def __init__(self, m, n: int) -> None:
        self._m = m
        self._n = n

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (m,n)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a qubit and a ququart.
        return (self._m,self._n)

    def _validate_args(self, qubits):
        return True 

    def create_block_matrix(self, m, n):
        if m < 1 or n < 1:
            raise ValueError("Both 'm' and 'n' must be at least 1.")
        
        unitary_matrix = np.eye(n,dtype=np.complex128)
        block_matrices = [np.roll(unitary_matrix.astype(np.complex128), i, axis=0) for i in range(m)]
        
        # Initialize the block matrix with zeros
        block_matrix = np.zeros((m*n, m*n),dtype=np.complex128)
        
        # Fill the diagonal with block matrices
        for i in range(m):
            block_matrix[i*n:(i+1)*n, i*n:(i+1)*n] = block_matrices[i]
        
        return block_matrix

    def _unitary_(self):
        return np.array(np.conjugate(self.create_block_matrix(self._m,self._n)).T,dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return 'o','[+*]'
    
class MIN(cirq.Gate):
    
    """A conditional gate that enacts the transformation SUM|m〉|n〉 = |m〉|n + m mod d〉.
    """
    
    def __init__(self, m, n: int) -> None:
        self._m = m
        self._n = n

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (m,n)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a qubit and a ququart.
        return (self._m,self._n)

    def _validate_args(self, qubits):
        return True 

    def create_block_matrix(self, m, n):
        if m < 1 or n < 1:
            raise ValueError("Both 'm' and 'n' must be at least 1.")
        
        unitary_matrix = np.eye(n,dtype=np.complex128)
        block_matrices = [np.roll(unitary_matrix.astype(np.complex128), -i, axis=0) for i in range(m)]
        
        # Initialize the block matrix with zeros
        block_matrix = np.zeros((m*n, m*n),dtype=np.complex128)
        
        # Fill the diagonal with block matrices
        for i in range(m):
            block_matrix[i*n:(i+1)*n, i*n:(i+1)*n] = block_matrices[i]
        
        return block_matrix

    def _unitary_(self):
        return np.array(self.create_block_matrix(self._m,self._n),dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return 'o','[-]'

class CShift(cirq.Gate):
    
    """A conditional gate that enacts the transformation SUM|m〉|n〉 = |m〉|n + 1 mod d〉 if m is not 0.
    """
    
    def __init__(self, m, n: int) -> None:
        self._m = m
        self._n = n

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (m,n)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a qubit and a ququart.
        return (self._m,self._n)

    def _validate_args(self, qubits):
        return True 

    def create_block_matrix_adapted(self, m, n):
        if m < 1 or n < 1:
            raise ValueError("Both 'm' and 'n' must be at least 1.")
        
        unitary_matrix = np.eye(n, dtype=np.complex128)
        first_block_matrix = np.roll(unitary_matrix.astype(np.complex128), 0, axis=0)  # No row roll for the first block
        rest_block_matrices = [np.roll(unitary_matrix.astype(np.complex128), 1, axis=0) for _ in range(m - 1)]
    
        # Initialize the block matrix with zeros
        block_matrix = np.zeros((m * n, m * n), dtype=np.complex128)
    
        # Fill the diagonal with block matrices
        block_matrix[:n, :n] = first_block_matrix
        for i in range(1, m):
            block_matrix[i * n:(i + 1) * n, i * n:(i + 1) * n] = rest_block_matrices[i - 1]
    
        return block_matrix


    def _unitary_(self):
        return np.array(self.create_block_matrix_adapted(self._m,self._n),dtype=np.complex128)
    
    def _circuit_diagram_info_(self, args):
        return 'o','[+]'
