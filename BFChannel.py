# -*- coding: utf-8 -*-
"""
Created on Tue Sep 12 17:20:51 2023

Qudit bitflip channel

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np

class BFd(cirq.Gate):
    
    """A bit flip channel for qudits
    """
    
    def __init__(self,p,d):
        self._p = p
        self._d = d

    def _qid_shape_(self):
        # By implementing this method this gate implements the
        # cirq.qid_shape protocol and will return the tuple (d,)
        # when cirq.qid_shape acts on an instance of this class.
        # This indicates that the gate acts on a qubit and a ququart.
        return (self._d,)
    
    def create_shift_matrix(self,d):
        if d < 2:
            raise ValueError("Dimension 'd' must be at least 2.")
        
        # Perform the shift operation (x -> x+1 mod d) on the matrix
        shift_matrix = np.roll(np.eye(d), -1, axis=0)
        
        return shift_matrix

    def _mixture_(self):
        ps = [1.0 - self._p, self._p]
        ops = [np.array(np.eye(self._d),dtype=np.complex128),np.array(self.create_shift_matrix(self._d),dtype=np.complex64)]
        return tuple(zip(ps, ops))

    def _has_mixture_(self) -> bool:
        return True

    def _circuit_diagram_info_(self, args) -> str:
        return f"BFd({self._p})"
