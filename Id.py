# -*- coding: utf-8 -*-
"""
Created on Tue Aug  6 10:45:14 2024

@author: James Keppens
Based on https://quantumai.google/cirq/build/qudits
"""
#Imports
import cirq
import numpy as np

class I(cirq.Gate):
    
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

    def _unitary_(self):
        # create the unitary matrix
        return np.array(np.eye(self._d))

    def _circuit_diagram_info_(self, args):
        return f'[I]'