"""
test for the SCF module
"""


import HF_SCF as hs
import pytest
import psi4


def test_scf():
    
    mol=psi4.geometry('''
        O
        H 1 1.1
        H 1 1.1 2 104''')
    mol.update_geometry()

     
    psi4.set_options({"scf_type":"pk"})
    psi4_energy = psi4.energy("SCF/sto-3g",molecule=mol)
    test_energy=hs.get_energy("sto-3g")
    print(test_energy)
    print(psi4_energy)
    assert  abs(test_energy-psi4_energy)<0.0001 
 

