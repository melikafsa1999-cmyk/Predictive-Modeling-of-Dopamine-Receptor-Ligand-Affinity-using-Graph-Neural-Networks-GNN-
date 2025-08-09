
"""
Data Featurization Functions
"""

import numpy as np
from rdkit import Chem

# ------------------------------
# 1. Data Featurization Functions
# ------------------------------
def featurize_atom(atom):
    """Extract feature vector for one atom."""
    features = []
    # (1) Atomic number
    features.append(atom.GetAtomicNum())
    # (2) Degree
    features.append(atom.GetDegree())
    # (3) Formal charge
    features.append(atom.GetFormalCharge())
    features.append(atom.GetTotalValence())
    features.append(atom.GetNumRadicalElectrons())
    # (4) Hybridization (one-hot encoding: SP, SP2, SP3, OTHER)
    hybrid_map = {
        Chem.rdchem.HybridizationType.SP: [1, 0, 0, 0],
        Chem.rdchem.HybridizationType.SP2: [0, 1, 0, 0],
        Chem.rdchem.HybridizationType.SP3: [0, 0, 1, 0]
    }
    hybridization = atom.GetHybridization()
    if hybridization in hybrid_map:
        features.extend(hybrid_map[hybridization])
    else:
        features.extend([0, 0, 0, 1])
    # (5) Is aromatic
    features.append(1 if atom.GetIsAromatic() else 0)
    return np.array(features, dtype=float)


def featurize_bond(bond):
    """Extract feature vector for one bond."""
    features = []
    # (1) Bond type (one-hot: SINGLE, DOUBLE, TRIPLE, AROMATIC)
    bond_map = {
        Chem.rdchem.BondType.SINGLE: [1, 0, 0, 0],
        Chem.rdchem.BondType.DOUBLE: [0, 1, 0, 0],
        Chem.rdchem.BondType.TRIPLE: [0, 0, 1, 0],
        Chem.rdchem.BondType.AROMATIC: [0, 0, 0, 1]
    }
    bond_type = bond.GetBondType()
    if bond_type in bond_map:
        features.extend(bond_map[bond_type])
    else:
        features.extend([0, 0, 0, 0])
    # (2) Conjugation
    features.append(1 if bond.GetIsConjugated() else 0)
    # (3) Ring membership
    features.append(1 if bond.IsInRing() else 0)
    return np.array(features, dtype=float)
