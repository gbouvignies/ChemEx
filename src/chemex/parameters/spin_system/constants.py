# Conversion dictionary from 3-letter to 1-letter amino-acid convention
AAA_TO_A = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


# Conversion dictionary to correct different ways to spell nuclei
CORRECT_ATOM_NAME = {"HN": "H", "C'": "C", "CO": "C"}

# fmt: off
STANDARD_ATOM_NAMES_PROTEIN = {
    "C", "CA", "CB", "CD", "CD1", "CD2", "CE", "CE1", "CE2", "CE3", "CG", "CG1", "CG2",
    "CH2", "CQD", "CQE", "CQG", "CZ", "CZ2", "CZ3", "H", "H2", "H3", "HA", "HA2",
    "HA3", "HB", "HB1", "HB2", "HB3", "HD", "HD1", "HD11", "HD12", "HD13", "HD2",
    "HD21", "HD22", "HD23", "HD3", "HE", "HE1", "HE2", "HE21", "HE22", "HE3", "HG",
    "HG1", "HG11", "HG12", "HG13", "HG2", "HG21", "HG22", "HG23", "HG3", "HH", "HH1",
    "HH11", "HH12", "HH2", "HH21", "HH22", "HZ", "HZ1", "HZ2", "HZ3", "MB", "MD", "MD1",
    "MD2", "ME", "MG", "MG1", "MG2", "MZ", "N", "ND1", "ND2", "NE", "NE1", "NE2", "NH",
    "NH1", "NH2", "NQH", "NZ", "QA", "QB", "QD", "QD2", "QE", "QE2", "QG", "QG1", "QH1",
    "QH2", "QMD", "QMG", "QQH", "QR", "QZ",
}

STANDARD_ATOM_NAMES_DNA_RNA = {
    "C1'", "C2", "C2'", "C3'", "C4", "C4'", "C5", "C5'", "C6", "C7", "C8",
    "H1", "H1'", "H2", "H2'", "H2'1", "H2'2", "H21", "H22", "H3", "H3'", "H4'", "H41",
    "H42", "H5", "H5'1", "H5'2", "H6", "H61", "H62", "H71", "H72", "H73", "H8",
    "M7", "N1", "N2", "N3", "N4", "N6", "N7", "N9",
}
# fmt: on

STANDARD_ATOM_NAMES = STANDARD_ATOM_NAMES_PROTEIN | STANDARD_ATOM_NAMES_DNA_RNA
