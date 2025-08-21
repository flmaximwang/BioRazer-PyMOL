HYDROPHOBIC_AAS = ['ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PHE', 'TRP', 'PRO', 'GLY']

THREE_TO_ONE_MAP = {
    "ALA": 'A',
    "ARG": 'R',
    "ASN": 'N',
    "ASP": 'D',
    "CYS": 'C',
    "GLN": 'Q',
    "GLU": 'E',
    "GLY": 'G',
    "HIS": 'H',
    "ILE": 'I',
    "LEU": 'L',
    "LYS": 'K',
    "MET": 'M',
    "PHE": 'F',
    "PRO": 'P',
    "SER": 'S',
    "THR": 'T',
    "TRP": 'W',
    "TYR": 'Y',
    "VAL": 'V'
}

ONE_TO_THREE_MAP = {v: k for k, v in THREE_TO_ONE_MAP.items()}