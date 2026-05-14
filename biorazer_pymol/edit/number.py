from pymol import cmd


def renumber_residue_continuously(selection: str, start_resi: int | None = None):
    """
    DESCRIPTION

    Renumber residues in a selection so that residue IDs become continuous.

    The input selection may contain discontinuous residue numbering, such as
    1-10 and 15-20, but it must come from a single chain.

    USAGE

    renumber_residue_continuously selection[, start_resi]

    ARGUMENTS

    selection = str
        PyMOL selection containing residues to renumber.

    start_resi = int or None
        Starting residue number for renumbering. If omitted, numbering starts
        from the smallest residue number present in the selection.

    EXAMPLE

    renumber_residue_continuously chain A and resi 1-10+15-20
    renumber_residue_continuously polymer.protein and chain B, 101

    NOTES

    Only selections from a single chain are supported.
    """

    chain_list = []

    def check_chain(atom):
        if not atom.chain in chain_list:
            chain_list.append(atom.chain)

    cmd.iterate(selection=selection, expression=check_chain)
    if len(chain_list) > 1:
        raise ValueError("Selection does not come from 1 chain")

    resi_list = []

    def collect_resi(atom):
        resi_list.append((atom.index, int(atom.resi)))

    cmd.iterate(selection=selection, expression=collect_resi)

    renumbered_resi_list = []
    mapping = {}
    if start_resi is None:
        next_id = min(list(map(lambda x: x[1], resi_list)))
    else:
        next_id = start_resi
    for index, resi in resi_list:
        if resi not in mapping:
            mapping[resi] = next_id
            next_id += 1
        renumbered_resi_list.append((index, mapping[resi]))

    for index, resi in renumbered_resi_list:

        cmd.alter(
            selection=f"({selection}) and index {index}",
            expression=f"resi=str({resi})",
        )


cmd.extend("renumber_residue_continuously", renumber_residue_continuously)
