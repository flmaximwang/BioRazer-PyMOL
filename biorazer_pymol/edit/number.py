from pymol import cmd


def renumber_residue_continuously(selection: str):
    """
    Renumber the selection continuously.
    - Selection should come from a single chain
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
    next_id = min(list(map(lambda x: x[1], resi_list)))
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
