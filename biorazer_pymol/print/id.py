from pymol import cmd


def print_res_id(selection="all"):
    model = cmd.get_model(selection)
    res_ids = set()
    for atom in model.atom:
        res_ids.add(int(atom.resi))
    print(sorted(res_ids))


cmd.extend("print_res_id", print_res_id)
