from pymol import cmd

def calc_hydrophobic_sasa(selection='all', state=1, quiet=0):
    """
    DESCRIPTION

    Calculate the hydrophobic solvent accessible surface area (SASA) of a
    selection.

    USAGE

    calc_hydro_sasa [selection [, state [, quiet]]]

    ARGUMENTS

    selection = string: a selection-expression or object-name {default: all}

    state = integer: state number {default: 1}

    quiet = 0/1: suppress output {default: 1}

    NOTES

    The hydrophobic SASA is calculated using the Lee-Richards algorithm with
    a probe radius of 1.4 Angstroms.

    EXAMPLES

    calc_hydro_sasa

    calc_hydro_sasa (chain A)

    calc_hydro_sasa (resi 123)

    SEE ALSO

    calc_sasa
    """
    dot_solvent = cmd.get('dot_solvent')
    cmd.set('dot_solvent', 1)
    
    # Calculate the SASA of the selection
    sasa = cmd.get_area(selection, load_b=1, state=state, quiet=1)

    # Calculate the hydrophobic SASA
    hydro_sasa = cmd.get_area(selection + ' and (elem C or elem H)', load_b=1,
                              state=state, quiet=1)

    # Calculate the hydrophobic SASA fraction
    hydro_sasa_fraction = hydro_sasa / sasa

    # Print the results
    if not quiet:
        print('SASA: %.2f' % sasa)
        print('Hydrophobic SASA: %.2f' % hydro_sasa)
        print('Hydrophobic SASA fraction: %.2f' % hydro_sasa_fraction)
    
    cmd.set('dot_solvent', dot_solvent)

    return hydro_sasa_fraction

cmd.extend('calc_hydrophobic_sasa', calc_hydrophobic_sasa)