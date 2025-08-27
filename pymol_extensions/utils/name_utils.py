from pymol import cmd

if cmd.get_version()[1] < 1.2:
    def get_unused_name(name):
        import random
        return name + '%04d' % random.randint(0, 1000)
    STATE = 1

else:
    from pymol.cmd import get_unused_name
    STATE = -1