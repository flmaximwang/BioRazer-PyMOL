from pymol import cmd


def _normalize_fmt(fmt):
    fmt = fmt.lower()
    if fmt in ["fasta", "fa"]:
        return "fasta"
    else:
        raise ValueError(f"Unsupported format: {fmt}")


def print_seq(selection="all", format="fasta", **kwargs):
    format = _normalize_fmt(format)
    if format == "fasta":
        # 对每个 Chain 输出一个 Fasta 条目
        for chain in cmd.get_chains(selection):
            seq = cmd.get_fastastr(f"{selection} and chain {chain}")
            print(f">{chain}\n{seq}")
    else:
        raise ValueError(f"Unsupported format: {format}")


cmd.extend("print_seq", print_seq)
