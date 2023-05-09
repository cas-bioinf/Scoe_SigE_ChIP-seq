import tempfile
import pandas as pd


def parse_motif_sites(f):
    """This is coarse motif sites parser
    use with caution - check if meme parsing in Biopython is finally supported before using this
    """
    r_start = "sites sorted by position p-value"
    for line in f:
        if r_start in line:
            motif_name = line.strip()[:-(len(r_start)+1)].split(' ')[-1]

            for _ in range(2):
                t = next(f)

            pointer_line = next(f)

            line = next(f)
            tab = []
            while not char_all(line.strip(), '-'):
                tab.append(line)
                line = next(f)

            yield motif_name, tab, pointer_line


def char_all(l, char):
    """
    check if line consists of chars only
    :param l: str
    :return: bool
    """
    ls = set(l)
    if len(ls) == 1 and list(ls)[0] == char:
        return True
    else:
        return False


def read_file(file, no_strand=False):
    columns = ['motif_id', 'sequence_name', 'strand', 'start', 'P-value', '_site_left', 'site', '_site_right']
    if no_strand:
        columns.remove('strand')
    with tempfile.TemporaryFile(mode="r+") as t, open(file, 'r') as f:
        for motif_name, table, pointer in parse_motif_sites(f):
            for line in table:
                t.write("{} {}".format(motif_name, line))
        
        t.seek(0)
        return pd.read_table(t, names=columns, header=None, delim_whitespace=True)
        
