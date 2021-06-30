
from difflib import SequenceMatcher


class Base(object):
    A = 'A'
    C = 'C'
    G = 'G'
    T = 'T'
    N = 'N'
    all = [A, C, G, T]
    ix = {A: 0, C: 1, G: 2, T: 3}
    dtypes = {A: int, C: int, G: int, T: int}


def sequence_to_bin(seq):
    """
    Convert string of bases ('A', 'C', 'G', 'T', 'N') to binary format.

    :param seq: string
                sequence of bases ('A', 'C', 'G', 'T', 'N')
    :return: list of 4 strings of the binary presentation of each base.
    for example: the sequence 'AC'
                 the return value: ['2', '1', '0', '0']
                 (A -> 10 -> 2, C -> 01 -> 1, G -> 00 -> 0, T -> 00 -> 0)
    """
    base_dict = {Base.A: '0', Base.G: '0', Base.C: '0', Base.T: '0', Base.N: '0'}
    bin_dict = {}
    for base in Base.all:
        curr_dict = base_dict.copy()
        curr_dict[base] = '1'
        bin_str = ''.join(map(lambda b: curr_dict[b], seq))
        bin_dict[base] = str(int(bin_str, 2))
    return bin_dict


def bin_to_sequence(bin_A, bin_C, bin_G, bin_T, seq_len):
    """
    :param bin_A: binary string representing the 'A's in the sequence
    :param bin_C: binary string representing the 'C's in the sequence
    :param bin_G: binary string representing the 'G's in the sequence
    :param bin_T: binary string representing the 'T's in the sequence
    :param seq_len: the sequence length
    :return: the sequence, length seq_len.
    """
    str_A = "{0:b}".format(int(bin_A)).zfill(seq_len)
    str_C = "{0:b}".format(int(bin_C)).zfill(seq_len)
    str_G = "{0:b}".format(int(bin_G)).zfill(seq_len)
    str_T = "{0:b}".format(int(bin_T)).zfill(seq_len)
    res = []
    for sa, sc, sg, st in zip(str_A, str_C, str_G, str_T):
        if sa == '1':
            res += 'A'
        elif sc == '1':
            res += 'C'
        elif sg == '1':
            res += 'G'
        else:
            res += 'T'
    res = "".join(res)
    return res


def similar(seq_a, seq_b):
    return SequenceMatcher(None, seq_a, seq_b).ratio()