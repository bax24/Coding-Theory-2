import numpy as np
from suppl.utils import *
import collections


def buildTree(Pd):
    Pd.sort(key=lambda t: t[1])
    while len(Pd) > 1:
        leastTwo = tuple(Pd[0:2])  # get the 2 to combine
        theRest = Pd[2:]  # all the others
        combFreq = leastTwo[0][1] + leastTwo[1][1]  # the branch points freq
        Pd = theRest + [(leastTwo, combFreq)]  # add branch point to the end
        Pd.sort(key=lambda t: t[1])

    return Pd[0]  # sort it into place


def assignCodes(node, codes, pat=''):
    assignCodes.codes = codes
    if type(node) == type(""):
        codes[node] = pat  # A leaf. set its code
    else:  #
        assignCodes(node[0], codes, pat + "0")  # Branch point. Do the left branch
        assignCodes(node[1], codes, pat + "1")


def trimTree(tree):
    # Trim the freq counters off, leaving just the letters
    p = tree[0]  # ignore freq count in [0]
    if type(p) == type(""):
        return p  # if just a leaf, return it
    else:
        return trimTree(p[0]), trimTree(p[1])  # trim left then right and recombine


def LZ77_search(search, look_ahead):
    ls = len(search)
    llh = len(look_ahead)

    if ls == 0:
        return 0, 0, look_ahead[0]

    if llh == 0:
        return -1, -1, ""

    best_length = 0
    best_offset = 0
    buf = search + look_ahead

    search_pointer = ls
    redundancy = False
    for i in range(0, ls):
        length = 0
        while buf[i + length] == buf[search_pointer + length]:
            redundancy = True
            length = length + 1
            if search_pointer + length == len(buf):
                length = length - 1
                break
            if i + length >= search_pointer:
                break
        if length >= best_length:
            if redundancy:
                best_offset = ls - i
                best_length = length
            else:
                best_offset = 0

    return best_offset, best_length, buf[search_pointer + best_length]


def LZ77(input_string, size_l):
    output = []
    window = ''
    i = 0
    while input_string:
        out = LZ77_search(window, input_string)
        output.append(out)
        length = out[1]
        i = i + length + 1
        while length >= 0:
            if len(window) >= size_l:
                window = window[1:]
            window = window + input_string[0]
            input_string = input_string[1:]
            length = length - 1
    return output


def get_expected_length(huff_codes, prob):
    codes_lengths = 0
    for key in huff_codes:
        codes_lengths += prob[key]*len(huff_codes[key])

    return codes_lengths


def get_marg_prob(f):
    prob = {}
    total = len(f) / 3
    while f:
        if f[0:3] in prob:
            prob[f[0:3]] += 1
        else:
            prob[f[0:3]] = 1
        f = f[3:]

    prob = {k: v / total for k, v in prob.items()}
    return prob


if __name__ == "__main__":

    Huffman, Lempel_zev, LZ77_algo = False, False, True
    Q5, Q6, Q7 = False, False, False

    if Huffman:
        # Exercise 7 verification
        Pd = {'0': 0.05,
              '1': 0.10,
              '2': 0.15,
              '3': 0.15,
              '4': 0.20,
              '5': 0.35}

        Pd_list = list(Pd.items())
        Tree = buildTree(Pd_list)
        trim_tree = trimTree(Tree)
        codes = {}
        assignCodes(trim_tree, {**codes})
        print('Huffman codes is : ' + str(assignCodes.codes))

    if LZ77_algo:
        # Example verification
        string = 'abracadabrad'
        sliding_size = 7
        _output = LZ77(string, sliding_size)
        print('LZ77 code is : ' + str(_output))

    if Q5:
        f = load_text_sample(name="suppl/genome.txt")
        if len(f) > 0:
            print('Text successfully loaded (starts with {})'.format(f[:10]))

        marg_prob = get_marg_prob(f)
        print(marg_prob)
        print(sum(marg_prob.values(), 0.0)) # verification

        marg_prob_list = list(marg_prob.items())
        Tree = buildTree(marg_prob_list)
        trim_tree = trimTree(Tree)
        codes = {}
        assignCodes(trim_tree, {**codes})
        print('Huffman code for genome : ' + str(assignCodes.codes))

    if Q6:
        out = get_expected_length(assignCodes.codes, marg_prob)
        print('expected average length is : ' + str(out))

