from project import get_marg_prob, assignCodes, trimTree, buildTree, get_expected_length, encoder
from suppl.utils import *

if __name__ == '__main__':
    f = load_text_sample(name="suppl/genome.txt")
    if len(f) > 0:
        print('Text successfully loaded (starts with {})'.format(f[:10]))

    marg_prob = get_marg_prob(f)
    print(marg_prob)
    print(sum(marg_prob.values(), 0.0))  # verification

    marg_prob_list = list(marg_prob.items())
    Tree = buildTree(marg_prob_list)
    trim_tree = trimTree(Tree)
    codes = {}
    assignCodes(trim_tree, {**codes})
    print('Huffman code for genome : ' + str(assignCodes.codes))

    out = encoder(f, assignCodes.codes)

