from project import *
from suppl.utils import *

if __name__ == '__main__':
    f = load_text_sample(name="suppl/genome.txt")
    if len(f) > 0:
        print('Text successfully loaded (starts with {})'.format(f[:10]))

    marg_prob = get_marg_prob(f)
