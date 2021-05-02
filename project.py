import numpy as np
from suppl.utils import *
import math
import matplotlib.pyplot as plt
from collections import OrderedDict
from functools import reduce
import operator as op


def buildTree(Pd):
    Pd.sort(key=lambda t: t[1])
    while len(Pd) > 1:
        leastTwo = tuple(Pd[0:2])  # get the 2 to combine
        theRest = Pd[2:]  # all the others
        combFreq = leastTwo[0][1] + leastTwo[1][1]  # the branch points freq
        Pd = theRest + [(leastTwo, combFreq)]  # add branch point to the end
        Pd.sort(key=lambda t: t[1])

    return Pd[0]  # sort it into place


# Hello adri

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
        codes_lengths += prob[key] * len(huff_codes[key])

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


def get_length_genome():
    file = load_binary_text_sample()
    if len(file) > 0:
        print('Binary text successfully loaded (starts with {})'.format(f[:10]))

    return len(file)


def encoder(message, code):
    encoding = ''
    while message:
        element = message[0:3]
        encoding += code[element]
        message = message[3:]
    return encoding


def ACGT_to_binary():
    
    chain_length = len(f)
    ACGT_bin = ""
    
    for i in range(0, chain_length):
        
        if f[i] == "A" :
            ACGT_bin += "01000001"
        elif f[i] == "C" :
            ACGT_bin += "01000011"
        elif f[i] == "G" :
            ACGT_bin += "01000111"
        elif f[i] == "T" :
            ACGT_bin += "01010100"
        else :
            print("Error : Unknown symbol \'",f[i],"\'.")
            
    return ACGT_bin
    

def binarize(v, nb=0):
    if nb == 0:
        return ''
    else:
        return np.binary_repr(v,width=nb)
    
    
def lempel_ziv_encoder(online):
    
# =============================================================================
# online can take 2 values :
#     0 for a fixed sized memory address
#     1 for a memory address depending on the size of the current dictionary
# =============================================================================

# We consider here an ACTG chain, exprimed in a ASCII (8-bits) representation
    
    dictio = OrderedDict()
    chain_length = len(ACGT_bin)
    cc = "" # for Current Character(s)
    nb = 0  # For Number of Bit(s)
    code = ""
    idx = int(1)
    dictio[''] = ''
    
    if not online :
        # cod, dic = lempel_ziv_encoder(True)
        # print("code_length : ",str(len(cod))," | ratio : ",str(len(ACGT_bin)/len(cod)))
        # int(np.ceil(np.log2(len(dic))))
        nb = 3
    
    for i in range(0,chain_length):
        cc += ACGT_bin[i]
        
        if cc in dictio :
            continue 
        
        elif i == chain_length :
            code += list(dictio.keys()).index(cc)
        
        else :
            
            if online :
                nb = int(np.ceil(np.log2(idx)))
                
            new_val = binarize(list(dictio.keys()).index(cc[:-1]),nb)
            dictio[cc] = new_val + cc[-1]
            code += new_val + cc[-1]
            idx += 1 ; cc = ""
        
    return code, dictio

# This function loads the sound signal (.wav)
def load_wav():
	rate, data = read('suppl//sound.wav')
	return rate, data

# This function save the sound signal (.wav)
def save_wav(filename,rate,data):
    write(filename, rate, data)
    
    
def sound_to_binary():
    
    bin_sound = ""
    sound_length = len(s)
    
    for i in range(0,sound_length):
        bin_sound += '{0:08b}'.format(s[i])
        
    return bin_sound


def simulate_channel(sound_in):
    
    signal_out = ""
    signal_length = len(sound_in)
    cnt = 0 
    
    for i in range(0,signal_length):
        
        if np.random.randint(1,100) == 1 :
            cnt += 1
            if sound_in[i] == "0" :
                signal_out += "1" 
            else :
                signal_out += "0"
                
        else :
            signal_out += sound_in[i]
            
    # print("Percentage of randomly flipped bits : ",cnt/signal_length)
            
    return signal_out


def binary_to_sound(sound_in):
    
    length_sound = len(sound_in)
    length_out = int(length_sound/8)
    sound_out = np.zeros(length_out, np.uint8)
    
    for i in range(0,length_sound, 8):
        sound_out[int(i/8)] = int(sound_in[i:i+8],2)
        
    return sound_out


def hamming_7_4_encoder(extend=True):

    length_sound = len(bin_sound)
    hamming_out = ""
    
    for i in range(0,length_sound,4):
        
        code = np.zeros(8,int)
        code[3], code[5], code[6], code[7] = bin_sound[i:i+4]
        if sum(code) == 0 :
            idx = int(0)
        else :
            idx = reduce(op.xor, [j for j, bit in enumerate(code) if bit])
        bits = binarize(idx,3)
        code[4], code[2], code[1] = bits
         
        if extend :
            code[0] = sum(code[0:8])%2
            
        for j in range(0,8):
            hamming_out += str(code[j])     
            
    return hamming_out


def hamming_7_4_decoder(ham_sound):
    
    length_sound = len(ham_sound)
    hamming_in = "" ; curr_bits = np.zeros(8,int)
    flip_bit = "0000" ; is_part_1 = True
    seq = "" ; pass_loop = False
    
    for i in range(0,length_sound,8):
        
        if pass_loop :
            pass_loop = False
            continue 
        
        for j in range(0,8):
            curr_bits[j] = int(ham_sound[i+j])
            
        bits_sum = sum(curr_bits)
            
        if bits_sum == 0 :
            idx = int(0)
        else :
            idx = reduce(op.xor, [k for k, bit in enumerate(curr_bits) if bit])
        bits = binarize(idx,3)
        
        if bits_sum%2 != 0 : # Probably one error
            
            if bits == '011' :      # bit 3 has been flipped
                flip_bit = "1000"
                    
            elif bits == '101' :    # bit 5 has been flipped
                flip_bit = "0100"
                
            elif bits == '110' :    # bit 6 has been flipped
                flip_bit = "0010"
                
            elif bits == '111' :    # bit 7 has been flipped
                flip_bit = "0001"
                
            else :                  # a parity bit has been flipped
                flip_bit = "0000"
                
        elif bits_sum%2 == 0 and bits != '000' : # 2 erros
        
            if i >= 8 :
                if is_part_1 :
                    hamming_in += hamming_in[int(i/2)-8:int(i/2)]
                    pass_loop = True
                    
                else :
                    hamming_in += hamming_in[int(i/2)-8:int(i/2)-4]
                    is_part_1 = True
                            
                continue 
                
        else :
            flip_bit = "0000"
                
        seq = ham_sound[i+3] + ham_sound[i+5] + ham_sound[i+6] + ham_sound[i+7]          
        hamming_in += str(binarize(int(seq,2) ^ int(flip_bit,2),4))
        is_part_1 = not is_part_1
            
    return hamming_in
        

def get_error(sound1,sound2,diff=False):
    
    lg1, lg2 = len(sound1), len(sound2)
    
    if lg1 != lg2 :
        print("Error : sounds are of different lengths !")
        return int(-1)
        
    cnt = int(0) ; diff_cnt = 0
        
    for i in range(0,lg1):
        if sound1[i] != sound2[i]:
            if diff :
                diff_cnt += abs(sound1[i] - sound2[i])
                print(diff_cnt)
                
            cnt += 1
           
    print("cnt " ,cnt)
    return cnt/lg1, diff_cnt/cnt


if __name__ == "__main__":

    Huffman, Lempel_ziv, LZ77_algo = False, False, False
    Q5, Q6, Q7 = False, False, False
    Q10 = False
    
    Q15, Q16, Q17, Q18, Q19 = True, True, True, True, True

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
        
    if Lempel_ziv :
        f = load_text_sample()
        if len(f) > 0 :
            print('Text successfully loaded (starts with {})'.format(f[:10]))
        
        ACGT_bin = "1011010100010" # ACGT_to_binary()
        code, final_dictio = lempel_ziv_encoder(online=False)
        
        print("Length of the initial chain : ",str(len(ACGT_bin)))
        print("Length of the compressed chain : ",str(len(code)))
        print("Compression ration = ",str(len(ACGT_bin)/len(code)))

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
        print(sum(marg_prob.values(), 0.0))  # verification

        marg_prob_list = list(marg_prob.items())
        Tree = buildTree(marg_prob_list)
        trim_tree = trimTree(Tree)
        codes = {}
        assignCodes(trim_tree, {**codes})
        print('Huffman code for genome : ' + str(assignCodes.codes))
        encoded = encoder(f, assignCodes.codes)
        # length_genome = get_length_genome()
        length_genome = 8627012
        compression_rate = length_genome / len(encoded)
        print('The length of the encoded genome is : ' + str(len(encoded)))
        print('The compression rate is : ' + str(compression_rate))

    if Q6:
        out = get_expected_length(assignCodes.codes, marg_prob)
        print('expected average length is : ' + str(out))

    if Q10:
        f = load_text_sample(name="suppl/genome.txt")
        if len(f) > 0:
            print('Text successfully loaded (starts with {})'.format(f[:10]))

        sliding_size = 7
        _output = LZ77(f, sliding_size)
        print('LZ77 code is : ' + str(_output))
        length_genome = 8627012
        compression_rate = length_genome / len(_output)
        print('The length of the encoded genome is : ' + str(len(_output)))
        print('The compression rate is : ' + str(compression_rate))
        
        
    if Q15 :
        r,s = load_wav()
        plt.plot(s)
        plt.title("Sound.wav signal")
        
        
    if Q16 :
        bin_sound = sound_to_binary()
        
    if Q17 :
        sound_through_channel = simulate_channel(bin_sound)
        signal_out = binary_to_sound(sound_through_channel)
        error, _ = get_error(s,signal_out)
        print("Error : ",error)
        plt.figure()
        plt.plot(signal_out)
        plt.title("Noisy_sound.wav signal")
        save_wav("suppl//Noisy_sound.wav",r,signal_out)
        
    if Q18 :
        ham_encod_sound = hamming_7_4_encoder(extend=True)
        
    if Q19 :
        ham_through_chan = simulate_channel(ham_encod_sound)
        ham_decod = hamming_7_4_decoder(ham_through_chan)
        ham_sound_out = binary_to_sound(ham_decod)
        error, average_err = get_error(s,ham_sound_out,diff=True)
        print("Hamming error : ",error)
        plt.figure()
        plt.plot(ham_sound_out)
        plt.title("Enhance Hamming_sound.wav signal")
        save_wav("suppl//Hamming_sound.wav",r,ham_sound_out)
        

