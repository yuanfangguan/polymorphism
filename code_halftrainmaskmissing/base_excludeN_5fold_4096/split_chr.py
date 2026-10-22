import random
import sys

chromosomes = ['NC_000001.11', 'NC_000003.12', 'NC_000005.10', 'NC_000007.14', 'NC_000009.12', 'NC_000011.10',
               'NC_000013.11', 'NC_000015.10', 'NC_000017.11', 'NC_000019.10', 'NC_000021.9', 'NC_000002.12',
               'NC_000004.12', 'NC_000006.12', 'NC_000008.11', 'NC_000010.11', 'NC_000012.12', 'NC_000014.9',
               'NC_000016.10', 'NC_000018.10', 'NC_000020.11', 'NC_000022.11']

def split_train_test(chromosomes, seed):
    random.seed(seed)
    random.shuffle(chromosomes)
    mid_point = len(chromosomes) // 2
    train_set, test_set = chromosomes[:mid_point], chromosomes[mid_point:]
    return train_set, test_set

seed = sys.argv[1] #int(input("Enter a random seed: "))
train_set, test_set = split_train_test(chromosomes, seed)

with open('train.chr', 'w') as f:
    for chr in train_set:
        f.write("%s\n" % chr)

with open('test.chr', 'w') as f:
    for chr in test_set:
        f.write("%s\n" % chr)
