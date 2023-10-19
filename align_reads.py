## Question 1

## Importing the required modules
import numpy as np
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


########################################################################################################################

## function called “assemble reads” which takes as input a file with reads and
## outputs a file with a FASTA sequence of contigs labeled numerically


def assemble_reads(file_path, k):
    # """
    # Function to get contig from a file of sequence reads
    # :param file_path: path to file
    # :param k: minimum length of overlap, with asssumption to align reads if k or more overlap
    # :return: fasta file with numbered contigs
    # """

    # reading the sequence reads using helper function get reads
    # list of reads stored in variable seq_lis
    seq_lis = _get_reads(file_path)


    ## getting the list of the contigs with k overlap(assumption k=10, or more) using helper function get contigs
    ## and check contigs

    final_contigs = _get_con(seq_lis, k=k)

    ## writing the contigs to a fasta file using seqio.write function from biopython module after craeting
    # seq objects and seq recors

    # list to create seq records/seq objects
    cont_lis = []
    for contig in range(len(final_contigs)):
        cont = Seq(final_contigs[contig])
        rec = SeqRecord(cont, f'Contig{contig+1}', description=f"this is contig {contig + 1}")
        cont_lis.append(rec)

    SeqIO.write(cont_lis, "Teena_Assignment1.fa", "fasta")

    ## also returning lsit of contigs to used for visualization
    return final_contigs








###########################################

## helper functions

## 1.  Function to read the file with path for the file as the parameter

def _get_reads(path):
    """
    Function to read read sequences from the give file
    :param path: path to the file
    :return: list of reads
    """
    seq_file = open(path, "r")

    seq_list = []

    for line in seq_file:
        line = line.rstrip()
        seq_list.append(line)

    return seq_list




## 2. function to merge contigs with a specific number of base overlaps

def _get_contigs(sort_seq, min_k):
    """
    function to merge reads with specified overlap
    :param sort_seq: list of sequence reads
    :param min_k: minimum over lap specified in argument
    :return: mergerd reads into contigs for further processing
    """

    contig_list = []
    i = 0

    j = len(sort_seq) - 1

    while i != len(sort_seq):
        if j == i:
            contig_list.append(sort_seq[i])
            #print("yes")
            i = i + 1
            j = len(sort_seq) - 1
            continue

        ## checking the overlap at the end and the start of the reads
        v = sort_seq[i][-min_k:]
        e = sort_seq[j][:min_k]

        q = sort_seq[i][:min_k]

        a = sort_seq[j][-min_k:]

        if (v == e):
            contig_list.append(sort_seq[i] + sort_seq[j][min_k:])
            i = i + 1
            del (sort_seq[j])
            j = len(sort_seq) - 1
            continue

        if (q == a):
            contig_list.append(sort_seq[j] + sort_seq[i][min_k:])
            i = i + 1
            del (sort_seq[j])
            j = len(sort_seq) - 1
            continue

        j = j - 1

    return (contig_list)



## 3. helper function to check the further overlaps until there are none with itearation with same number of overlap of bases
## specified in the argument

def _check_contig(s_list, k):
    """
    function to check if reads can be merged with same k again
    :param s_list: read list
    :param k: overlap of bases at the end and start of the reads
    :return: the final contig for each k specified
    """
    num = 0
    new_num = -1

    while new_num < num :
        num_lis = _get_contigs(s_list, k)
        num = len(num_lis)

        s_list = num_lis

        new_lis = _get_contigs(s_list, k)
        new_num = len(new_lis)
        s_list = new_lis


    if num == new_num:
        return new_lis



## 4. function to get final contig list with the folllowing assumptions
# ## k =10. base pair overlap k or greater implies reads align
# # No divergent overlaps
# # A->B; A->C but B!=C

def _get_con(seq_lis, k):
    """
    function to get contigs with specified k or more overlaps
    :param seq_lis: list of reads
    :param min_num:
    :return:
    """
    check_len = len(max(seq_lis))
    i = check_len
    while i != k-1:
        a = _get_contigs(seq_lis, i)
        new_list = _check_contig(a, i)
        seq_lis = new_list

        i -= 1

    if i == k-1:
        return seq_lis






####################################
## running the main function assemble_reads with k=10

path = "/home/teena/Desktop/binf6400/Assignment 1/seqReadFile.txt"

print(assemble_reads(path, k=10))









###########################################
## creating the visualization for mapping the reads on the contigs to check the coverage

## checking the mapping of reads for all the contigs to get coverage
# list of contigs
contig_list = assemble_reads(path, 10)

# list of seq reads
seq_reads = _get_reads(path)


# initializing dictionary to store the read coverage for all the the contigs for further visualization
cont_dict = {}


for contg in contig_list:
    empty_list = [0] * len(contig_list[0])

    # array to store coverage
    j = np.array(empty_list)
    for read in seq_reads:
        if read in contg:
            # getting the first position of reads
            initial = contg.index(read)
            end = initial + len(read)

        # incrementig coverage if there is a read
            j[initial:end] += 1

            cont_dict[contg] = j

##print(cont_dict)


## visualization
# x-axis with length of the contig
conti_lis = list(range(len(contig_list[0])))

# loop to show coverage for all the contigs
for cont in range(len(cont_dict)):
    plt.plot(conti_lis, cont_dict[contig_list[cont]])

plt.xticks(conti_lis, rotation=90, fontsize = 'xx-small')

plt.show()

## visualization of coverage in two contig craeted by k=10
#fig, axes = plt.subplot(2)

plt.plot(conti_lis, cont_dict[contig_list[0]])
plt.plot(conti_lis, cont_dict[contig_list[1]])

plt.plot(conti_lis, cont_dict[contig_list[2]])
plt.plot(conti_lis, cont_dict[contig_list[3]])

plt.xticks(conti_lis, rotation=90, fontsize = 'xx-small')

plt.show()












########################################################################################################################
## Part 2



# effect of increasing the value of k

increase_k = {}

for k in range(11, 21):
    klen = len(assemble_reads(path, k))
    increase_k[f'k{k}'] = klen


print(increase_k)

## the number of contigs generated by incraesing the number of overlap incraeses gradually with significant incrase
## after overlap of 15 or more bases


## effect of decreasing the value of k
decrease_k = {}

for k in range(1, 10):
    klen = len(assemble_reads(path, k))
    decrease_k[f'k{k}'] = klen


print(decrease_k)

## the number of contigs generated with overlap from 5-10 or above till 21 are almost same, but the number of contigs decreases
## after overlap of 4 or less


## effect on reads coverage with incraesed k and decreased k

# checking the coverage with k of 2 where 9 contigs are created


contig_k2 = assemble_reads(path, 2)

# list of seq reads
seq_reads = _get_reads(path)


## to store coverage for contigs with k 2
cont_k2 = {}


for contg in contig_k2:
    empty_list = [0] * len(contig_k2[0])

    # array to store coverage
    j = np.array(empty_list)
    for read in seq_reads:
        if read in contg:
            # getting the first position of reads
            initial = contg.index(read)
            end = initial + len(read)

        # incrementig coverage if there is a read
            j[initial:end] += 1

            cont_k2[contg] = j




## visualization with decreasing k
conti_lis = list(range(len(contig_k2[0])))

# loop to show coverage for all the contigs
#for cont in range(len(cont_k2)):
plt.plot(conti_lis, cont_k2[contig_k2[0]], label="contig1")
plt.plot(conti_lis, cont_k2[contig_k2[1]], label="contig2")
plt.plot(conti_lis, cont_k2[contig_k2[2]], label="contig3")


plt.legend(loc="upper right")
plt.xticks(conti_lis, rotation=90, fontsize = 'xx-small')

plt.show()

##################################

# checking the coverage with k of 15


contig_k15 = assemble_reads(path, 15)

# list of seq reads
seq_reads = _get_reads(path)


## to store coverage for contigs with k 2
cont_k15 = {}


for contg in contig_k15:
    empty_list = [0] * len(contig_k15[0])

    # array to store coverage
    j = np.array(empty_list)
    for read in seq_reads:
        if read in contg:
            # getting the first position of reads
            initial = contg.index(read)
            end = initial + len(read)

        # incrementig coverage if there is a read
            j[initial:end] += 1

            cont_k15[contg] = j




## visualization with increasing k
conti_lis = list(range(len(contig_k15[0])))

# loop to show coverage for all the contigs
#for cont in range(len(cont_k2)):
for i in range(0, 9):
    plt.plot(conti_lis, cont_k15[contig_k15[i]], label=f'contig{i}')


plt.xticks(conti_lis, rotation=90, fontsize = 'xx-small')
plt.legend(loc="upper right")
plt.show()