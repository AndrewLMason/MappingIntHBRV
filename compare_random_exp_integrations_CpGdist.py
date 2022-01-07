import random
import sys
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.stats.api as sms
import matplotlib.ticker as mticker
from matplotlib.ticker import FixedLocator, FixedFormatter
import re
import statistics

'''
    USAGE: python compare_random_exp_integrations_CpGdist.py <CpG_ref_file> <exp_integrations_file> <database>
'''  


def main():
    # Length of hg19 chromosomes
    chrom_length = {
        'chr1': 248956422,
        'chr2': 242193529,
        'chr3': 198295559,
        'chr4': 190214555,
        'chr5': 181538259,
        'chr6': 170805979,
        'chr7': 159345973,
        'chr8': 145138636,
        'chr9': 138394717,
        'chr10': 133797422,
        'chr11': 135086622,
        'chr12': 133275309,
        'chr13': 114364328,
        'chr14': 107043718,
        'chr15': 101991189,
        'chr16': 90338345,
        'chr17': 83257441,
        'chr18': 80373285,
        'chr19': 58617616,
        'chr20': 64444167,
        'chr21': 46709983,
        'chr22': 50818468,
        'chr23': 156040895,
        'chr24': 57227415
    }

    # Define limits of each chromosome as the length of
    # all previous chromosomes plus the length of the actual
    # chromosome. i.e. limit of chr2 is chr1 length +
    # chr2 length
    limits = {}

    for key in chrom_length:
        if key == 'chr1':
            limits[key] = 248956422
        else:
            prev_key = key.replace('chr', '')
            prev_key = int(prev_key) - 1
            prev_key = f'chr{prev_key}'

            limits[key] = chrom_length[key] + limits[prev_key]

    # Import CpG table and exp. integrations
    cpg   = open(sys.argv[1], 'r')
    integ = open(sys.argv[2], 'r')
    db    = sys.argv[3]

    if db == 'all':
        number_int = 5006
    elif db == 'invitro':
        number_int = 3386
    elif db == 'invivo':
        number_int = 1620

    # Put CpG reference coordinates in a dictionary
    cpg_dict = file2dict(cpg, 'ref')

    # Put exp integrations in a dictionary
    integ_dict = file2dict(integ, 'exp')


    gen_length = 3088269832

    bootstrapping_list = []

    # Handle random integrations
    iter_counter = 0
    for i in range(2):
        iter_counter += 1
        print(f'Iteration:{iter_counter}')
        # Generate random points along the whole hg19 genome
        randomlist_gen = lambda gen_length: random.sample(range(1, gen_length), number_int)
        random_integrations   = randomlist_gen(gen_length)

        # Determine nearest CpG
        selected_distances = []
        for integration in random_integrations:
            chrom, coord  = locate_on_chromosomes(integration, limits)
            nearest_cpg = find_nearest_feature(cpg_dict, chrom, coord)
            selected_distances.append(nearest_cpg)

        out_sim = open("distances_sim_to_CpG_per_bin_onlyPosDist.tsv", 'w')
        lists_dict_sim = sort_by_bin(selected_distances)

        for key in lists_dict_sim:
            out_sim.write(f'{key}\n')
            for integ in lists_dict_sim[key]:
                out_sim.write(f'{integ}\n')

        out_sim.close()

        hist, bin_edges = np.histogram(selected_distances, density = False, bins = 20, range=[0, 20000])    
        bootstrapping_list.append(hist)
    
    list_of_averages = list(map(lambda x: sum(x)/len(x), zip(*bootstrapping_list)))

    # compute a 95%CI modelling the distribution as a Gaussian
    list_of_stdev    = list(map(statistics.stdev, zip(*bootstrapping_list)))
    list_of_lower_CI = [mean - 1.95996398 * stdev  for mean, stdev in zip(list_of_averages, list_of_stdev)]
    list_of_upper_CI = [mean + 1.95996398 * stdev  for mean, stdev in zip(list_of_averages, list_of_stdev)]

    # Handle experimental integrations
    selected_distances_exp = []

    for key in integ_dict:
        fields = key.split('_')
        chrom  = fields[0]
        coord  = fields[1]
        strand = fields[2]
        if strand == '-':
            coord = chrom_length[chrom] - coord
        nearest_cpg = find_nearest_feature(cpg_dict, chrom, coord)
        selected_distances_exp.append(nearest_cpg)

    out_exp = open("distances_exp_to_CpG_per_bin_onlyPosDist.tsv", 'w')
    lists_dict_exp = sort_by_bin(selected_distances_exp)

    for key in lists_dict_exp:
        out_exp.write(f'{key}\n')
        for integ in lists_dict_exp[key]:
            out_exp.write(f'{integ}\n')

    out_exp.close()


    fig, ax = plt.subplots() 
    kwargs = dict(histtype='bar', color='green', density=False, bins=20,  range=[0, 20000],
        edgecolor='black', linewidth=1)
    ax.hist(selected_distances_exp, **kwargs)

    kwargs = dict(histtype='step', color='white', density=False, bins=20,  range=[0, 20000],
        edgecolor='black', linewidth=1, alpha=0.9)
    
    ax.hist(range(0, 20000, 1000), weights = list_of_averages, **kwargs)
    ax.hist(range(0, 20000, 1000), weights = list_of_lower_CI, **kwargs)
    ax.hist(range(0, 20000, 1000), weights = list_of_upper_CI, **kwargs)

    new_ticks_labels = ['0', '5', '10', '15', '20']
    x_locator   = FixedLocator([0, 5000, 10000, 15000, 20000])
    x_formatter = FixedFormatter(new_ticks_labels)
    ax.xaxis.set_major_locator(x_locator)
    ax.xaxis.set_major_formatter(x_formatter)
    plt.xlabel("Nearest CpG Island (Kb)")
    plt.ylabel("Frequency")

    fig.savefig("retrovirus_integrations_nearestCpG_1000-simulations_onlyPosDist.png")

def sort_by_bin(integrations_list):
    dict_of_lists = {}
    for i in range (-200000, 200000, 10000): 
        minimum = i
        maximum = i + 10000
        key = f'{minimum}-{maximum}'
        dict_of_lists[key] = []

        for integ in integrations_list:
            if integ > minimum and integ < maximum + 1:
                dict_of_lists[key].append(integ)

    return dict_of_lists


def find_nearest_feature(dictionary, chrom, coord):
    nearest_cpg = 248956422
    for key in dictionary:
        fields = key.split('_')
        my_chrom = fields[0]
        my_coord = fields[1]

        if my_chrom == chrom:
            dist = float(my_coord) - int(coord)
            if abs(dist) < abs(nearest_cpg):
                nearest_cpg = dist
    return abs(nearest_cpg)

def rename_chromosome (chr_name):
    if chr_name == 'chrX':
        chr_name = 'chr23'
    elif chr_name == 'chrY':
        chr_name = 'chr24'
    return chr_name


def file2dict(integrations_file, int_type):
    integrations_file.readline()
    # int_type can be 'ref' or 'exp'
    my_dict = {}
    while True:
        line = integrations_file.readline()
        if not line:
            break

        fields   = line.split('\t')
        if int_type == 'ref':
            cpg_mean = (int(fields[2]) + int(fields[3])) / 2
            chrom = rename_chromosome(fields[1])
            my_key = f'{chrom}_{cpg_mean}'
        elif int_type == 'exp':
            chrom = rename_chromosome(fields[0])
            my_key = f'{chrom}_{fields[1]}_{fields[2]}'
        else:
            print("Please input a valid integration type: 'ref' or 'exp'")

        if not my_key in my_dict:
            my_dict[my_key] = 1
        else:
            print("Repeated key")

    return my_dict

def locate_on_chromosomes(gen_position, limits): 

    if gen_position <= 248956422:
        chrom = 'chr1'
        coord = gen_position
        return chrom, coord

    for key in limits:
        if key != 'chr1':
            prev_key = key.replace('chr', '')
            prev_key = int(prev_key) - 1
            prev_key = f'chr{prev_key}'
            if gen_position > int(limits[prev_key]) and gen_position < int(limits[key]):
                chrom = key
                coord = gen_position - int(limits[prev_key])
                return chrom, coord
    
if __name__ == '__main__': 
    main()
