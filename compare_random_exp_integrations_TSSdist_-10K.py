'''
USAGE sys.argv[0] <TSS_BED_file> <integrations_list> <mode> <number_of_simulations>

'''

import random
import sys
import matplotlib.pyplot as plt
import numpy as np
import statsmodels.stats.api as sms
import matplotlib.ticker as mticker
from matplotlib.ticker import FixedLocator, FixedFormatter
import re
import statistics
import pandas as pd
import xlsxwriter

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

    # Import TSS table and exp. integrations
    tss      = open(sys.argv[1], 'r')
    integ    = open(sys.argv[2], 'r')
    db       = sys.argv[3]
    num_iter = int(sys.argv[4])

    if db == 'all':
        number_int = 5006
        output_file_name_prefix = f'all_integrations_dist_nearestTSS_{num_iter}-simulations_10K'
    elif db == 'invitro':
        number_int = 3386
        output_file_name_prefix = f'invitro_integrations_dist_nearestTSS_{num_iter}-simulations_10K'
    elif db == 'invivo':
        number_int = 1620
        output_file_name_prefix = f'invivo_integrations_dist_nearestTSS_{num_iter}-simulations_10K'
    
    tss.readline()

    # Put TSS in a dictionary
    tss_dict = file2dict(tss, 5)

    # Put exp integrations in a dictionary
    integ_dict = file2dict(integ, 2)

    gen_length = 3088269832

    bootstrapping_list = []

    iter_number = 0
    # Handle random integrations
    for i in range(num_iter):
        iter_number += 1
        print(f'Iteration:{iter_number}')
        # Randomly pick a sequence strand
        number_picker = lambda: random.sample(range(1,3), 1) 
        number = number_picker()
        number = number[0]
        if number == 1:
            strand = '+'
        else:
            strand = '-'

        # Generate n random points along the whole hg19 genome
        randomlist_gen = lambda gen_length: random.sample(range(1, gen_length), number_int)
        random_integrations   = randomlist_gen(gen_length)

        # Determine nearest TSS
        selected_distances = []
        for integration in random_integrations:
            chrom, coord  = locate_on_chromosomes(integration, limits)
            if strand == 2:
                coord = chrom_length[chrom] - coord

            if strand == '+':
                word = 'plus'
            elif strand == '-':
                word = 'minus'

            new_key = f'{chrom}_{word}'
            nearest_tss = 248956422
            for coordinate in tss_dict[new_key]:
                    dist = int(coordinate) - coord
                    if abs(dist) < abs(nearest_tss):
                        nearest_tss = dist 

            selected_distances.append(nearest_tss)


        sim_distances_file = f'sim_{output_file_name_prefix}.tsv'
        out_sim = open(sim_distances_file, 'w')
        lists_dict_sim = sort_by_bin(selected_distances)

        for key in lists_dict_sim:
            out_sim.write(f'{key}\n')
            for integ in lists_dict_sim[key]:
                out_sim.write(f'{integ}\n')

        out_sim.close()

        hist, bin_edges = np.histogram(selected_distances, density = False, bins = 10, range=[-10000, 0])    
        bootstrapping_list.append(hist)
        ### bootstrapping_list.append(truncated_list)
    
    list_of_averages = list(map(lambda x: sum(x)/len(x), zip(*bootstrapping_list)))
    # compute a 95%CI modelling the distribution as a Gaussian
    list_of_stdev = list(map(statistics.stdev, zip(*bootstrapping_list)))
    list_of_lower_CI = [mean - 1.95996398 * stdev  for mean, stdev in zip(list_of_averages, list_of_stdev)]
    list_of_upper_CI = [mean + 1.95996398 * stdev  for mean, stdev in zip(list_of_averages, list_of_stdev)]
    ### hist, bin_edges = np.histogram(selected_distances, density=False, bins=40)

    # Handle experimental integrations
    selected_distances_exp = []

    for key in integ_dict:
        for integration in integ_dict[key]:
            if re.search('chr.*plus', key):
                chrom = key.replace('_plus', '')
            if re.search('chr.*minus', key):
                chrom = key.replace('_minus', '')
            coord = integration

            nearest_tss = 248956422
            for coordinate in tss_dict[key]:
                dist = int(coordinate) - int(coord)
                if abs(dist) < abs(nearest_tss):
                    nearest_tss = dist

            selected_distances_exp.append(nearest_tss)

    exp_distances_file = f'exp_{output_file_name_prefix}.tsv'
    out_exp = open(exp_distances_file, 'w')
    lists_dict_exp = sort_by_bin(selected_distances_exp)

    for key in lists_dict_exp:
        out_exp.write(f'{key}\n')
        for integ in lists_dict_exp[key]:
            out_exp.write(f'{integ}\n')

    out_exp.close()


    hist_exp, bin_edges_exp = np.histogram(selected_distances_exp, density = False, bins = 10, range=[-10000, 0])

    list_of_lists= []
    for i in range(10):
        my_list = [list_of_averages[i], list_of_stdev[i], list_of_lower_CI[i], list_of_upper_CI[i], hist_exp[i]]
        list_of_lists.append(my_list)

    df = pd.DataFrame.from_records(list_of_lists, columns=['Average_sim','StDev_sim','CI_lower_limit','CI_upper_limit','Freq_exp'])
    stats_file = f'{output_file_name_prefix}.xlsx'
    printer(stats_file, df)


    fig, ax = plt.subplots() 
    kwargs = dict(histtype='bar', color='green', density=False, bins=10,  range=[-10000, 0],
        edgecolor='black', linewidth=1)
    ax.hist(selected_distances_exp, **kwargs)

    kwargs = dict(histtype='step', color='white', density=False, bins=10,  range=[-10000, 0],
        edgecolor='black', linewidth=1, alpha=0.9)
    ### ax.hist(list_of_averages, **kwargs)
    ax.hist(range(-9700, 0, 1000), weights = list_of_averages, **kwargs)
    ax.hist(range(-9700, 0, 1000), weights = list_of_lower_CI, **kwargs)
    ax.hist(range(-9700, 0, 1000), weights = list_of_upper_CI, **kwargs)

    new_ticks_labels = ['-10', '-9', '-8', '-7', '-6', '-5', '-4', '-3', '-2', '-1', '0']
    x_locator   = FixedLocator([-10000, -9000, -8000, -7000, -6000, -5000, -4000, -3000, -2000, -1000, 0])
    x_formatter = FixedFormatter(new_ticks_labels)
    ax.xaxis.set_major_locator(x_locator)
    ax.xaxis.set_major_formatter(x_formatter)
    plt.xlabel("Nearest TSS (Kb)")
    plt.ylabel("Frequency")

    fig.savefig(f'{output_file_name_prefix}.png')

def sort_by_bin(integrations_list):
    dict_of_lists = {}
    for i in range (-10000, 0, 1000):
        minimum = i
        maximum = i + 1000
        key = f'{minimum}-{maximum}'
        dict_of_lists[key] = []

        for integ in integrations_list:
            if integ > minimum and integ < maximum + 1:
                dict_of_lists[key].append(integ)

    return dict_of_lists

def printer(filename, df):
  writer = pd.ExcelWriter(filename, engine='xlsxwriter')
  df.to_excel(writer, index=False)
  writer.save()

def file2dict(infile, strand_col):
    infile.readline()
    my_dict = {}

    while True:
        line = infile.readline().strip()

        if not line:
            break

        fields = line.split('\t')
        if fields[0]   == 'chrX':
            fields[0]   = 'chr23'
        elif fields[0] == 'chrY':
            fields[0]   = 'chr24'
        if fields[strand_col]   == '+':
            word = 'plus'
        elif fields[strand_col] == '-':
            word = 'minus'
        my_key = f'{fields[0]}_{word}'
        if not my_key in my_dict:
            my_dict[my_key] = []
        my_dict[my_key].append(fields[1])

    infile.close()
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
