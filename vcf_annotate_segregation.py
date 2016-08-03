#!/usr/bin/env python3

import argparse
import sys
import gzip
import math
import re
import logging

'''
Genotypes are encoded as single integers.
0/0 = 0, 0/1 = 1, 1/1 = 2, 0/2 = 3, etc.
'''

def encode_genotype(allele1, allele2):
    try:
        allele1 = int(allele1)
        allele2 = int(allele2)
    except ValueError:
        return -1

    if allele1 > allele2:
        sys.exit("allele1 must always be less than allele2")
    i = allele2 - 1
    return allele1 + int((i + 1) * (i + 2) / 2)

def decode_genotype(gt):
    i = 0
    while ((i + 1) * (i + 2) / 2) <= gt:
        i += 1
    return (int(i - (((i + 1) * (i + 2) / 2) - gt - 1)), i)

def get_mendel_probs(fa, mo, mi_prob=0.005):
    # Basic probabilites #
    probs = [0.0, 0.0, 0.0]
    fa_alleles = [(fa & 2) >> 1, (fa & 1 | ((fa & 2) >> 1))]
    mo_alleles = [(mo & 2) >> 1, (mo & 1 | ((mo & 2) >> 1))]
    for i in range(2):
        for j in range(2):
            probs[fa_alleles[i] + mo_alleles[j]] += 0.25

    # Accounting for mendelian inconsistencies #
    for i in range(3):
        if probs[i] < 0.1: # will cause problems with mi_prob > 0.1
            if i == 0:
                if (probs[1] > 0.1):
                    probs[0] += mi_prob
                    probs[1] -= mi_prob / 2
                    probs[2] -= mi_prob / 2
                else:
                    probs[0] += mi_prob / 2
                    probs[2] -= mi_prob / 2
            elif (i == 1):
                if (probs[0] > 0.1):
                    probs[1] += mi_prob
                    probs[0] -= mi_prob
                else:
                    probs[1] += mi_prob
                    probs[2] -= mi_prob
            else:
                if (probs[1] > 0.1):
                    probs[2] += mi_prob
                    probs[0] -= mi_prob / 2
                    probs[1] -= mi_prob / 2
                else:
                    probs[2] += mi_prob / 2
                    probs[0] -= mi_prob / 2
    return probs

def update_likelihoods(ml, ol, prob, n_children, n_expected):
    ol += n_children * math.log(prob) - math.lgamma(n_children + 1)
    ml += n_expected * math.log(prob) - math.lgamma(n_expected + 1)
    return ml, ol

def phred_segregation_score(genotypes, ped_map, min_num_informative, mendel_probs):
    maximum_likelihood = math.log(1.0)
    observed_likelihood = math.log(1.0)
    num_informative = 0
    
    # Count all possible combinations at biallelic sites #
    genotype_combinations = []
    for i in range(3):
        genotype_combinations.append([])
        for j in range(3):
            genotype_combinations[i].append([])
            for k in range(3):
                genotype_combinations[i][j].append(0)

    for i in range(len(genotypes)):
        if i in ped_map:
            father_gt = genotypes[ped_map[i][0]]
            mother_gt = genotypes[ped_map[i][1]]
            proband_gt = genotypes[i]
            if father_gt == -1 or mother_gt == -1 or proband_gt == -1:
                continue
            num_informative += 1
            genotype_combinations[father_gt][mother_gt][proband_gt] += 1

    if num_informative < min_num_informative:
        return (0.0, num_informative)
    logging.debug(str(genotype_combinations))

    # Iterate over all possible genotype combinations updating likelihoods #
    for i in range(3): # father genotype
        for j in range(3): # mother genotype
            n_parents = sum(genotype_combinations[i][j])
            probs = mendel_probs[i][j]

            observed_likelihood += math.lgamma(n_parents + 1)
            maximum_likelihood += math.lgamma(n_parents + 1)

            expected_children = [int(n_parents * prob + 0.4999) for prob in probs]
            if sum(expected_children) != n_parents:
                expected_children[probs.index(max(probs))] += 1
            if sum(expected_children) != n_parents:
                logging.error("Number of expected and observed children do not match")
                sys.exit(1)

            for k in range(3): # proband genotype
                maximum_likelihood, observed_likelihood = update_likelihoods(
                    maximum_likelihood, observed_likelihood, probs[k], 
                    genotype_combinations[i][j][k], expected_children[k])

    return (((observed_likelihood - maximum_likelihood) / math.log(10) * -10), num_informative)

def write_out(file_handle, text, args):
    if args.gzip_outfile:
        file_handle.write(text.encode("ascii"))
    else:
        file_handle.write(text)

def process_args():
    parser = argparse.ArgumentParser(description="Annotate a vcf file with variant segregation scores")
    parser.add_argument("--infile")
    parser.add_argument("--infile_gzipped", action="store_true")
    parser.add_argument("--outfile")
    parser.add_argument("--gzip_outfile", action="store_true")
    parser.add_argument("--mendelian_inconsistency", type=float, default=0.005, help="The probability of encountering a mendelian inconsistency (0.0 - 0.1)")
    parser.add_argument("--debug", action="store_true")
    parser.add_argument("--min_geno_qual", type=int, default=30, help="The minimum quality genotype to include in the metric")
    parser.add_argument("--min_num_informative", type=int, default=5, help="The minimum number of informative segregation events to output a score")
    parser.add_argument("ped_file", help="A file containing pedigree information")
    return parser.parse_args()

def main(args):
    if not args:
        args = process_args()

    log_format = "%(filename)s::%(funcName)s [%(levelname)s] %(message)s"
    if args.debug:
        logging.basicConfig(level=logging.DEBUG,
                            format=log_format)
    else:
        logging.basicConfig(level=logging.INFO,
                            format=log_format)
        
    if args.mendelian_inconsistency <= 0.0 or args.mendelian_inconsistency >= 0.1:
        sys.exit("Error: --mendelian_inconsistency must be between 0.0 and 0.1 (exclusive)")
    if not args.infile and args.infile_gzipped:
        sys.exit("Error: cannot read gzipped stdin. Please decompress stdin using zcat")
    if not args.outfile and args.gzip_outfile:
        sys.exit("Error: cannot compress stdout. Please compress stdout using gzip")

    sample_ped_map = {}
    with open(args.ped_file) as f:
        for line in f:
            line = line.rstrip().split()
            if not line or not line[0]:
                continue
            if line[2] != '0' and line[3] != '0':
                sample_ped_map[line[1]] = (line[2], line[3])

    if args.infile:
        if args.infile_gzipped:
            args.infile = gzip.open(args.infile, 'rb')
        else:
            args.infile = open(args.infile, 'r')
    else:
        args.infile = sys.stdin

    if args.outfile:
        if args.gzip_outfile:
            args.outfile = gzip.open(args.outfile, 'wb')
        else:
            args.outfile = open(args.outfile, 'w')
    else:
        args.outfile = sys.stdout

    mendel_probs = []
    for i in range(3):
        mendel_probs.append([])
        for j in range(3):
            mendel_probs[i].append(get_mendel_probs(i, j, args.mendelian_inconsistency))

    logging.debug(str(mendel_probs))

    sample_col_map = {}
    col_ped_map = {}
    genotypes = []
    for line in args.infile:
        if args.infile_gzipped:
            line = line.decode("ascii")
        
        if line.startswith("##"):
            write_out(args.outfile, line, args)
            continue
        
        if line.startswith("#"):
            write_out(args.outfile, "##INFO=<ID=SEG,Number=1,Type=Float,Description=\"Phred-scaled likelihood of the genotypes being consistant with mendelian segregation\">\n", args)
            write_out(args.outfile, line, args)
            line = line.rstrip().split()
            for i, sample in enumerate(line[9:]):
                sample_col_map[sample] = i
            for sample, parents in sample_ped_map.items():
                col_ped_map[sample_col_map[sample]] = (sample_col_map[parents[0]], sample_col_map[parents[1]])
            genotypes = [0] * len(line[9:])
            continue

        line = line.split('\t')
        
        if ',' in line[4]:
            write_out(args.outfile, '\t'.join(line), args)
            continue

        try:
            gq_idx = line[8].split(':').index("GQ")
        except ValueError:
            continue

        for i, sample_gt in enumerate(line[9:]):
            sample_gt = sample_gt.split(':')
            if len(sample_gt) <= gq_idx:
                genotypes[i] = -1
            elif sample_gt[gq_idx] == '.':
                genotypes[i] = -1
            elif int(sample_gt[gq_idx]) < args.min_geno_qual:
                genotypes[i] = -1
            else:
                gt = re.split("[/|]", sample_gt[0])
                genotypes[i] = encode_genotype(gt[0], gt[1])

        segregation_score, num_informative = phred_segregation_score(genotypes, col_ped_map, args.min_num_informative, mendel_probs)
        if num_informative < args.min_num_informative:
            continue
        line[7] += ";SEG={0:.3f}".format(segregation_score + 0.000001)
        write_out(args.outfile, '\t'.join(line), args)
            
if __name__ == "__main__":
    main(None)
