# calculate mutation frequency of covsonar outputs #

import pandas as pd
import sys
import argparse

def extract_profiles(df, profile):
    allowed_dna = ["DNA", "dna", "DNS", "dns"]
    allowed_as = ["aa", "as", "AA", "AS", "aminoacid"]
    if profile in allowed_dna:
        mut_profiles = df["dna_profile"]
    elif profile in allowed_as:
        mut_profiles = df["aa_profile"]
    else:
        sys.exit("no such profile: use aa or dna")

    return mut_profiles

def calc_freq(df, profile, threshold):

    frequencys = dict()
    consensus_mutations = dict()

    mut_profiles = extract_profiles(df, profile)

    n_entrys = len(mut_profiles)

    for profile in mut_profiles:
        mutations = profile.split(" ")
        for mut in mutations:
            if mut in frequencys:
                frequencys[mut] += 1
            else:
                frequencys[mut] = 1

    print("# n of sequences: ", n_entrys)
    print("mutation", "frequency", sep = "\t")
    for mut in frequencys:
        frequencys[mut] = frequencys[mut]/n_entrys
        if threshold > 1 or threshold < 0:
            sys.exit("threshold must be between 0 and 1")
        if frequencys[mut] >= threshold:
            print(mut, round(frequencys[mut],2), sep = "\t")

def filter_lineages(df, lineages):
    to_drop = []
    for i,entry in df.iterrows():
        if entry["lineage"] not in lineages:
            to_drop.append(i)
    df = df.drop(index = to_drop)
    if len(df) == 0:
        sys.exit("none of the specified lineages found")

    return df

def calculate_consensus_mutations(csv, lineages = "", combine_lineages = False, profile = "dna", threshold = 0.75):

    covsonar_df = pd.read_csv(csv)

    if lineages:
        separators = ["; ",", "," ","\t",";",","," ","\t"]
        for sep in separators:
            if sep in lineages:
                lineages = lineages.split(sep)
                lineages = set(lineages)
                break
        else:
            lineages = set([lineages])

    if lineages and not combine_lineages:

        covsonar_df = filter_lineages(covsonar_df, lineages)

        for lineage in lineages:
            covsonar_df_temp = covsonar_df[covsonar_df["lineage"] == lineage]
            print("", "", sep="\t")
            print("lineage:", lineage, sep = "\t")
            calc_freq(covsonar_df_temp, profile, threshold)

    elif lineages and combine_lineages:
        covsonar_df = filter_lineages(covsonar_df, lineages)
        calc_freq(covsonar_df, profile, threshold)
    else:
        calc_freq(covsonar_df, profile, threshold)




if __name__ == "__main__":

    parser = argparse.ArgumentParser()

    parser.add_argument("infile", help="Filepath to covsonar.csv")
    parser.add_argument("-l", "--lineages", type = str, default = "", help="lineages to calculate consensus mutations")
    parser.add_argument("-p", "--profile", type = str, default="dna", help="dna or aa")
    parser.add_argument("-t", "--threshold", type = float, default=0.75, help="threshold to consider as consensus")
    parser.add_argument("-c", "--combine-lineages", type = bool, default = False, help="one single consensus for defined lineages")

    args = parser.parse_args()

    calculate_consensus_mutations(args.infile, lineages = args.lineages, combine_lineages = args.combine_lineages, profile = args.profile, threshold = args.threshold)
