#! /usr/bin/python
import argparse
import random
import allel
import msprime
import numpy as np
import pyslim
import os


parser = argparse.ArgumentParser(
    description='Input a SLiM TreeSeq, add mutations, downsample to correct population '
                'sampling configuration and make an input file ready for 3pclr'
                'need to provide -ts slim'
)
parser.add_argument(
    "-indir",
    required=True,
    help="name of the input dir we want to read from"
)
parser.add_argument(
    "-slim_ts",
    required=True,
    help="name of the treeSeq file after slim forward simulation"
)
# parser.add_argument(
#     "-sim_sample_config",
#     nargs="*",  # 0 or more values expected => creates a list
#     type=int,
#     default=[1, 36, 38, 20, 22], # default if nothing is provided,
#     help = "haploid sample configuration to produce. Default is 0, 36, 38, 20, 22 for B,C,E,N,W"
# )
parser.add_argument(
    "-outdir",
    required=True,
    help="name of the output dir we want to write to"
)
parser.add_argument(
    "-output",
    help="name of the output TreeSeq file we want to create"
)
args = parser.parse_args()


def downsample_tree_seq(ts, sim_sample_config):
    pops = list(range(len(sim_sample_config)))
    #pops = [1, 2, 3, 4]
    all_pop_sampled_nodes = []
    for j, pop in enumerate(pops):
        node_list = []
        for n in ts.nodes():
            if n.population == pop and n.time == 0.0:
                node_list.append(n.id)
        sampled_nodes = random.sample(node_list, sim_sample_config[j])
        all_pop_sampled_nodes.extend(sampled_nodes)
    all_pop_sampled_nodes.sort()
    ts_reduced = ts.simplify(samples=all_pop_sampled_nodes)
    return ts_reduced


def get_tree_seq_add_mutations_and_downsample(sim_sample_config, infile):
    # lets read in the slim TreeSeq, and add some mutations!
    ts = pyslim.load(infile).simplify()
    mutated = pyslim.SlimTreeSequence(msprime.mutate(ts, rate=1.2e-8, keep=True))  # no random seed argument, wrap in psl
    ts = downsample_tree_seq(ts=mutated, sim_sample_config=sim_sample_config)
    return ts

## check position of selected site....
# for s in ts.sites():
#     print(s)

def threep_clr_output(ts, outfile):
    header = 'chr\tphyspos\tgenpos\tmBonobo\tnBonobo\tmCentral\tnCentral\tmEastern\tnEastern\tmNigeria\tnNigeria\tmWestern\tnWestern'
    #print("outfile")
    print(header, file = open(outfile, "w"))
    chr = 23
    v = np.zeros((ts.get_num_mutations(), ts.get_sample_size()),
                 dtype=np.int8)
    for variant in ts.variants():
        v[variant.index] = variant.genotypes
    pops = [pop.id for pop in ts.populations()]
    pop_indices = [ts.samples(population=pop) for pop in pops]
    # create a haplotype array in allele from numpy arrays of 0s/1s
    # as the dataset is reasonably small there is no need to use the chunking functionality of allele
    haplotypes = allel.HaplotypeArray(v)
    # create a list of allele counts for each population
    allele_counts = [haplotypes.count_alleles(max_allele=None, subpop=pop_idx) for pop_idx in pop_indices]
    # need to keep track of sites that are duplicated.....
    previous_sites = []
    for j, site in enumerate(ts.sites()):
        physpos = int(round(site.position))
        while physpos in previous_sites:
            physpos = physpos + 1
        previous_sites.append(physpos)
        genpos = round(physpos / 1e6 * 0.96, 6)
        outline = [chr, physpos, genpos] + [[row[1], row[0] + row[1]] for row in allele_counts[0]][j] + \
                  [[row[1], row[0] + row[1]] for row in allele_counts[1]][j] + \
                  [[row[1], row[0] + row[1]] for row in allele_counts[2]][j] + \
                  [[row[1], row[0] + row[1]] for row in allele_counts[3]][j] + \
                  [[row[1], row[0] + row[1]] for row in allele_counts[4]][j]
        print('\t'.join(str(e) for e in outline), file = open(outfile, "a"))
        #print('\t'.join(str(e) for e in outline))


def main():
    sim_sample_config = [1, 36, 38, 20, 22]
    #tree_seq = get_tree_seq_add_mutations_and_downsample(sim_sample_config= sim_sample_config, infile=os.path.join(format(args["indir"]), format(args["slim_ts"])))
    #threep_clr_output(tree_seq, outfile=os.path.join(format(args["outdir"]), format(args["outfile"])))
    tree_seq = get_tree_seq_add_mutations_and_downsample(sim_sample_config=sim_sample_config,
                                                         infile=os.path.join(format(args.indir), format(args.slim_ts)))
    threep_clr_output(tree_seq, outfile=os.path.join(format(args.outdir), format(args.output)))


if __name__ == "__main__":
    main()