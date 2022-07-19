#! /usr/bin/python
import msprime
import math
import pyslim
import os
import argparse


parser = argparse.ArgumentParser(description='generate a coalescent prehistory of chimp demography')
parser.add_argument("-sim_file_name", required=True,
help="name of the output file we want to create")
parser.add_argument("-outdir", required=True,
help="name of the output dir we want to write to")
parser.add_argument("-size", required=True,
help="size in Mbp, of the locus to simulate")
args = vars(parser.parse_args())


def truncated_pan_demography(sim_locus_length = round(float(args["size"]) * 1.0e6)):
    # we only have to have two demes: bonobo and chimp_common_ancestor
    # in the coalescent simulation the first specification of the bonobo Ne
    # is n_421kya_bonobo: int = round(0.299481445247702 * n_ms)  # 2995
    # this means that bonobo size from 421 back to the split with chimpanzees is
    # 2995. Note that 421kya is after the split of chimp lineages into central/eastern and
    # nigeria/western, so in our truncated coalescent, this is the initial size of bonobo....
    # for the common ancestor of all chimps, there is a bottleneck associated with the split with bonobo
    # bottleneck is for 2000 years (1598kya back to 1600kya)
    # t_1598kya: int = 1598800 / generation_time
    # n_1598kya_all_chimps_CA: int = round(0.00336130452736601 * n_ms) # 34
    # there is then an expansion, so that the size of the chimp common ancestor at the time of the split
    # is 16692
    # t_482kya: int = 482625 / generation_time
    # n_482kya_all_chimps_ca: int = round(1.66920782430592 * n_ms) # 16692
    # there is migration between the common ancestor of chimps and bonobo.
    # the merger of bonobo and chimps is at 1600800, and the size of the pan ancestor is 14711
    # t_1600kya: int = 1600800 / generation_time
    # n_1600kya_all_pan_CA: int = round(1.47105091660349 * n_ms) # 14711
    n_ms: int = 10000
    generation_time: int = 25
    m_scale = 4 * n_ms
    # need a scalar to shift times as we have a truncated simulaiton
    #
    t_shift = 482625 / generation_time
    n_present_bonobo: int = round(0.299481445247702 * n_ms) # 2995
    n_present_all_chimps_ca: int = round(1.66920782430592 * n_ms) # 16692
    t_pan_split_all_chimps_bottleneck: int = (1598800 / generation_time) - t_shift
    n_pan_split_all_chimps_ca: int = round(0.00336130452736601 * n_ms) # 34
    n_pan_ancestor: int = round(1.47105091660349 * n_ms) # 14711
    t_pan_split: int = (1600800 / generation_time) - t_shift
    migration_matrix = [
            [0, 0.241282075772286/ m_scale],
            [0.0101771164248256/ m_scale, 0]
        ]
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=n_present_bonobo * 2, initial_size=n_present_bonobo),
        msprime.PopulationConfiguration(
            sample_size=n_present_all_chimps_ca * 2, initial_size=n_present_all_chimps_ca)
        ]
    demographic_events = [
        # define all demographic events into the past....
        # Ne decreases in chimp ancestor 2000 years/80 generations after the pan split.
        msprime.PopulationParametersChange(
            time=t_pan_split_all_chimps_bottleneck, initial_size=n_pan_split_all_chimps_ca, growth_rate=0, population_id=1),
        # merge bonobo to chimp ancestor
        msprime.MassMigration(
            time=t_pan_split, source=0, destination=1, proportion=1.0),
        # turn migration off for now only 1 lineage
            msprime.MigrationRateChange(time=t_pan_split, rate=0),
        # set size of pan ancestor from pan split backwards in time
            msprime.PopulationParametersChange(
                time=t_pan_split, initial_size=n_pan_ancestor)
    ]
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    #dd = msprime.DemographyDebugger(
    #    population_configurations = population_configurations,
    #    migration_matrix = migration_matrix,
    #    demographic_events = demographic_events)
    #dd.print_history()
    tree_seq: object = msprime.simulate(
        population_configurations=population_configurations,
        migration_matrix=migration_matrix,
        demographic_events=demographic_events,
        length=sim_locus_length,
        recombination_rate=0.96e-8,
        mutation_rate=0
    )
    return(tree_seq)

def main():
    tree_seq = truncated_pan_demography()
    # num_replicates=num_replicates)
    slim_tree_seq = pyslim.annotate_defaults(tree_seq, model_type="WF", slim_generation=1)
    suffix = 'trees'
    slim_tree_seq.dump(os.path.join(format(args["outdir"]), format(args["sim_file_name"]) + "." + suffix))


if __name__ == "__main__":
    main()






