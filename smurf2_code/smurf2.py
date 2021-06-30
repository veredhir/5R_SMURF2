#!/usr/bin/env python2


INRTO = \
"""
SMURF2: Short MUltiple Regions Framework for Reconstruction of Genes from the Environment using
Expectation-Maximization Iterative method.

Some of the code is based on EMIRGE:
EMIRGE: Expectation-Maximization Iterative Reconstruction of Genes from the Environment
Copyright (C) 2010-2016 Christopher S. Miller  (christopher.s.miller@ucdenver.edu)
https://github.com/csmiller/EMIRGE
"""

USAGE = \
    """usage: %prog DIR <required_parameters> [options]

This version of SMURF2 (%prog) attempts to reconstruct rRNA SSU genes
from Illumina amplicon data based on the 5R protocol.

DIR is the working directory to process data in.
Use --help to see a list of required and optional arguments

Fuks G, Elgart M, Amir A, Zeisel A, Turnbaugh PJ, Soen Y, Shental N. 
Combining 16S rRNA gene variable regions enables high-resolution microbial community profiling. 
Microbiome. 2018 Jan 26;6(1):17. doi: 10.1186/s40168-017-0396-x. PMID: 29373999; PMCID: PMC5787238.

Miller CS, Baker BJ, Thomas BC, Singer SW, Banfield JF (2011)
EMIRGE: reconstruction of full-length ribosomal genes from microbial community short read sequencing data.
Genome biology 12: R44. doi:10.1186/gb-2011-12-5-r44.

Miller CS, Handley KM, Wrighton KC, Frischkorn KR, Thomas BC, Banfield JF (2013)
Short-Read Assembly of Full-Length 16S Amplicons Reveals Bacterial Diversity in Subsurface Sediments.
PloS one 8: e56018. doi:10.1371/journal.pone.0056018.
"""

import sys
import os

import shutil
from optparse import OptionParser, OptionGroup, SUPPRESS_HELP
from smurf2_iteration import Smurf2Iteration
from post_process import convert_to_smurf_format
from smurf2_utills import *


class EM(object):
    """
    driver class for EM algorithm
    """

    def __init__(self,
                 working_dir,
                 is_debug):
        self.working_dir = working_dir
        self.is_debug = is_debug

    @time_it
    def execute(self,
                max_iter,
                fastq,
                fastq_reversed,
                primers_path,
                fasta,
                read_length,
                number_of_region,
                update_weight_using_the_reads,
                max_changed_bases_rate_for_split,
                min_coverage_for_split,
                min_minor_prob_for_split,
                min_similar_bases_rate_for_merge,
                max_priors_diff_for_stability_test,
                allow_split,
                taxa_path,
                matlab_path):

        first_subdir = os.path.join(self.working_dir, "iter.%02d" % 0)

        os.mkdir(first_subdir)

        curr_emirge_iteration = Smurf2Iteration(first_subdir,
                                                fastq_path=fastq,
                                                reversed_fastq_path=fastq_reversed,
                                                fasta_path=fasta,
                                                primers_path=primers_path,
                                                read_len=read_length,
                                                n_regions=number_of_region,
                                                update_weight_using_the_reads=update_weight_using_the_reads,
                                                max_changed_bases_rate_for_split=max_changed_bases_rate_for_split,
                                                min_coverage_for_split=min_coverage_for_split,
                                                min_minor_prob_for_split=min_minor_prob_for_split,
                                                min_similar_bases_rate_for_merge=min_similar_bases_rate_for_merge,
                                                max_priors_diff_for_stability_test=max_priors_diff_for_stability_test,
                                                allow_split=allow_split,
                                                debug_mode=self.is_debug)
        curr_emirge_iteration.do_iteration(False)
        subdirs = []
        subdirs.append(first_subdir)
        for i in range(1, max_iter):
            subdirs.append(os.path.join(self.working_dir, "iter.%02d" % i))
            os.mkdir(subdirs[i])
            logging.info("Create iteration {} working dir = {}".format(i, subdirs[i]))
            prev_iteration = curr_emirge_iteration
            curr_emirge_iteration = Smurf2Iteration(subdirs[i], prev_smurf2_iteration=prev_iteration)

            # we only need the first iteration and the current one.
            # once we used the previous iteration to initialize the current one, we don't need the data any more.
            if i > 2 and not self.is_debug:
                shutil.rmtree(subdirs[i - 2])

            is_max_iteration = (i == (max_iter - 1))
            if curr_emirge_iteration.do_iteration(is_max_iteration):
                # stable state
                logging.info("DONE SMURF2, iterations = #{}".format(i))
                result_path = curr_emirge_iteration.paths.final_results
                SMURF2_results_path = os.path.join(self.working_dir, "SMURF2_results.csv")
                shutil.copyfile(result_path, SMURF2_results_path)
                shutil.rmtree(subdirs[i])

                # create file structure to support taxonomy script
                sample_name = "sample"

                if not taxa_path:
                    taxa_path = os.path.join(os.getcwd(), "taxonomy/Header_uni_smurf2.csv")

                logging.info("post process: # reads = {}".format(curr_emirge_iteration.n_reads))
                convert_to_smurf_format(SMURF2_results_path,
                                        fasta,
                                        self.working_dir,
                                        sample_name,
                                        curr_emirge_iteration.n_reads,
                                        number_of_region,
                                        taxa_path)

                res_path = self.working_dir + "/"
                if matlab_path is not None:
                    matlab_exe = os.path.join(os.getcwd(), "taxonomy/run_main_taxonomy_single_SMURF2.sh")
                    taxonomy_data = os.path.join(os.getcwd(), "taxonomy/name_calls_ravid_split_genus_split_species.mat")
                    taxonomy_cmd = "{} {} {} {}".format(matlab_exe, matlab_path, res_path, taxonomy_data)
                    logging.info("Running taxonomy command: {}".format(taxonomy_cmd))
                    os.system(taxonomy_cmd)
                    logging.info("Done running taxonomy command: {}".format(taxonomy_cmd))

                if not self.is_debug:
                    for iter_dir in subdirs:
                        if os.path.isdir(iter_dir):
                            shutil.rmtree(iter_dir)
                    os.remove(os.path.join(self.working_dir, "sample_sample_results.mat"))
                break


def main(argv=sys.argv[1:]):
    """
    command line interface to emirge

    """
    parser = OptionParser(USAGE)

    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "These flags are all required to run EMIRGE, and may be supplied in any order.")

    group_reqd.add_option("-1", dest="fastq_reads_1", metavar="reads_1.fastq[.gz]",
                          type="string",
                          help="path to fastq file with \\1 (forward) reads from paired-end."
                               " File must be unzipped for mapper.")
    group_reqd.add_option("-2", dest="fastq_reads_2", metavar="reads_2.fastq",
                          type="string",
                          help="path to fastq file with \\2 (reverse) reads from paired-end run."
                               " File must be unzipped for mapper.")
    group_reqd.add_option("--pr", dest="primers_path", default="primers.csv",
                          type="string",
                          help="path to the primer.csv file."
                               "The file should contain one column for each region, "
                               "the header should be the regions, "
                               "the regions should match the regions indexes in the fasta file.")
    group_reqd.add_option("-f", "--fasta_db",
                          type="string",
                          help="path to fasta files directory of candidate SSU sequences, expected file for each region. "
                               "region index should be include in the fasta file name.")
    group_reqd.add_option("-l", "--read_length",
                          type="int",
                          help="""length of read in input data""")
    group_reqd.add_option("-r", "--regions",
                          type="int",
                          help="""number of regions""")
    parser.add_option_group(group_reqd)

    # THRESHOLDS
    group_thresholds = OptionGroup(parser, "Thresholds parameters",
                                   "Defaults should normally be fine for these options in order to run EMIRGE")
    group_thresholds.add_option("-m", "--min_minor_prob_for_split",
                                type="float", default="0.1",
                                help="minimum probability of second most probable base at a site "
                                     "required in order to call site a variant."
                                     "See also split_threshold. (default: %default)")
    group_thresholds.add_option("-t", "--split_threshold",
                                type="float", default="0.05",
                                help="limit the amount of bases change in one split to this threshold, the split"
                                     "sequence can't be too different from the original sequence (default: %default; valid range: [0.0, 1.0] ) ")
    group_thresholds.add_option("-j", "--join_threshold",
                                type="float", default="0.999",
                                help="If two candidate sequences share >= this fractional identity over their bases with mapped reads, "
                                     "then merge the two sequences into one for the next iteration.  (default: %default; valid range: [0.95, 1.0] ) "
                                     "See also --min_coverage_for_split, --min_minor_prob_for_split (default: %default)")
    group_thresholds.add_option("-c", "--min_coverage_for_split",
                                type="float", default="0.01",
                                help="the split candidate must have estimated coverage >= this threshold"
                                     "See also --join_threshold. (default: %default)")
    group_thresholds.add_option("-s", "--stability_threshold",
                                type="float", default="0.005",
                                help="In stable state the difference between the previous iteration priors and the current priors "
                                     "should be less then this threshold. (default: %default)")
    parser.add_option_group(group_thresholds)

    # OPTIONAL
    group_opt = OptionGroup(parser, "Optional parameters",
                            "Defaults should normally be fine for these options in order to run EMIRGE")
    group_opt.add_option("-n", "--iterations",
                         type="int", default=40,
                         help="""Number of iterations to perform.  It may be necessary to use more iterations for more
                      complex samples (default=%default)""")
    group_opt.add_option("--weight",
                         action="store_true", default=False,
                         help="update the reference db with the number "
                              "of the region that were amplified according to the reads, "
                              "if this option is off, use the weight from the reference DB")
    group_opt.add_option("--debug",
                         action="store_true", default=False,
                         help="Add debug info")
    group_opt.add_option("--split",
                         action="store_true", default=False,
                         help="Allow splitting sequences")
    group_opt.add_option("--taxa",
                         type="string", default=None,
                         help="Store the taxa path")
    group_opt.add_option("--matlab",
                         type="string", default=None,
                         help="Path to the directory where version 9.7 of the MATLAB Runtime is installed on the machine.")

    parser.add_option_group(group_opt)

    # ACTUALLY PARSE ARGS
    (options, args) = parser.parse_args(argv)

    # minimal sanity checking of input
    if len(args) != 1:
        parser.error(
            "DIR is required, and all options except DIR should have a flag associated with them (options without flags: %s)" % args)
    if options.join_threshold < 0.95 or options.join_threshold > 1:
        parser.error(
            "join_threshold must be between [0.95, 1.0].  You supplied %.3f. (see --help)" % options.join_threshold)

    for filename_option_string in ["fastq_reads_1", "fastq_reads_2", "fasta_db"]:
        filename_option = getattr(options, filename_option_string)
        if filename_option is not None:
            if not os.path.exists(filename_option):
                parser.error("file not found for --%s: %s" % (filename_option_string, filename_option))

    working_dir = os.path.abspath(args[0])

    sys.stdout.write(INRTO)

    sys.stdout.write("Command:\n")
    sys.stdout.write(' '.join([__file__] + argv))
    sys.stdout.write('\n\n')
    sys.stdout.flush()

    required = ["fastq_reads_1", "fasta_db", "fastq_reads_2", "read_length"]

    for o in required:
        if getattr(options, o) is None or getattr(options, o) == 0:
            parser.error("--%s is required, but is not specified (try --help)" % (o))
    if options.fastq_reads_2.endswith('.gz'):
        parser.error("Read 2 file cannot be gzipped (see --help)")

    if not os.path.exists(working_dir):
        os.mkdir(working_dir)
    else:
        if len(os.listdir(working_dir)) > 1:  # allow 1 file in case log file is redirected here.
            print >> sys.stderr, os.listdir(working_dir)
            parser.error("Directory not empty: %s\nIt is recommended you run emirge in a new directory each run; "
                         "delete this directory or specify a new one." % working_dir)

    # clean up options to be absolute paths
    for o in ["fastq_reads_1", "fastq_reads_2", "fasta_db"]:
        current_o_value = getattr(options, o)
        if current_o_value is not None:
            setattr(options, o, os.path.abspath(current_o_value))

    emirge_log_name = os.path.join(working_dir, "smurf2.log")

    if options.debug is True:
        define_logger(logging.DEBUG, emirge_log_name)
    else:
        define_logger(file_name=emirge_log_name)

    # CREATE EM OBJECT
    em = EM(working_dir=working_dir,
            is_debug=options.debug)

    # BEGIN ITERATIONS
    em.execute(options.iterations,
               fastq=options.fastq_reads_1,
               fastq_reversed=options.fastq_reads_2,
               primers_path=options.primers_path,
               fasta=options.fasta_db,
               read_length=options.read_length,
               number_of_region=options.regions,
               update_weight_using_the_reads=options.weight,
               max_changed_bases_rate_for_split=options.split_threshold,
               min_coverage_for_split=options.min_coverage_for_split,
               min_minor_prob_for_split=options.min_minor_prob_for_split,
               min_similar_bases_rate_for_merge=options.join_threshold,
               max_priors_diff_for_stability_test=options.stability_threshold,
               allow_split=options.split,
               taxa_path=options.taxa,
               matlab_path=options.matlab)

    return


if __name__ == '__main__':
    main()
