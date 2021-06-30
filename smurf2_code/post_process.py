
from Bio import SeqIO
import numpy as np
from smurf2_headers import *
from smurf2_utills import *
import os
# import seaborn as sns
# import mat4py
import scipy.io as sio

# import matplotlib
# matplotlib.use('Agg')
#
# import matplotlib.pyplot as plt
import pandas as pd
# from pandas.tools.plotting import table
from optparse import OptionParser, OptionGroup
import glob
import sys
# from oct2py import octave as oc

SIMILAR = 'similar'
CONTAIN = 'contain'
RECALL = 'Recall'
PRECISION = 'Precision'

VALIDATION_THRESHOLD_THRESHOLD = 0.00001
EXPECTED_RES_FILE_NAME = "expected_res.csv"
ACTUAL_RES_FILE_NAME = "emirge_smurf_WFalseSTrue.csv"


class Header():
    ref_id = 'ref'
    prior = 'prior'
    sequence = 'sequence'
    region = 'region'
    weight = 'weight'
    is_changed = 'is_changed'
    new_id = 'new_id'


def write_bacterium_to_fasta(new_bacteria_path,
                             new_bacteria):
    """
    :param new_bacteria_path:
    :param new_bacteria: dictionary, keys are the ids and values are the sequences
    :return:
    """
    with open(new_bacteria_path, 'w') as mock_fasta1:
        for id in new_bacteria.keys():
                mock_fasta1.write(">{}\n".format(id))
                mock_fasta1.write(new_bacteria[id] + '\n')
                mock_fasta1.write('\n')

def validate_priors(df, threshold=VALIDATION_THRESHOLD_THRESHOLD):
    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])
    df.prior = df.prior.apply(lambda r: r / sum_prior)

    sum_prior = sum(df.drop_duplicates(Header.ref_id)[Header.prior])

    # logging.debug("sum of priors is {}".format(sum_prior))
    if abs(sum_prior - 1) > threshold:
        raise Exception("sum of prior is not 1")


def write_new_bacterium_to_fasta(new_bacteria_df,
                                 regions,
                                 fasta_dir):
    regions = range(1, regions+1)
    for region in regions:
        if region not in new_bacteria_df.columns:
            continue
        fasta_name = "new_bacteria_{}.fasta".format(region)
        fasta_file = os.path.join(os.path.join(fasta_dir, fasta_name))
        with open(fasta_file, 'w') as mock_fasta1:
            for _, row in new_bacteria_df.iterrows():
                id = row[Header.new_id]
                seq = row[region]
                if seq != '':
                    mock_fasta1.write(">{}\n".format(id))
                    mock_fasta1.write(seq + '\n')
                    mock_fasta1.write('\n')


class Found_bacteria():
    def __init__(self):
        self.frequency=[]
        self.assigned_reads=[]
    def get(self):
        return {'frequency': self.frequency,
                'assigned_reads': self.assigned_reads}



class BacteriaMetaGroup():
    def __init__(self):
        self.db_ind=[]
        self.comb_vector=[]
        self.new_ind=[]
        self.is_out_of_db=[]

    def append(self, db_ind, comb_vector, new_ind, is_out_of_db):
        self.db_ind.append(db_ind)
        self.comb_vector.append(comb_vector)
        self.new_ind.append(new_ind)
        self.is_out_of_db.append(is_out_of_db)

    def get(self):
        return {'comb_vec': np.array(self.comb_vector),
               'db_ind': np.array(self.db_ind),
               'new_ind': self.new_ind,
               'is_out_of_db':self.is_out_of_db}



def read_db_ind_map(path_to_db_ind):
    db_map = pd.read_csv(path_to_db_ind, index_col=None, header=None)
    db_map['fasta_id'] = db_map[0].apply(lambda x: str(int(x)))
    db_map['index'] = range(1, len(db_map)+1)
    res = db_map.groupby('fasta_id')['index'].apply(list).reset_index(name='indices')
    res = res[['fasta_id','indices' ]].drop_duplicates('fasta_id')
    return res


def fasta_id_to_indices(fasta_ids, db_map):
    for fasta_id in fasta_ids:
        res = db_map[db_map['fasta_id'] == str(int(fasta_id))]['indices']
        if len(res) < 1:
            print("ERROR - size of db map for id [{}] = {}".format(fasta_id, len(res)))
            continue
        return [int(r) for r in res.iloc[0]]
    return []


def to_smurf_format(df, overall_reads, num_regions, mat_dest, taxa_path):
    """
    found_bacteria: frequency, assigned_reads
    bateriaMetaGroups: db_ind, comb_vec, new_ind, is_out_of_db
    :return:
    """
    bacteriaMetaGroups=BacteriaMetaGroup()
    found_bacteria=Found_bacteria()
    assigned_reads_frac = []
    db_map = read_db_ind_map(taxa_path)

    ref_groups = df.groupby(Header.ref_id)
    for ref_id, ref_df in ref_groups:
        regions = 5*[0]
        regions_list = ref_df[Header.region].unique()
        for region in regions_list:
            regions[int(region)-1] = 1
        prior = ref_df[Header.prior].iloc[0]
        new_id = ref_df[Header.new_id].iloc[0]
        is_out_of_db = ref_df[Header.is_changed].iloc[0]
        found_bacteria.frequency.append(prior)
        assigned_reads_frac.append(prior*sum(regions)/num_regions)
        fasta_id = ref_df['Reference_id'].tolist()
        garys_ids = fasta_id_to_indices(fasta_id, db_map)
        bacteriaMetaGroups.append(garys_ids, regions, new_id, is_out_of_db)

    norm_factor = sum(assigned_reads_frac)
    found_bacteria.assigned_reads = [overall_reads*reads/norm_factor for reads in assigned_reads_frac]
    data = {'found_bacteria': found_bacteria.get(), 'bactMetaGroups':bacteriaMetaGroups.get()}
    sio.savemat(mat_dest, data)
    print("save mat to {}".format(mat_dest))
    #
    # oc.unpackStruct(mat_dest)
    # oc.save(mat_dest)


def are_similar(seq1, db_seq):
    length = len(seq1)/2
    seq2 = db_seq[:length] + db_seq[-1*length:]
    if seq2 == seq1:
        return True
    if 'N' in seq2:
        N_pos = [pos for pos, char in enumerate(seq2) if char == 'N']
        for pos in N_pos:
            seq1 = seq1[0:pos] + 'N' + seq1[pos + 1:]
        if seq2 == seq1:
            return True
    return False


def count_changes(seq1, db_seq):
    length = len(seq1) / 2
    seq2 = db_seq[:length] + db_seq[-1*length:]
    counter=0
    for s1, s2 in zip(seq1, seq2):
        if s1 == s2:
            continue
        if s1 == 'N' or s2 == 'N':
            continue
        counter+=1
    return counter


def convert_to_smurf_format(path, fasta_dir, output_dir, sample_name, overall_reads, num_regions, taxa_path):
    """
    smurf2 to smurf format
    :param output_dir:
    :param path: path to 'final_results.csv produced by smurf2.py
    :return: df hold the final results
    """
    df = pd.read_csv(path, index_col=False)
    df = df.rename(columns = {'Sequence': Header.sequence,
                              HeadersFormat.Region: Header.region,
                              HeadersFormat.Priors: Header.prior,
                              'Unique_Reference_id': Header.ref_id})
    df['Header'] = df[Header.ref_id].apply(lambda r: round(r, 2))
    new_id_ix=1

    regions = df[Header.region].unique()
    for region in regions:
        df[region] = (df[Header.region] == region) * df[Header.sequence]

    full_df = df.drop_duplicates([Header.ref_id, Header.region])
    full_df.loc[:, Header.is_changed] = False
    full_df.loc[:, Header.new_id] = 0

    fasta_files = get_fasta_files(fasta_dir, num_regions)
    for fasta_file in fasta_files:
        records = SeqIO.index(fasta_file.path, "fasta")
        curr_region = fasta_file.region
        for index, row in full_df.iterrows():
            if curr_region not in row.index or row[int(curr_region)]=='':
                continue
            id = row['Reference_id']
            if float(id)-int(id) != 0:
                id = str(int(id))
            try:
                seq = records[id].seq
                db_seq = seq.__str__()
            except Exception as ex:
                logging.error("error extracting: ex = {}, id={}, region = {}, prior={}".format(ex, id, curr_region, row[Header.prior]))
                db_seq = ''
            curr_seq = row[curr_region]
            if not are_similar(curr_seq, db_seq):
                full_df.loc[index, Header.is_changed] = True
                full_df.loc[index, Header.new_id] = new_id_ix
                new_id_ix+=1
        records.close()

    validate_priors(df)
    write_new_bacterium_to_fasta(full_df[full_df[Header.is_changed] == True], num_regions, output_dir)

    df = df.merge(full_df[[Header.is_changed, Header.ref_id, Header.new_id]], on=Header.ref_id)
    results_path = os.path.join(output_dir, sample_name + "_sample_results.mat")
    to_smurf_format(df, overall_reads, num_regions, results_path, taxa_path)

    return results_path


class FastaFile(object):
    def __init__(self, directory, file_name, total_num_of_regions=None, region=None):
        self.path = os.path.join(directory, file_name)
        self.file_name = file_name
        self.region = region
        self.total_num_of_regions = total_num_of_regions

    def initialize(self):
        if self.region == None:
            if self.total_num_of_regions is None:
                raise Exception("Missing 'total number of regions'")
            self.region = self.get_region_from_file_name()

    def get_region_from_file_name(self):
        region = -1
        for i in range(1, self.total_num_of_regions + 1):
            if str(i) in self.file_name:
                region = i
        return region


def get_fasta_files(fasta_dir, total_num_of_regions):
    fasta_files = []
    for f_name in os.listdir(fasta_dir):
        if f_name.endswith('fasta'):
            fasta_file = FastaFile(fasta_dir, f_name, total_num_of_regions)
            fasta_file.initialize()
            fasta_files.append(fasta_file)
    return fasta_files


def main(argv = sys.argv[1:]):
    """
    command line interface to emirge

    """
    parser = OptionParser("Convert SMURF2 results format to SMURF results format")

    # REQUIRED
    group_reqd = OptionGroup(parser, "Required flags",
                             "")
    group_reqd.add_option("-r", "--smurf2_results",
                      type="string",
                      help="path to smurf2 results csv")
    group_reqd.add_option("-f", "--fasta",
                          type="string",
                          help="path to regions fasta directory")
    group_reqd.add_option("-o", "--output",
                          type="string",
                          help="path to output directory")
    group_reqd.add_option("-s", "--sample_name",
                          type="string",
                          help="sample name")
    group_reqd.add_option("-l", "--overall_reads",
                      type="int",
                      help="""overall mapped reads""")
    group_reqd.add_option("-t", "--taxa_path",
                          type="string", default="/home/vered/EMIRGE/data/reference_db/Header_uni_forVered.csv",
                          help="taxa_and_head file path")

    parser.add_option_group(group_reqd)
    (options, args) = parser.parse_args(argv)

    path=options.smurf2_results
    fasta_dir=options.fasta
    output_dir=options.output
    sample_name=options.sample_name
    overall_reads=options.overall_reads

    convert_to_smurf_format(path, fasta_dir, output_dir, sample_name, overall_reads, 5, options.taxa_path)




if __name__ == "__main__":
    main()



