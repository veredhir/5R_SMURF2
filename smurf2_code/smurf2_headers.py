from seqBIn import Base

class HeadersFormat:
    Group_id = 'Reads_group_id'
    Region = 'Region'
    Quals = 'Quals'
    Count = 'Count'
    Read_id = 'Read_id'
    Weight = '#_of_amplified_references-Weight'
    Priors = 'Priors'
    ProbN = 'ProbN'
    Reference_id = 'Reference_id'
    Ref_ids_list = 'Reference_Unique_ids_list'
    Unique_Ref_id = 'Unique_Reference_id'
    Likelihood = 'Likelihood'
    Posterior = 'Posterior'
    Map_weight = 'Map_weight'


class UnqiueRefToRefFormat:
    Unique_id = HeadersFormat.Unique_Ref_id
    Ref_id = HeadersFormat.Reference_id


class PosteriorsFormat:
    Ref_id = HeadersFormat.Unique_Ref_id
    Read_id = HeadersFormat.Read_id
    Likelihood = HeadersFormat.Likelihood
    Posterior = HeadersFormat.Posterior


class ReferenceFormat(object):
    Original_Id = HeadersFormat.Reference_id
    Region = HeadersFormat.Region
    ReferenceIdsList = HeadersFormat.Ref_ids_list
    Bases = Base()
    Ref_Id = HeadersFormat.Unique_Ref_id
    Group_id = HeadersFormat.Group_id
    temp_all = [Original_Id, Region] + Bases.all
    header = [Region, Ref_Id] + Bases.all

    def __init__(self):
        self.dtypes = {self.Original_Id: int, self.Region: int}
        self.dtypes.update(self.Bases.dtypes)
        self.ix = {}
        for key in self.Bases.ix:
            self.ix[key] = self.Bases.ix.get(key) + 2

    def get_ref_dict(self,title, region_ix, base_dict):
        ref_dict = {self.Original_Id: title, self.Region: region_ix}
        ref_dict.update(base_dict)
        return ref_dict

    @staticmethod
    def get_format(seq_id, region, bases_list):
        return [seq_id, region] + bases_list


class ReadsFullDataFormat:
    Id = HeadersFormat.Read_id
    Region = HeadersFormat.Region
    Quals = HeadersFormat.Quals
    Count = HeadersFormat.Count
    Group_id = HeadersFormat.Group_id
    Bases = Base()
    all = [Id, Region, Quals] + Bases.all
    unique_reads_all = [Id, Region] + Bases.all

    def __init__(self):
        self.dtypes = {self.Id: int}
        self.dtypes.update(self.Bases.dtypes)
        self.ix = {}
        for key in self.Bases.ix:
            self.ix[key] = self.Bases.ix.get(key) + 1

    @staticmethod
    def get_format(read_id, region, bases_list, quals):
        return [read_id, region, quals] + bases_list

    def get_dict(self, read_id, region, base_dict, quals):
        read_dict = {self.Id: read_id, self.Region: region, self.Quals: quals}
        read_dict.update(base_dict)
        return read_dict


class ReadsFormat:
    Id = HeadersFormat.Read_id
    Quals = HeadersFormat.Quals
    Group_id = HeadersFormat.Group_id


#maps GROUP of reads to reference ids
class MappingForamt(object):
    Group_id = HeadersFormat.Group_id
    Region = HeadersFormat.Region
    Bases = Base()
    Count = HeadersFormat.Count
    Ref_id = HeadersFormat.Unique_Ref_id
    Map_weight = HeadersFormat.Map_weight
    full_header = [Group_id, Region, Count, Ref_id, Map_weight] + Bases.all


# each reference sequence state
class CurrentStateFormat(object):
    Reference_id = HeadersFormat.Unique_Ref_id
    Region = HeadersFormat.Region
    Weight = HeadersFormat.Weight
    Priors = HeadersFormat.Priors
    Bases = Base()
    ProbN = HeadersFormat.ProbN