import errno
import os
import pandas as pd
from qiime2 import Artifact
from qiime2 import Metadata
from qiime2.plugins import emperor
from qiime2.plugins.diversity.pipelines import core_metrics_phylogenetic
from qiime2.plugins.diversity.visualizers import alpha_group_significance


def taxon2fasta(taxonomy, sequences, taxon, path):
    '''
    taxonomy is an artifact of type FeatureData[Taxonomy]
    sequences is an artifact of type FeatureData[Sequence]
    taxon is the annotated OTU we are interested in. input string
    path is where to export the fasta files. input string
    '''
    # convert FeatureData[Taxonomy] to pandas dataframe
    df_taxon = taxonomy.view(pd.DataFrame)
    
    # filter ASV that were annotated to 'taxon'
    df_taxon = df_taxon.loc[(df_taxon.loc[:,'Taxon'] == taxon)]
    
    # convert FeatureData[Sequence] to pandas series
    ser = sequences.view(pd.Series)
    
    # filter seqs that were annotated to 'taxon'
    ser_taxon = ser[df_taxon.index]
    
    # covert filtered seqs to artifact
    taxon_seq = Artifact.import_data('FeatureData[Sequence]', ser_taxon)
    
    # export fasta files to given path
    taxon_seq.export_data(path) 
    

def mkdir_(name):
    try:
        os.mkdir(name)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise
        pass    
    
    
def core_metrics2qzv(table, phylogeny, metadata, name, n_jobs=1):
    '''
    table is an artifact of type FeatureTable[Frequency]
    phylogeny is an artifact of type Phylogeny[Rooted]
    metadata is the metadata imported by qiime2.Metadata.load()
    name is the name of a folder where to export .qzv files. input string
    '''
    # TODO
    # input type check
    # phylogenetic metrics / non-phylogenetic metrics

    df_table = table.view(pd.DataFrame)

    # core metrics calculation
    sampling_depth = int(min(df_table.sum(axis=1)))
    core_metrics = core_metrics_phylogenetic(table = table,
                                             phylogeny = phylogeny,
                                             n_jobs = n_jobs,
                                             sampling_depth = sampling_depth,
                                             metadata = metadata)


    mkdir_(name)
    mkdir_(name + '/alpha_diversity')
    mkdir_(name + '/PCOA')

    # alpha diversity
    alpha_metrics = ['faith_pd_vector', 'evenness_vector', 'observed_otus_vector', 'shannon_vector' ]
    for attr in alpha_metrics:
        significance = alpha_group_significance(getattr(core_metrics, attr),
                                                metadata = metadata)
        significance.visualization.save(project_name + '/alpha_diversity/' + attr[:-7])

        alpha_div = (getattr(core_metrics, attr).view(Metadata)).to_dataframe()
        df = pd.concat([df, alpha_div], axis=1) 

    df.to_csv(name + '/alpha_diversity/alpha_diversity.tsv', sep='\t')

   
    # beta diversity PCOA
    pcoa_results = ['bray_curtis_pcoa_results', 'unweighted_unifrac_pcoa_results',\
                    'weighted_unifrac_pcoa_results', 'jaccard_pcoa_results' ]
    for attr in pcoa_results:
        pcoa_plot = emperor.actions.plot(getattr(core_metrics, attr),
                                         metadata)
        # export emperor PCoA plot
        pcoa_plot.visualization.save(name + '/PCOA/' +  attr[:-8])
        # export raw data for PCoA plot 
        pcoa_result = (getattr(core_metrics, attr)).view(skbio.OrdinationResults)
        coordinate = pcoa_result.samples
        coordinate.loc['proportion_explained'] = pcoa_result.proportion_explained
        coordinate.to_csv(name + '/PCOA/' + attr[:-8] + '.tsv', sep='\t')

def core_ASV(table, abu, prev):
    '''
    table is an artifact of type FeatureTable[Frequency]
    abu is the threshold of abundance.  Abundance of ASV lower than this will be filtered out. [0-1]
    prev is the threshold of prevalence. [0-1]
    '''
    df = table.view(pd.DataFrame) # convert artifact to pd.DataFrame
    df = df.div(df.sum(axis = 1), axis = 0) # covert counts of ASV to abundance
    
    # get ASVs with prevelance higher than threshold. When calculate the prevalence of given ASV,
    # only those abundance greater than threshold will be considered.
    core = df.loc[:,(df>abu).sum()>(len(df)*prev)]
    
    return core

def table2lefse(table, metadata, delimiter, c, sub_c=False):
    '''
    table is an artifact of type FeatureTable[Frequency]
    metadata is the metadata imported by qiime2.Metadata.load()
    delimiter is the delimiter in taxon to separate different levels
    c: class in lefse format
    sub_c: sub class in lefse format
    '''
    df_table = table.view(pd.DataFrame)
    col_sep = df_table.columns.str.split(delimiter).str.join('|')
    df_table.columns = col_sep
    df_table = df_table.div(df_table.sum(axis = 1), axis = 0)
    
    df_meta = (metadata.get_column(c)).to_dataframe()
    df_meta = df_meta.loc[df_table.index]

    # class
    loc=0
    df_table.insert(loc=loc, column=c, value=df_meta[c])
    loc += 1
    
    # sub-class
    if sub_c != False:
        df_meta = (metadata.get_column(sub_c)).to_dataframe()
        df_meta = df_meta.loc[df_table.index]
        df_table.insert(loc=loc, column=sub_c, value=df_meta[sub_c])
        loc += 1
    # ID
    df_table.insert(loc=loc, column='Sample_ID', value=df_table.index)
    
    df_table = df_table.transpose()  
    df_table.to_csv('LEfSe_input.tsv', sep='\t', index=True, header=False)
