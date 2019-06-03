import pandas as pd
from qiime2 import Artifact

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