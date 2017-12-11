import glob
import os
import io
import zipfile
import pandas as pd
from pandas.api.types import is_numeric_dtype
from scipy.stats import zscore
import requests


# from pathlib import Path

# IDs to download the human brain data from Allen API
well_known_file_IDs = [178238387, 178238373, 178238359, 178238316, 178238266,
                      178236545]

# raw data pathlib
raw_data_location = './data/raw'


def download_HBA_adult_human_data(file_id, path):
    url = 'http://api.brain-map.org/api/v2/'
    url += 'well_known_file_download/{}'.format(file_id)
    r = requests.get(url)
    z = zipfile.ZipFile(io.BytesIO(r.content))
    return z.extractall(path)


def read_expression_file(file_name):
    expression_df = pd.read_csv(file_name, index_col=0, header=None)
    expression_df.index.rename('probe_id', inplace=True)
    return expression_df


def read_samples_file(samples_file):
    sample_df = pd.read_csv(samples_file)
    sample_df.set_index(sample_df.index+1, inplace=True)
    sample_df.index.rename('sample_id', inplace=True)
    return sample_df


def get_probes_data(probes_file, probes_strategy='default'):
    strats = ['default', 'reannotator', 'qc_filter', 'qc_scale']
    assert probes_strategy in strats

    # depending on strategy, may merge in diff tables to update probes info
    probes_df = pd.read_csv(probes_file)

    # rename columns for consistency between adult and fetal brain datasets
    if 'probeset_name' in probes_df.columns:
        probes_df.rename(columns={'probeset_name': 'probe_name',
                                  'probeset_id': 'probe_id'}, inplace=True)
    cols = ['probe_id', 'probe_name', 'gene_symbol']
    probes_df = probes_df.loc[:, cols]

    if probes_strategy == 'reannotator':
        reannotations = get_probe_reannotations('./data/raw/gene_symbol_annotations/AllenInstitute_custom_Agilent_Array.txt')
        # drop the original gene_symbol column
        probes_df.drop('gene_symbol', axis=1, inplace=True)
        # merge in the reannotated gene_symbols
        probes_df = probes_df.merge(reannotations, on='probe_name')

    elif probes_strategy in ['qc_filter', 'qc_scale']:
        # qc_filt = get_probe_qc_filter('/Users/derek_howard/Downloads/12864_2013_7016_MOESM8_ESM.xlsx')
        qc_filt = get_probe_qc_filter('./data/raw/gene_symbol_annotations/12864_2013_7016_MOESM8_ESM.xlsx')
        probes_df = probes_df.merge(qc_filt, left_on='probe_name', right_on='probe')
        probes_df = probes_df[probes_df.qc_filter == True]
        assert is_numeric_dtype(probes_df.m)
        assert is_numeric_dtype(probes_df.b)

        print('After getting probes_df which merged qc data, shape is {}'.format(probes_df.shape))

    probes_df.set_index('probe_id', inplace=True)

    return probes_df


def get_probe_reannotations(re_annotations_file):
    # pre-processing function to prepare the probe reannotations file to be merge with probes df
    re_annotations = pd.read_table(re_annotations_file,
                                   usecols=['#PROBE_ID', 'Gene_symbol'])
    re_annotations.rename(columns={'#PROBE_ID': 'probe_name'}, inplace=True)
    re_annotations.dropna(inplace=True)
    re_annotations.set_index('probe_name', inplace=True)
    # split gene_symbols which have multiple genes associated with a single
    # probe_name creates a new row for each of the gene_symbols
    re_annotations = (re_annotations.Gene_symbol.str.split(';', expand=True)
                                    .stack()
                                    .reset_index()
                                    .drop('level_1', axis=1)
                                    .rename(columns={0: 'gene_symbol'}))
                                    #.rename(columns={0: 'reannotated_gene_symbol'}))
    return re_annotations


def get_probe_qc_filter(qc_file):
    # pre-processes the Allen qc filter file
    df = pd.read_excel(qc_file)
    df = df.iloc[:, [0, 1, 2, 3, 7]].drop(0)  # .rename(columns=[col_names])
    col_names = ['probe', 'gene_symbol', 'm', 'b', 'qc_filter']
    df.columns = col_names
    df[['m', 'b']] = df[['m', 'b']].apply(pd.to_numeric)
    return df.drop('gene_symbol', axis=1)


def scale_expression_vals(exp_df, probes_qc):
    print('Expression df shape before scaling: {}'.format(exp_df.shape))
    df = exp_df.merge(probes_qc.loc[:, ['m', 'b']], left_index=True, right_index=True)
    df.m = df.m.astype(float)
    df.b = df.b.astype(float)
    try:
        df = df.apply(lambda x: df.m * x + df.b)
        #df = df.apply(lambda x: np.log2(df.m * 2**x + df.b), axis=1)
    except Exception as e:
        print('Exception occurred on row when applying scale fxn, skip row')
        print(e)
    df = df.drop(['m', 'b'], axis=1)
    # assert
    df.dropna(inplace=True)
    print('Expression df shape after scaling: {}'.format(df.shape))
    return df


def get_donor_data(donor_file_list, probes_strategy):
    # to work with both fetal and adult metadata
    probe_file_strings = ['Probes', 'rows_meta']
    samples_file_strings = ['Sample', 'columns_meta']
    expression_file_strings = ['Expression', 'expression']

    for file in donor_file_list:
        if any(string in file for string in probe_file_strings):
            probes_df = get_probes_data(file, probes_strategy=probes_strategy)
        elif any(string in file for string in samples_file_strings):
            samples_df = read_samples_file(file)
        elif any(string in file for string in expression_file_strings):
            exp_df = read_expression_file(file)
        else:
            continue

    return exp_df, samples_df, probes_df


def get_mean_expression_by_brain_area(exp_df, samples_df):
    assert(exp_df.T.shape[0] == samples_df.shape[0])

    # merge in metadata (brain area of sample)
    annotated_df = exp_df.T.merge(samples_df[['structure_name']],
                                  left_index=True, right_index=True)

    # get mean expression level for samples within a brain area
    expression_by_structure = annotated_df.groupby('structure_name').mean()
    expression_by_structure.T.index.rename('gene_symbol', inplace=True)

    return expression_by_structure.T


def get_exp_by_genes(exp_df, probes_df):
    """
    input is exp_df and probes metadata
    output is exp_df grouped by gene_symbols and averaged
    """
    annotated_df = exp_df.merge(probes_df[['gene_symbol']],
                                left_index=True, right_index=True)

    # remove probes that don't have an associated gene symbol
    #annotated_df = annotated_df[~((annotated_df.gene_symbol.str.startswith('CUST')) | (annotated_df.gene_symbol.str.startswith('A_')))]

    exp_by_genes = (annotated_df.groupby('gene_symbol')
                                .mean())
                                # .drop('na'))
    return exp_by_genes


def strip_left_right(structure_name):
    brain_area_fragments = structure_name.split(',')
    clean_fragments = []

    for frag in brain_area_fragments:
        if frag.strip() not in ['left', 'right']:
            clean_fragments.append(frag)

    clean_structure_name = ','.join(clean_fragments)
    return clean_structure_name


def get_single_donor_tidy_df(exp_df, samples_df, probes_df, donor_id, probes_strategy):
    # remove left/right from brain structure_names
    samples_df.structure_name = samples_df.structure_name.apply(strip_left_right)

    if probes_strategy == 'qc_scale':
        exp_df = scale_expression_vals(exp_df, probes_df)

        print('size of df after scaling: {}'.format(exp_df.shape))
    # merge in probes metadata and get expression by gene_symbol
    expression_by_genes = get_exp_by_genes(exp_df, probes_df)
    print('size of df after merging probe info, grouping by gene: {}'.format(expression_by_genes.shape))
    # merge in sample metadata and get expression by brain area
    exp_brain_area_by_genes = get_mean_expression_by_brain_area(expression_by_genes, samples_df)

    # change these processing steps for debugging SSMIs
    ranked_exp_by_area = exp_brain_area_by_genes.rank(ascending=True)
    zscored_exp_by_area = ranked_exp_by_area.apply(zscore, axis=1)

    melted = pd.melt(zscored_exp_by_area.reset_index(),
    #melted = pd.melt(exp_brain_area_by_genes.reset_index(),
                     id_vars='gene_symbol',
                     var_name='brain_area')
    melted['donor_id'] = donor_id

    return melted


def generate_HBA_dataset(dataset, probes_strategy):
    assert dataset in ['adult', 'fetal']
    if dataset == 'adult':
        print('--- Generating adult human HBA dataset ---')
        data_path = './data/raw/allen_HBA'
    elif dataset == 'fetal':
        print('--- Generating fetal human HBA dataset ---')
        data_path = './data/raw/allen_human_fetal_brain'

    hba_donor_folders = glob.glob(os.path.join(data_path, '*'))
    all_donors = []
    for donor_folder in hba_donor_folders:
        donor_id = donor_folder.split('/')[-1].split('_')[-1]
        donor_files = glob.glob(os.path.join(donor_folder, '*'))
        expression, samples, probes = get_donor_data(donor_files,
                                                     probes_strategy)

        tidy_donor = get_single_donor_tidy_df(expression,
                                              samples,
                                              probes,
                                              donor_id=donor_id,
                                              probes_strategy=probes_strategy)
        all_donors.append(tidy_donor)

    all_donors_long = pd.concat(all_donors)

    structure_genes_exp_matrix = pd.pivot_table(all_donors_long,
                                                values='value',
                                                index='gene_symbol',
                                                columns='brain_area')
    return structure_genes_exp_matrix


def get_dataset(dataset, probes_strategy):
    # to reduce code duplication between get_HBA_dataset and
    # get_fetal_HBA_dataset
    strats = ['reannotator', 'qc_filter', 'qc_scale', 'default']
    datasets = ['adult', 'fetal']
    assert probes_strategy in strats
    assert dataset in datasets

    filename = '{0}_brainarea_vs_genes_exp_{1}.tsv'.format(dataset,
                                                           probes_strategy)
    HBA_data_out_path = os.path.join('data', 'processed', filename)

    if os.path.exists(HBA_data_out_path):
        print(('Processed HBA brain dataset found locally. ' +
               'Loading from {}'.format(HBA_data_out_path)))
        structure_genes_exp_matrix = pd.read_csv(HBA_data_out_path,
                                                 index_col='gene_symbol',
                                                 sep='\t')

    else:
        os.makedirs('./data/processed', exist_ok=True)
        structure_genes_exp_matrix = generate_HBA_dataset(dataset,
                                                          probes_strategy)
        print('-- Writing data to ' + HBA_data_out_path + ' -- ')
        structure_genes_exp_matrix.to_csv(HBA_data_out_path, sep='\t')
    return structure_genes_exp_matrix
