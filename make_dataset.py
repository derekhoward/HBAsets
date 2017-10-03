# %load allen_HBA.py
import pandas as pd
import glob
import os
from scipy.stats import zscore
import requests
import zipfile
import io

# IDs to download the human brain data from Allen API
well_know_file_IDS = [178238387, 178238373, 178238359, 178238316, 178238266,
                      178236545]


def download_HBA_adult_human_data(file_id):
    url = 'http://api.brain-map.org/api/v2/'
    url += 'well_known_file_download/{}'.format(file_id)
    r = requests.get(url)
    z = zipfile.ZipFile(io.BytesIO(r.content))
    # should make sure to find right location to save extracted data
    return z.extractall()


def read_expression_file(file_name):
    expression_df = pd.read_csv(file_name, index_col=0, header=None)
    expression_df.index.rename('probe_id', inplace=True)
    return expression_df


def read_samples_file(samples_file):
    sample_df = pd.read_csv(samples_file)
    sample_df.set_index(sample_df.index+1, inplace=True)
    sample_df.index.rename('sample_id', inplace=True)
    return sample_df


def get_probes_data(probes_file, reannotate=False):
    probes_df = pd.read_csv(probes_file)

    # rename columns for consistency between adult and fetal brain datasets
    if 'probeset_name' in probes_df.columns:
        probes_df.rename(columns={'probeset_name': 'probe_name',
                                  'probeset_id': 'probe_id'}, inplace=True)

    if reannotate:
        reannotations = get_probe_reannotations('./data/raw/ \
            gene_symbol_annotations/AllenInstitute_custom_Agilent_Array.txt')
        probes_df = merge_reannotations(probes_df, reannotations)

    probes_df.set_index('probe_id', inplace=True)

    return probes_df


def get_probe_reannotations(re_annotations_file):
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
                                    .rename(columns={0: 'reannotated_gene_symbol'}))
    return re_annotations


def merge_reannotations(probes_df, re_annotations_df):
    merged_probes = probes_df.merge(re_annotations_df, on='probe_name')
    return merged_probes


def get_donor_data(donor_file_list, reannotate=False):
    probe_file_strings = ['Probes', 'rows_meta']
    samples_file_strings = ['Sample', 'columns_meta']
    expression_file_strings = ['Expression', 'expression']

    for file in donor_file_list:
        if any(string in file for string in probe_file_strings):
            probes_df = get_probes_data(file, reannotate=reannotate)
        elif any(string in file for string in samples_file_strings):
            samples_df = read_samples_file(file)
        elif any(string in file for string in expression_file_strings):
            exp_df = read_expression_file(file)
        else:
            continue

    return exp_df, samples_df, probes_df


def get_mean_expression_by_brain_area(exp_df, samples_df):
    assert(exp_df.T.shape[0] == samples_df.shape[0])

    # merge brain area info from which expression_data was taken
    annotated_df = exp_df.T.merge(samples_df[['structure_name']],
                                  left_index=True, right_index=True)

    # get mean expression level for samples within a brain area
    expression_by_structure = annotated_df.groupby('structure_name').mean()
    expression_by_structure.T.index.rename('probe_id', inplace=True)

    return expression_by_structure.T


def add_genesymbols2expression(exp_df, probes_df, reannotate=False):
    # assert(exp_df.shape[0] == probes_df.shape[0])
    if reannotate:
        annotated_df = exp_df.merge(probes_df[['entrez_id',
                                               'reannotated_gene_symbol']],
                                    left_index=True, right_index=True)
        annotated_df.rename(columns={'reannotated_gene_symbol': 'gene_symbol'},
                            inplace=True)

    else:
        annotated_df = exp_df.merge(probes_df[['gene_symbol', 'entrez_id']],
                                    left_index=True, right_index=True)

    # remove probes that don't have an associated gene symbol
    # annotated_df = annotated_df[~((annotated_df.gene_symbol.str.startswith('CUST')) | (annotated_df.gene_symbol.str.startswith('A_')))]

    return annotated_df


def get_expression_by_struct_genesymbol(exp_by_structure_df):
    """
    input is a df of (rows: probe_ids) vs (columns: brain_areas + gene_ids), where values are expression levels per probe_id
    output is a df where (rows: gene_symbols) vs (columns: brain_areas), where values are mean expression levels of gene_symbol in a brain_area
    """
    exp_struct_genes = (exp_by_structure_df.drop(['entrez_id'], axis=1)
                                           .groupby('gene_symbol')
                                           .mean())
                                           # .drop('na'))
    return exp_struct_genes


def strip_left_right_from_brain_areas(df):
    stripped_col_names = []
    for col in df:
        brain_area_fragments = col.split(',')
        clean_fragments = []
        for fragment in brain_area_fragments:
            if fragment.strip() not in ['left', 'right']:
                clean_fragments.append(fragment)
        clean_col_name = ','.join(clean_fragments)
        stripped_col_names.append(clean_col_name)
    return set(stripped_col_names)


def merge_left_right_brain_areas(df):
    # input is a dataframe with columns that have separate columns for left/right brain areas
    # creates a list where each element is a series where the values are averaged expression values
    # from left/right brain areas
    # output is a dataframe which concats the list of series
    stripped_col_names = strip_left_right_from_brain_areas(df)
    merged_cols = []
    for stripped_col in stripped_col_names:
        stripped_col_fragments = stripped_col.split(',')
        similar_df = df.copy()
        for fragment in stripped_col_fragments:
            similar_df = similar_df.filter(like=fragment)
        mean_col = similar_df.mean(axis=1).rename(stripped_col)
        merged_cols.append(mean_col)

    output_df = pd.concat(merged_cols, axis=1)

    return output_df.sort_index(axis=1)


def get_single_donor_tidy_df(exp_df, samples_df, probes_df, donor_id, reannotate=False):
    expression_by_structure = get_mean_expression_by_brain_area(exp_df, samples_df)
    expression_by_structure_noLR = merge_left_right_brain_areas(expression_by_structure)
    annotated_expression_by_structure = add_genesymbols2expression(expression_by_structure_noLR, probes_df, reannotate=reannotate)
    # the following function will remove an 'na' row from the df
    gene_symbol_expression_by_structure = get_expression_by_struct_genesymbol(annotated_expression_by_structure)
    ranked_exp_by_area = gene_symbol_expression_by_structure.rank(ascending=True)
    zscored_exp_by_area = ranked_exp_by_area.apply(zscore, axis=1)
    melted = pd.melt(zscored_exp_by_area.reset_index(),
                     id_vars='gene_symbol',
                     var_name='brain_area')
    melted['donor_id'] = donor_id

    return melted


def make_HBA_dataset(reannotate=False):
    if reannotate:
        HBA_data_out_path = './data/processed/brainarea_vs_genes_exp_w_reannotations.tsv'
    else:
        HBA_data_out_path = './data/processed/brainarea_vs_genes_exp.tsv'

    if os.path.exists(HBA_data_out_path):
        print('Processed HBA brain dataset found locally. Loading from {}'.format(HBA_data_out_path))
        structure_genes_exp_matrix = pd.read_csv(HBA_data_out_path,
                                                 index_col='gene_symbol',
                                                 sep='\t')

    else:
        data_path = './data/raw/allen_HBA'
        hba_donor_folders = glob.glob(os.path.join(data_path, '*'))

        all_donors = []

        for donor_folder in hba_donor_folders:
            donor_id = donor_folder.split('/')[-1].split('_')[-1]
            donor_files = glob.glob(os.path.join(donor_folder, '*'))
            expression, samples, probes = get_donor_data(donor_files,
                                                         reannotate=reannotate)
            tidy_donor = get_single_donor_tidy_df(expression,
                                                  samples,
                                                  probes,
                                                  donor_id=donor_id,
                                                  reannotate=reannotate)
            all_donors.append(tidy_donor)

        all_donors_long = pd.concat(all_donors)

        structure_genes_exp_matrix = pd.pivot_table(all_donors_long,
                                                    values='value',
                                                    index='gene_symbol',
                                                    columns='brain_area')
        structure_genes_exp_matrix.to_csv(HBA_data_out_path, sep='\t')
    return structure_genes_exp_matrix


def make_fetal_HBA_dataset(reannotate=False):
    if reannotate:
        fetal_data_out_path = './data/processed/fetal_brainarea_vs_genes_exp_w_reannotations.tsv'
    else:
        fetal_data_out_path = './data/processed/fetal_brainarea_vs_genes_exp.tsv'

    # fetal_data_out_path = './data/processed/fetal_brainarea_vs_genes_exp.tsv'

    if os.path.exists(fetal_data_out_path):
        print('Processed fetal brain dataset found locally. Loading from {}'.format(fetal_data_out_path))
        structure_genes_exp_matrix = pd.read_csv(fetal_data_out_path,
                                                 index_col='gene_symbol',
                                                 sep='\t')

    else:
        raw_fetal_data_path = './data/raw/allen_human_fetal_brain'
        fetal_donor_folders = glob.glob(os.path.join(raw_fetal_data_path, '*'))

        all_donors = []

        for donor in fetal_donor_folders:
            donor_id = donor.split('/')[-1].split('_')[-1]
            donor_files = glob.glob(os.path.join(donor, '*'))
            expression, samples, probes = get_donor_data(donor_files,
                                                         reannotate=reannotate)
            tidy_donor = get_single_donor_tidy_df(expression,
                                                  samples,
                                                  probes,
                                                  donor_id=donor_id,
                                                  reannotate=False)
            all_donors.append(tidy_donor)

        all_donors_long = pd.concat(all_donors)
        structure_genes_exp_matrix = pd.pivot_table(all_donors_long,
                                                    values='value',
                                                    index='gene_symbol',
                                                    columns='brain_area')
        print('-- Writing data to ' + fetal_data_out_path + ' -- ')
        structure_genes_exp_matrix.to_csv(fetal_data_out_path, sep='\t')
    return structure_genes_exp_matrix
