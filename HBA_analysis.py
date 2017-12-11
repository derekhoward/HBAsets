import pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.sandbox.stats.multicomp
from sklearn import metrics


def calc_raw_pvalues(exp_df, gene_list):
    """
    Gets uncorrected p_values for whether a geneset is enriched in a structure.

    Parameters
    ----------
    exp_df: dataframe
        expression matrix: rows->genes ; columns->brain_areas
    gene_list: series
        list of gene symbols of interest
    Returns
    -------
    results
        a series of P-values for each brain area
    """
    results = []
    for col in exp_df:
        X = exp_df[col][exp_df[col].index.isin(gene_list)]
        y = exp_df[col][~exp_df[col].index.isin(gene_list)]
        # make sure that you are comparing the gene list to all other genes
        # assert(X.shape[0] + y.shape[0] == exp_df.shape[0])
        results.append(mannwhitneyu(X, y, alternative='two-sided')[1])
    return pd.Series(results)


def calc_AUC_for_ranked_genelist(exp_df, gene_list):
    """
    Calculates AUROC to determine whether genes of interest are expressed at
    levels greater than chance within the brain area

    Parameters
    ----------
    exp_df: dataframe
        expression matrix: rows->genes ; columns->brain_areas
    gene_list: series
        list of gene symbols of interest
    Returns
    -------
    results
        a series of AUC values for each brain area
    """
    results = []
    for col in exp_df:
        y_true = exp_df[col].index.isin(gene_list)
        y_score = exp_df[col].values
        results.append(metrics.roc_auc_score(y_true, y_score))
    return pd.Series(results)


def generate_stats_table(exp_df, gene_list):
    """
    Creates a table of summary stats for each brain area

    Parameters
    ----------
    exp_df : dataframe
        expression matrix: rows->genes ; columns->brain_areas
    gene_list : series
        list of gene symbols of interest
    Returns
    -------
    table
        dataframe with brain areas as index

    """
    count = len(gene_list)
    n_genes_in_matrix = gene_list.isin(exp_df.index).sum()
    genes_not_found = gene_list[~gene_list.isin(exp_df.index)].values

    print('You submitted a gene list with {} genes.\n\
{} of those genes are present in the HBA dataset.\n\
Genes not found in our reference data: {}'.format(
            count, n_genes_in_matrix, genes_not_found))

    raw_pvalues = calc_raw_pvalues(exp_df, gene_list)
    corrected_p_vals = statsmodels.sandbox.stats.multicomp.multipletests(
            raw_pvalues, method="fdr_bh")[1]
    corrected_p_vals = pd.Series(corrected_p_vals)
    auc = calc_AUC_for_ranked_genelist(exp_df, gene_list)
    table = pd.concat([raw_pvalues, corrected_p_vals, auc],
                      keys=['raw p_values', 'corrected_p', 'AUC'], axis=1)
    table.set_index(exp_df.columns, inplace=True)
    return table.sort_values('AUC', ascending=False)
