import pandas as pd
from scipy.stats import mannwhitneyu
import statsmodels.sandbox.stats.multicomp
from sklearn import metrics
import matplotlib.pyplot as plt
import seaborn as sns
import plotly
import plotly.figure_factory as ff


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


def make_ROC_plot(exp_df, brain_area, gene_list):
    """
    Creates static ROC plot for genes of interest in specified brain area

    Parameters
    ----------
    exp_df: dataframe
        expression matrix: rows->genes ; columns->brain_areas
    brain_area: string
        name of brain area testing within
    gene_list: series
        list of gene symbols of interest
    Returns
    -------
    ax
        axis containing ROC plot
    """
    y_true = exp_df[brain_area].index.isin(gene_list)
    y_score = exp_df[brain_area].values
    # sns.distplot(y_score)
    fpr, tpr, thresholds = metrics.roc_curve(y_true, y_score, pos_label=True)
    brainarea_auc = metrics.auc(fpr, tpr)

    plt.figure(figsize=(8, 6), dpi=80)
    ax = plt.gca()
    plt.plot(fpr, tpr, color='darkorange',
             label='ROC (area = {:0.6f})'.format(brainarea_auc))
    plt.plot([0, 1], [0, 1], color='navy', linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('{} ROC'.format(brain_area))
    plt.legend(loc="lower right")

    return ax


def make_dist_plot(exp_matrix, brain_area, gene_list):
    """
    Plots the distributions of zscored levels of expression of genes of
    interest compared to the rest of genes

    Parameters
    ----------
    exp_df : dataframe
        expression matrix: rows->genes ; columns->brain_areas
    gene_list : series
        list of gene symbols of interest
    Returns
    -------

    """
    brain_area_exp = exp_matrix.loc[:, brain_area].reset_index()
    brain_area_exp['In gene list'] = brain_area_exp.gene_symbol.isin(gene_list)
    g = sns.FacetGrid(brain_area_exp, row='In gene list', aspect=4)
    g.map(sns.distplot, brain_area, norm_hist=True, hist=True, rug=False)
    g.set_xlabels('Mean expression value')


def raster(gene_list, **kwargs):
    """
    Creates a raster plot

    Parameters
    ----------
    gene_list: iterable
        a list of genes of interest
    color: string
        color of vlines
    Returns
    -------
    ax
        an axis containing the raster plot
    """
    plt.figure(figsize=(8, 1), dpi=80)
    ax = plt.gca()

    i = 0
    for position, value in enumerate(gene_list[::-1]):
        # somehow this doesnt plot if == is switched to "is" ...?
        if value:
            plt.vlines(position, ymin=0, ymax=1, **kwargs)
            i += 1
    print('Number of genes of interest: {}'.format(i))
    plt.ylim(0, 1)
    plt.xlim(0, len(gene_list))
    ax.set_yticklabels([])
    return ax


def interactive_distplot(exp_df, brain_area, gene_list, disease_label=None, filename='Distplot with Normal Curve'):
    brain_area_exp = exp_df.loc[:, brain_area].reset_index()
    brain_area_exp['hit'] = brain_area_exp.gene_symbol.isin(gene_list)

    x1 = brain_area_exp[brain_area_exp['hit'] == True].loc[:, brain_area].values
    x2 = brain_area_exp[brain_area_exp['hit'] == False].loc[:, brain_area].values

    hist_data = [x1, x2]
    group_labels = ['{} hits'.format(disease_label), 'Background']
    # colors = ['#94F3E4', '#3A4750']

    rug_text_hits = brain_area_exp[brain_area_exp['hit'] == True].loc[:, 'gene_symbol']
    rug_text_background = brain_area_exp[brain_area_exp['hit'] == False].loc[:, 'gene_symbol']
    rug_text = [rug_text_hits, rug_text_background]

    fig = ff.create_distplot(hist_data, group_labels, bin_size=.25, curve_type='normal', rug_text=rug_text)#, colors=colors)
    plotly.offline.iplot(fig, filename=filename)


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
