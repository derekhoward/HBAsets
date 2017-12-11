import matplotlib.pyplot as plt
import seaborn as sns
import plotly
import plotly.figure_factory as ff


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


def make_dist_plot(exp_matrix, brain_area, gene_list, ax, gene_text=False):
    """
    Plots the distributions of zscored levels of expression of genes of
    interest compared to the rest of genes

    Parameters
    ----------
    exp_matrix : dataframe
        expression matrix: rows->genes ; columns->brain_areas
    brain_area : str
        select brain area to be plotted
    gene_list : series
        list of gene symbols of interest
    Returns
    -------

    """
    brain_area_exp = exp_matrix.loc[:, brain_area].reset_index()
    
    # get the data for genes of interest
    brain_area_exp['In gene list'] = brain_area_exp.gene_symbol.isin(gene_list)
    plot_info = brain_area_exp[brain_area_exp['In gene list']].loc[:, ['gene_symbol', brain_area]]
    disease_vals = plot_info.loc[:, brain_area]
    disease_names = plot_info.loc[:, 'gene_symbol']
    names = []
    exp_vals = []
    for ix, info in plot_info.iterrows():
        names.append(info[0])
        exp_vals.append(info[1])

    background_vals = brain_area_exp[~brain_area_exp['In gene list']].loc[:, brain_area]

    kde = gaussian_kde(background_vals)

    # these are the values over which your kernel will be evaluated
    dist_space = np.linspace(min(background_vals), max(background_vals), 100)
    
    # plot the results
    with sns.plotting_context('notebook', font_scale = 1.25):
        distplot = ax.plot(dist_space, kde(dist_space), color='k', lw=1.5)
        ax.set_ylabel('Density')

        # plot lines over the distribution for the disease genes
        palette = itertools.cycle(sns.color_palette())
        i = 0

        for name, position, colour in zip(names, exp_vals, palette):
            ax.vlines(position, ymin=-0.1, ymax=1, lw=1.5, colors=colour)
            if gene_text is True:
                plt.text(6, 5.25-(i/10), name, fontsize=16, color=colour)
            i+=1
    sns.despine()
    return distplot


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
        if value:
            plt.vlines(position, ymin=0, ymax=1, **kwargs)
            i += 1
    print('Number of genes of interest: {}'.format(i))
    plt.ylim(0, 1)
    plt.xlim(0, len(gene_list))
    ax.set_yticklabels([])
    return ax


def make_violins(exp, brain_areas, gene_list):
    subset = exp.loc[:, brain_areas]
    subset['in_gene_list'] = subset.index.isin(gene_list)
    tidy = subset.reset_index().melt(id_vars=['gene_symbol', 'in_gene_list'], var_name='brain_area', value_name='expression')
    
    with sns.plotting_context('notebook', font_scale=1.25):
        fig, ax = plt.subplots(figsize=(12, 9))
        sns.violinplot(y='brain_area', x='expression', edgecolor='black', hue='in_gene_list', palette={False:'#636364', True:'#D9D4D3'}, cut=2, split=True, inner='quartiles', data=tidy, ax=ax)
        ax.set_xlabel('Expression (z-scored)')
        ax.set_ylabel('')

        legend = ax.get_legend()
        legend.set_title('')
        legend._loc = 7
        legend_labels = {'False': 'Background genes', 'True': 'Disease genes'}
        for text, label in zip(legend.texts, legend_labels.items()): 
            text.set_text(label[1])

        sns.despine()


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
