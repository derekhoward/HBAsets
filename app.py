import dash
from dash.dependencies import Input, Output
import dash_core_components as dcc
import dash_html_components as html
import pandas as pd
from HBA_analysis import generate_stats_table

exp = pd.read_csv('./data/processed/brainarea_vs_genes_exp_w_reannotations.tsv',
                  index_col='gene_symbol', sep='\t')


def generate_table(dataframe, max_rows=10):
    return html.Table(
        # Header
        [html.Tr([html.Th(col) for col in dataframe.columns])] +

        # Body
        [html.Tr([
            html.Td(dataframe.iloc[i][col]) for col in dataframe.columns
        ]) for i in range(min(len(dataframe), max_rows))]
    )


app = dash.Dash()

app.layout = html.Div(children=[
    dcc.Textarea(id='input_genes',
                 placeholder='Input list of genes here',
                 value='PDGFRB PDGFB SLC20A2'),
    html.H2(children='HBASets'),
    # generate_table()
    html.Div(id='my_table')
])


@app.callback(
    Output(component_id='my_table', component_property='children'),
    [Input(component_id='input_genes', component_property='value')]
)
def update_table(input_value):
    gene_list = pd.Series(input_value.split())
    output_dataframe = generate_stats_table(exp, gene_list).reset_index()
    return generate_table(output_dataframe)


if __name__ == '__main__':
    app.run_server(debug=True)
