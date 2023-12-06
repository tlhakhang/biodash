from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash.dependencies import Output, Input, State
from functools import reduce
import scipy
import scipy.stats
from functools import reduce
import numpy 
import os
import pandas as pd
import dash_core_components as dcc
import dash_html_components as html
from sklearn.preprocessing import scale
import dash_table
from flask_login import current_user
from app.models import User, Post, Ownership
import base64
import plotly.express as px
import dash
import scipy.cluster.hierarchy as shc
import numpy as np
from sklearn.manifold import TSNE
import plotly.express as px
from umap import UMAP
import dash_dangerously_set_inner_html
from markdownify import markdownify
import scanpy as sc
import matplotlib.pyplot as plt
from io import BytesIO
from base64 import b64encode
#from gseapy.plot import barplot, dotplot
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import subprocess
import sys
import json


def fig_to_uri(in_fig, close_all=True, **save_args):
    # type: (plt.Figure) -> str
    """
    Save a figure as a URI
    :param in_fig:
    :return:
    """
    out_img = BytesIO()
    in_fig.savefig(out_img, format='png', **save_args)
    if close_all:
        in_fig.clf()
        plt.close('all')
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)





def get_mean_expression(input_matrix,chamber_classification):
    subset_matrix = input_matrix.loc[:, '# Cells':'Cell_Identifier']
    subset_matrix = subset_matrix.loc[subset_matrix['Cell_Identifier'] == chamber_classification]
    del subset_matrix["# Cells"]
    del subset_matrix["Cell_Identifier"]
    
    # final = pd.DataFrame(np.nan_to_num(np.log(subset_matrix),neginf=0))
    # final.index = subset_matrix.index
    # final.columns = subset_matrix.columns
    
    
    mean = pd.DataFrame(subset_matrix.mean(axis=0))
    mean['Chamber Classification'] = chamber_classification
    mean.columns = ['Mean Expression','classification']
    mean['Protein'] = mean.index
    return(mean)

        
def register_callbacks(dashapp):
    
     
    @dashapp.callback(
        dash.dependencies.Output('intermediate-value-illumina-number', 'children'),
        [dash.dependencies.Input('button', 'n_clicks')],
        [dash.dependencies.State('my-dropdown', 'value')])
    def return_transcriptome_cellbarcodes(n_clicks,value):
        
        path_to_run = value.split('/')
        path_to_run.pop()
        last_elem = path_to_run.pop()
        
        sample_id = last_elem.split('_')[0]
        
        path_to_run = '/'.join(path_to_run)
        run_path = path_to_run + '/'+last_elem
        
        
        matrix = run_path + '/Results/force5000_output/alevin/quants_mat.mtx.gz'
        
        adata = sc.read_mtx(matrix)
        
        number_of_cells = len(adata)
        string_value = 'Number of Cell Barcodes:'+ str(number_of_cells)

        return(string_value)
    

       # return([])

    
    
        
        
    
    @dashapp.callback(
        dash.dependencies.Output('intermediate-value-scanpy', 'data'),
        [dash.dependencies.Input('button', 'n_clicks')],
        [dash.dependencies.State('my-dropdown', 'value')])
    def get_alevin_scanpy_results(n_clicks,value):
        
        path_to_run = value.split('/')
        path_to_run.pop()
        last_elem = path_to_run.pop()
        
        sample_id = last_elem.split('_')[0]
        
        path_to_run = '/'.join(path_to_run)
        run_path = path_to_run + '/'+last_elem

        
        matrix = run_path + '/Results/force5000_output/alevin/quants_mat.mtx.gz'
        cell_bar = run_path + '/Results/force5000_output/alevin/quants_mat_rows.txt'
        genes = run_path + '/Results/force5000_output/alevin/quants_mat_cols.txt'
        

        adata = sc.read_mtx(matrix)
        cell_barcodes = pd.read_csv(cell_bar,header=None)
        cell_barcodes.columns = ['cell_barcodes']

        gene_names = pd.read_csv(genes,header=None)
        gene_names.columns = ['genesymbol']

        #cell_barcodes = cell_barcodes['cell_barcodes'].to_list()
        adata.obs = cell_barcodes
        adata.var_names = gene_names['genesymbol'].to_list()

        print(adata.obs)
        genes = pd.read_csv('/home/isodev/development/dev/scripts/gene.mapping.csv')
        #print(genes)
        #genes = pd.read_csv('/home/duomic/demo/dev/scripts/gene.mapping.csv')
        adata.var_names = genes['genesymbol']
        adata.var.index = adata.var.index.map(str)
        adata.var_names_make_unique()




        adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

        
        gene_metrics = adata.var.sort_values(by=['total_counts'],ascending=False)[1:100]
        gene_metrics['gene'] = gene_metrics.index
        gene_metrics_subset = gene_metrics.sort_values(by=['total_counts'], ascending=False)[1:50]
        
#         enr_res = gseapy.enrichr(gene_list=gene_metrics.index.to_list(),
#                      gene_sets='GO_Biological_Process_2018',
#                      description='pathway',
#                      cutoff = 0.1)
        
        
        
        fig = make_subplots(rows=1, cols=6)

        fig.add_trace(go.Box(y=adata.obs['n_genes_by_counts']),row=1, col=1)
        fig.add_trace(go.Box(y=adata.obs['total_counts']),row=1, col=2)
        fig.add_trace(go.Box(y=adata.obs['pct_counts_mt']),row=1, col=3)
        fig.add_trace(go.Bar(x=gene_metrics_subset['total_counts'], y=gene_metrics_subset['gene'],orientation='h'),row=1,col=4)
        
        fig.update_layout(yaxis={'categoryorder':'total ascending'}) # add only this line
        #fig.add_trace(go.Bar(x=gene_metrics_subset['total_counts'], y=gene_metrics_subset['gene'],orientation='h'),)
        
        observations  = run_path + '/Results/force5000_output/alevin/featureDump.txt'
        obs = pd.read_csv(observations,sep='\t')
        
        fig.add_trace(go.Violin(y=obs['DedupRate'],
                            box_visible=True,
                            name='Deduplication Rate',
                            meanline_visible=True,points='all'),row=1,col=5)
        
        fig.add_trace(go.Violin(y=obs['MappingRate'],
                            box_visible=True,
                            name='Mapping Rate',
                            meanline_visible=True,points='all'),row=1,col=6)
        #gene_metrics_subset = gene_metrics.sort_values(by=['total_counts'], ascending=False)[1:50]
        
        


        fig.update_xaxes(title_text="Number of genes per cell",row=1, col=1)
        fig.update_xaxes(title_text="Total Counts", row=1, col=2)
        fig.update_xaxes(title_text="Percent Mitochondrial", row=1, col=3)
        fig.update_xaxes(title_text="Top 50 Genes by Total UMI Counts", row=1, col=4)


        fig.update_layout(height=900, width=1800, title_text="Summary QC Metrics")
        #fig.update_traces(boxpoints=False) 
        fig.update_yaxes(tickfont_size=18)
        fig.update_layout(
            font=dict(
                family="Arial",
                size=14,
            ))
        
        
        
        

        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)

        sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    
        adata.raw = adata
        adata = adata[:, adata.var.highly_variable]




        #sc.pp.regress_out(adata, ['total_counts'])
        sc.pp.scale(adata, max_value=10)
        sc.tl.pca(adata, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40)
        sc.tl.umap(adata,n_components =2)
        sc.tl.leiden(adata)

        sc.tl.rank_genes_groups(adata, 'leiden', method='t-test')
        result = adata.uns['rank_genes_groups']
        groups = result['names'].dtype.names
            
        ranked_df = pd.DataFrame(
                {'Cluster ' + group + '_' + key[:1]: result[key][group]
                for group in groups for key in ['names', 'pvals']})[1:100]

        print(ranked_df.head())
        adata_obs = adata.obs


        #print(out_url)

        df = pd.DataFrame((adata.obsm['X_umap']))
        df['n_genes'] = adata.obs['n_genes_by_counts'].values
        df['cell'] = adata.obs['cell_barcodes'].to_list()
        df['leiden'] = adata.obs['leiden'].values
        df.columns = ['X','Y','n_genes_by_counts','cell','leiden']
        
        
        print(df.head())

        
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        
        return(df.to_json(date_format='iso', orient='split'),ranked_df.to_json(date_format='iso', orient='split'))

    @dashapp.callback(
    dash.dependencies.Output('scanpy-umap', 'figure'),
    [dash.dependencies.Input('intermediate-value-scanpy', 'data')])

    def get_umap(json_file):
        

        df = pd.read_json(json_file[0],orient='split')
        
        fig = px.scatter(df, x='X', y='Y',color='leiden',hover_data=['cell'])
        fig.update_layout(height=700, width=1800)
        fig.update_layout(
            title={
            'text' : 'UMAP - leiden',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_layout(
            font=dict(

                family="Arial",
                size=11,
            ))
        
        

        return(fig)

    @dashapp.callback(Output('marker_rank_table', 'children'),
        Input('intermediate-value-scanpy', 'data'))
    def update_output(json_file):
        #
        df = pd.read_json(json_file[1],orient='split')


               # print(df.head())
        table = html.Div([dash_table.DataTable(
            columns=[{"name": i, "id": i} for i in df.columns],
            style_table={
                'overflowX': 'scroll',
                'width':'100%',
                'margin':'auto'},
            data=df.to_dict('records'),
            virtualization=True,
            page_size=40,
            sort_action="native",
            style_header={'backgroundColor': 'rgb(30, 30, 30)'},

            style_filter={
                'backgroundColor': '#A5A3A6',
                'color': 'black'
            },  
            style_cell={
                'backgroundColor': 'rgb(0, 25, 51)',
                'color': 'white',
                'textAlign':'center',
                'fontSize':15, 
                'font-family':'sans-serif'},
            
            export_format='xlsx',
            export_headers='display',
            merge_duplicate_headers=True)]
        )
        
                    

        return(table)
                                    
                                    
                                    
    
    @dashapp.callback(
        dash.dependencies.Output('scanpy-qc', 'figure'),
        [dash.dependencies.Input('button', 'n_clicks')],
        [dash.dependencies.State('my-dropdown', 'value')])

    def return_transcriptome_qc(n_clicks,value):
        

        path_to_run = value.split('/')
        path_to_run.pop()
        last_elem = path_to_run.pop()
        
        sample_id = last_elem.split('_')[0]
        
        path_to_run = '/'.join(path_to_run)
        run_path = path_to_run + '/'+last_elem

        
        matrix = run_path + '/Results/force5000_output/alevin/quants_mat.mtx.gz'
        cell_bar = run_path + '/Results/force5000_output/alevin/quants_mat_rows.txt'
        genes = run_path + '/Results/force5000_output/alevin/quants_mat_cols.txt'
        

        adata = sc.read_mtx(matrix)
        cell_barcodes = pd.read_csv(cell_bar,header=None)
        cell_barcodes.columns = ['cell_barcodes']

        gene_names = pd.read_csv(genes,header=None)
        gene_names.columns = ['genesymbol']

        #cell_barcodes = cell_barcodes['cell_barcodes'].to_list()
        adata.obs = cell_barcodes
        adata.var_names = gene_names['genesymbol'].to_list()

        print(adata.obs)
        genes = pd.read_csv('/home/isodev/development/dev/scripts/gene.mapping.csv')
        #print(genes)
        #genes = pd.read_csv('/home/duomic/demo/dev/scripts/gene.mapping.csv')
        adata.var_names = genes['genesymbol']
        adata.var.index = adata.var.index.map(str)
        adata.var_names_make_unique()




        adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


        
    
        gene_metrics = adata.var.sort_values(by=['total_counts'],ascending=False)[1:100]
        gene_metrics['gene'] = gene_metrics.index
        gene_metrics_subset = gene_metrics.sort_values(by=['total_counts'], ascending=False)[1:50]
        
#         enr_res = gseapy.enrichr(gene_list=gene_metrics.index.to_list(),
#                      gene_sets='GO_Biological_Process_2018',
#                      description='pathway',
#                      cutoff = 0.1)
        
        
        
        fig = make_subplots(rows=1, cols=6)

        fig.add_trace(go.Box(y=adata.obs['n_genes_by_counts']),row=1, col=1)
        fig.add_trace(go.Box(y=adata.obs['total_counts']),row=1, col=2)
        fig.add_trace(go.Box(y=adata.obs['pct_counts_mt']),row=1, col=3)
        fig.add_trace(go.Bar(x=gene_metrics_subset['total_counts'], y=gene_metrics_subset['gene'],orientation='h'),row=1,col=4)
        
        fig.update_layout(yaxis={'categoryorder':'total ascending'}) # add only this line
        #fig.add_trace(go.Bar(x=gene_metrics_subset['total_counts'], y=gene_metrics_subset['gene'],orientation='h'),)
        
        observations  = run_path + '/Results/force5000_output/alevin/featureDump.txt'
        obs = pd.read_csv(observations,sep='\t')
        
        fig.add_trace(go.Violin(y=obs['DedupRate'],
                            box_visible=True,
                            name='Deduplication Rate',
                            meanline_visible=True,points='all'),row=1,col=5)
        
        fig.add_trace(go.Violin(y=obs['MappingRate'],
                            box_visible=True,
                            name='Mapping Rate',
                            meanline_visible=True,points='all'),row=1,col=6)
        #gene_metrics_subset = gene_metrics.sort_values(by=['total_counts'], ascending=False)[1:50]
        
        


        fig.update_xaxes(title_text="Number of genes per cell",row=1, col=1)
        fig.update_xaxes(title_text="Total Counts", row=1, col=2)
        fig.update_xaxes(title_text="Percent Mitochondrial", row=1, col=3)
        fig.update_xaxes(title_text="Top 50 Genes by Total UMI Counts", row=1, col=4)


        fig.update_layout(height=900, width=1800, title_text="Summary QC Metrics")
        #fig.update_traces(boxpoints=False) 
        fig.update_yaxes(tickfont_size=18)
        fig.update_layout(
            font=dict(
                family="Arial",
                size=11,
            ))

        return (fig)
    
    
    @dashapp.callback(
    dash.dependencies.Output('scanpy-umap-topgenes', 'figure'),
    [dash.dependencies.Input('scanpy-umap', 'selectedData')],
    [State('my-dropdown','value')])
    def select_top_gene_barplot_umap(selectedData,value):
        
        
        print(selectedData)
        
        cell_barcode_to_filter_list = [x['customdata'] for x in selectedData['points']]
        cell_barcode_to_filter = [item for sublist in cell_barcode_to_filter_list for item in sublist]
        
        
        path_to_run = value.split('/')
        path_to_run.pop()
        last_elem = path_to_run.pop()
        
        sample_id = last_elem.split('_')[0]
        
        path_to_run = '/'.join(path_to_run)
        run_path = path_to_run + '/'+last_elem

        
        matrix = run_path + '/Results/force5000_output/alevin/quants_mat.mtx.gz'
        cell_bar = run_path + '/Results/force5000_output/alevin/quants_mat_rows.txt'
        genes = run_path + '/Results/force5000_output/alevin/quants_mat_cols.txt'
        


        adata = sc.read_mtx(matrix)
        cell_barcodes = pd.read_csv(cell_bar,header=None)
        cell_barcodes.columns = ['cell_barcodes']

        gene_names = pd.read_csv(genes,header=None)
        gene_names.columns = ['genesymbol']

        #cell_barcodes = cell_barcodes['cell_barcodes'].to_list()
        adata.obs.index = cell_barcodes['cell_barcodes'].to_list()
        adata.var_names = gene_names['genesymbol'].to_list()

        
        genes = pd.read_csv('/home/isodev/development/dev/scripts/gene.mapping.csv')
        #print(genes)
        #genes = pd.read_csv('/home/duomic/demo/dev/scripts/gene.mapping.csv')
        adata.var_names = genes['genesymbol']
        adata.var.index = adata.var.index.map(str)
        adata.var_names_make_unique()

        adata =adata[cell_barcode_to_filter].copy()


        adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
        sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)




        gene_metrics = adata.var.sort_values(by=['total_counts'],ascending=False)[1:100]
        gene_metrics['gene'] = gene_metrics.index
        gene_metrics_subset = gene_metrics.sort_values(by=['total_counts'], ascending=False)[1:50]

        fig = make_subplots(rows=1, cols=1)

        fig.add_trace(go.Bar(x=gene_metrics_subset['total_counts'], y=gene_metrics_subset['gene'],orientation='h'),row=1,col=1)

        fig.update_layout(yaxis={'categoryorder':'total ascending'}) # add only this line

        fig.update_layout(
            title={
            'text' : 'Top Expressed Genes - UMAP Selection',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig.update_layout(
            font=dict(
                family="Arial",
                size=11,
            ))
        
        
        fig.update_layout(height=900, width=900, title_text="Top Expressed Genes - UMAP Selection")
        fig.update_layout(xaxis_title="Total UMI Counts", yaxis_title="Gene Symbol")

        return(fig)
    
    
    
    
    @dashapp.callback(
        dash.dependencies.Output('intermediate-value', 'children'),
        [dash.dependencies.Input('button', 'n_clicks')],
        [dash.dependencies.State('my-dropdown', 'value')])
    def return_transcriptome_alevin_report_5000(n_clicks,value):
        

        for root, dirs, files in os.walk(value):
            for file in files:
                if file.endswith('5000_output.html'):
                    p = os.path.join(root,file)
                    html_file = open(p,"r").read()
                    html_file_obj = html.Iframe(srcDoc=html_file,style={"height": 800, "width": 1100}) 
        
        return(html_file_obj)

    
    @dashapp.callback(
        dash.dependencies.Output('intermediate-value-2', 'children'),
        [dash.dependencies.Input('button', 'n_clicks')],
        [dash.dependencies.State('my-dropdown', 'value')])
    def return_transcriptome_alevin_report_default(n_clicks,value):
        

        for root, dirs, files in os.walk(value):
            for file in files:
                if file.endswith('default_output.html'):
                    p = os.path.join(root,file)
                    html_file = open(p,"r").read()
                    html_file_obj = html.Iframe(srcDoc=html_file,style={"height": 800, "width": 1100}) 
        
        return(html_file_obj)

        
        
    
    @dashapp.callback(
    dash.dependencies.Output('intermediate-value-fastqc-r1', 'children'),
    [dash.dependencies.Input('button', 'n_clicks')],
    [dash.dependencies.State('my-dropdown', 'value')])
    
    def return_transcriptome_r1_fastqc_report(n_clicks,value):
        
        print(value)
        import dash_dangerously_set_inner_html
        from markdownify import markdownify
        
        for root, dirs, files in os.walk(value+'/Raw/'):
            for file in files:
                if file.endswith('R1_fastqc.html'):
                    p = os.path.join(root,file)
                    html_file = open(p,"r").read()

        
        html_file_obj = html.Iframe(srcDoc=html_file,style={"height": 800, "width": 1100}) 
        #html_file_obj = dash_dangerously_set_inner_html.DangerouslySetInnerHTML(html_file)
        return(html_file_obj)
    
    
    @dashapp.callback(
    dash.dependencies.Output('intermediate-value-fastqc-r2', 'children'),
    [dash.dependencies.Input('button', 'n_clicks')],
    [dash.dependencies.State('my-dropdown', 'value')])
    
    def return_transcriptome_r2_fastqc_report(n_clicks,value):
        
        print(value)
        import dash_dangerously_set_inner_html
        from markdownify import markdownify
        
        for root, dirs, files in os.walk(value+'/Raw/'):
            for file in files:
                if file.endswith('R2_fastqc.html'):
                    p = os.path.join(root,file)
                    html_file = open(p,"r").read()

        
        html_file_obj = html.Iframe(srcDoc=html_file,style={"height": 800, "width": 1100}) 
        #html_file_obj = dash_dangerously_set_inner_html.DangerouslySetInnerHTML(html_file)
        return(html_file_obj)
    
    
    
    @dashapp.callback(dash.dependencies.Output('fastqscreen-html', 'children'),
    [dash.dependencies.Input('button', 'n_clicks')],
    [dash.dependencies.State('my-dropdown', 'value')])
    
    def return_transcriptome_r2_fastqscreen_report(n_clicks,value):
        
        print(value)
        import dash_dangerously_set_inner_html
        from markdownify import markdownify
        
        for root, dirs, files in os.walk(value+'/Trimmed/'):
            for file in files:
                if file.endswith('R2_screen.html'):
                    p = os.path.join(root,file)
                    html_file = open(p,"r").read()

        
        html_file_obj = html.Iframe(srcDoc=html_file,style={"height": 800, "width": 1400}) 
        #html_file_obj = dash_dangerously_set_inner_html.DangerouslySetInnerHTML(html_file)
        return(html_file_obj)
    
    

    
    @dashapp.callback(
    dash.dependencies.Output('protein_table_hidden', 'children'),
    [dash.dependencies.Input('button_protein', 'n_clicks')],
    [dash.dependencies.State('protein-dropdown', 'value')])
    def get_protein_concat_table(n_clicks,value):
        
        
        if n_clicks:
            input_files = []
            

            print(value)
            for root, dirs, files in os.walk(value):
                for file in files:
                    if file.endswith('Table.csv'):
                        cell_classifier = root.split('/')[6]
                        p = os.path.join(root,file)
                        table = pd.read_csv(p,index_col=0)
                        table['Cell_Identifier'] = cell_classifier
                        input_files.append(table)
                        
            print(input_files,"Input Files########")
            df = (pd.concat(input_files,ignore_index=True))
            
            #import os
            directory = value+'/Total/Tables/'
            if not os.path.exists(value+'/Total/Tables/'):
                os.makedirs(directory)
                df.to_csv(directory+'Protein.Total.csv',index=False)
            #print(df)
            return(df.to_json())

    @dashapp.callback(
    dash.dependencies.Output('protein_combined_table', 'children'),
    [dash.dependencies.Input('protein_table_hidden', 'children')],
    [State('protein_table_hidden','children')])
    def return_protein_concat_table(n_clicks, value):


        df = pd.read_json(value)

       # print(df.head())
        table = html.Div([dash_table.DataTable(
                columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records'),

                page_size=15,
                sort_action="native",
                style_header={'backgroundColor': 'rgb(30, 30, 30)'},

                style_filter={
                    'backgroundColor': '#A5A3A6',
                    'color': 'black'
                },  
                style_cell={
                    'backgroundColor': 'rgb(0, 25, 51)',
                    'color': 'white',
                    'textAlign':'center',
                    'fontSize':12, 
                    'font-family':'sans-serif'},
                export_format='xlsx',
                export_headers='display',
                merge_duplicate_headers=True)]
            )

        return([table])

    
    @dashapp.callback(
    dash.dependencies.Output('protein_barplot', 'figure'),
    [dash.dependencies.Input('protein_table_hidden', 'children')],
    [State('protein_table_hidden','children')])
    def return_protein_cell_barplot(input_1,value):
        


        df = pd.read_json(value)
        cell_table = pd.DataFrame((df['Cell_Identifier'].value_counts()))
        fig = px.bar(cell_table,x=cell_table.index,y='Cell_Identifier')
        fig.update_layout(
            title={
            'text' : 'Bar Plot - Chambers Identified vs. Classification',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig.update_layout(
            font=dict(
                family="Arial",
                size=11,
            ))
        

        return(fig)
    
    @dashapp.callback(
    dash.dependencies.Output('protein_scatterplot', 'figure'),
    [dash.dependencies.Input('protein_table_hidden', 'children'),
    dash.dependencies.Input('protein_umap', 'selectedData')],
    [State('protein_table_hidden','children')])
    def return_protein_location_scatter(input_1,selectedData,value):
        
        

        df = pd.read_json(value)
        df['Chamber_Column'] =df['Chamber Column'].astype(str)
        df['Chamber_Row'] = df['Chamber Row'].astype(str)
        df['Location'] = df['Chamber_Column'].str.cat(df[['Chamber_Row']], sep='|')
        
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        
        
        #points = (json.dumps(selectedData, indent=2))
        
        print(selectedData)
        
        location = [x['customdata'] for x in selectedData['points']]
        newlist = [item for items in location for item in items]
        new_df = df[df['Location'].isin(newlist)]
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        

        fig = px.scatter(new_df,x='Chamber Column',y='Chamber Row',color="Cell_Identifier")
        fig.update_layout(
            title={
            'text' : 'Protein Chambers - Chip Location',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig.update_yaxes(autorange="reversed")
        fig.update_layout(
            font=dict(
                family="Arial",
                size=11,
            ))
        

        return(fig)
    
    

    
    
    
    @dashapp.callback(
    dash.dependencies.Output('protein_barplot_mean', 'figure'),
    [dash.dependencies.Input('protein_table_hidden', 'children')],
    [State('protein_table_hidden','children')])
    def return_protein_expression_barplot(input_1,value):
        
        df = pd.read_json(value)
        


        zero = get_mean_expression(df,'Zero Cell')
        single  = get_mean_expression(df,'Single Cell')
        two  = get_mean_expression(df,'Two Cell')
        three  = get_mean_expression(df,'Multi Cell')
        
        combined = (pd.concat([zero,single,two,three]))
        
        fig = px.bar(combined, x="Mean Expression", y="Protein",
             color='classification', barmode='group',orientation='h',
             height=800,width=1600)
        
        fig.update_layout(
            title={
            'text' : 'Mean Expression vs. Chamber Classification',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig.update_yaxes(autorange="reversed")
        fig.update_layout(
            font=dict(
                family="Arial",
                size=17,
            ))
        

        return(fig)
        
        
        
    @dashapp.callback(
    dash.dependencies.Output('protein_boxplot', 'figure'),
    [dash.dependencies.Input('protein_table_hidden', 'children')],
    [State('protein_table_hidden','children')])
    def return_protein_chip_sum_scatterplot(input_1,value):
        
        df = pd.read_json(value)
        matrix = df.loc[:, '# Cells':'Cell_Identifier']
        del matrix['# Cells']
        del matrix['Cell_Identifier']
        
        cell_sum = pd.DataFrame(matrix.sum(axis=1))
        cell_sum.columns = ['Sum Signal']
        cell_sum['Sum Signal'] = cell_sum['Sum Signal']*10
        cell_sum['Chamber Column'] = df['Chamber Column']
        cell_sum['Chamber Row'] = df['Chamber Row']
        cell_sum['Cell_Identifier'] = df['Cell_Identifier']
        fig = px.scatter(cell_sum, x="Chamber Column", y="Chamber Row", color="Cell_Identifier",hover_data=['Sum Signal'],size=cell_sum["Sum Signal"].clip(0, 600000))



        fig.update_layout(
            title={
            'text' : 'Sum Chamber Intensities',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig['layout']['yaxis']['autorange'] = "reversed"
        #fig.update_yaxes(autorange="reversed")
        fig.update_layout(
            font=dict(
                family="Arial",
                size=11,
            ))
        

        return(fig)
    
    
    
    

    @dashapp.callback(
    dash.dependencies.Output('protein_umap', 'figure'),
    [dash.dependencies.Input('protein_table_hidden', 'children')],
    [State('protein_table_hidden','children')])
    def return_protein_umap(input_1,value):
        
        df = pd.read_json(value)
        
        
        matrix = df.loc[:, '# Cells':'Cell_Identifier']
        del matrix['# Cells']
        del matrix['Cell_Identifier']
        
        
        df['Chamber Column'] =df['Chamber Column'].astype(str)
        df['Chamber Row'] = df['Chamber Row'].astype(str)
        df['Location'] = df['Chamber Column'].str.cat(df[['Chamber Row']], sep='|')



        # for Python 3 (the unicode type does not exist and is replaced by str)
        matrix.index = matrix.index.map(str)


        final = pd.DataFrame(np.nan_to_num(np.log(matrix.T), neginf=0)).T
        final.index = matrix.index
        final.columns = matrix.columns

        
        
        
        
        umap_3d = UMAP(n_components=2, init='random', random_state=0)
        proj_3d = umap_3d.fit_transform(final)
        
        
        test = pd.DataFrame(proj_3d)
        test['Location'] = df['Location']
        test['Column'] = df['Chamber Column']
        fig = px.scatter(
            test, x=0, y=1,hover_data=['Location'],height=800,width=1600,
            color=test.Column,labels={'color': 'Column'})
        fig.update_traces(marker_size=9)




        
        # fig = px.scatter(
        # proj_3d, x=0, y=1,height=800,width=1600,
        # color=df.Cell_Identifier, text=df.Location,labels={'color': 'Cell_Identifier'})
        # fig.update_traces(marker_size=9)
        
        # fig = px.scatter_3d(
        # proj_3d, x=0, y=1, z=2,height=800,width=1600,
        # color=df.Cell_Identifier, labels={'color': 'Cell_Identifier'})
        # fig.update_traces(marker_size=5)


        fig.update_layout(
            title={
            'text' : 'UMAP',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig.update_yaxes(autorange="reversed")
        fig.update_layout(
            font=dict(
                family="Arial",
                size=11,
            ))
        

        return(fig)
    
    
    
    
    @dashapp.callback(
        dash.dependencies.Output('intermediate-value-bead-line', 'figure'),
        [dash.dependencies.Input('button_bead', 'n_clicks')],
        [dash.dependencies.State('bead-dropdown', 'value')])
    def return_zipcode_nucleotideRC_scatterplot(n_clicks,value):
        
        bead_sequence = pd.read_csv(value+'/IntegratorResults/PC/longest_barcode_beads.txt',sep='\t')
        
        
        test = bead_sequence['rc_barcode'].apply(lambda x: pd.Series(list(x)))
        ctabs = {}
        for column in test.columns:
            summary = test[column].value_counts()
            summary = (summary/len(test))
            ctabs[column] = summary
        final_T = pd.DataFrame(ctabs).T            
        final_T['Position'] = final_T.index
        aggregate_df = []
        for column in ['A','C','G','T']:
            subset_columns = [column,'Position']
            final = final_T[subset_columns]
            final['Base'] = column
            final.columns = ['Frequency','Position','Base']
            aggregate_df.append(final)
        df = pd.concat(aggregate_df)
        df = df.fillna(0)
        
        fig = px.line(df, x="Position", y="Frequency", color='Base')
        #fig = px.pie(match_non_match_df, values='Value', names='Match')
        fig.update_layout(
            title={
            'text' : 'Zipcode Reverse Complement - Per Base Sequence',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig.update_layout(
            yaxis_range=[0,1],
            font=dict(
                family="Arial",
                size=17,
            ))
        
        return(fig)
        
    @dashapp.callback(
    dash.dependencies.Output('intermediate-value-bead-line-2', 'figure'),
    [dash.dependencies.Input('button_bead', 'n_clicks')],
    [dash.dependencies.State('bead-dropdown', 'value')])
    def return_zipcode_nucleotide_scatterplot(n_clicks,value):
        
        bead_sequence = pd.read_csv(value+'/IntegratorResults/PC/longest_barcode_beads.txt',sep='\t')
        
        
        test = bead_sequence['barcode'].apply(lambda x: pd.Series(list(x)))
        ctabs = {}
        for column in test.columns:
            summary = test[column].value_counts()
            summary = (summary/len(test))
            ctabs[column] = summary
        final_T = pd.DataFrame(ctabs).T            
        final_T['Position'] = final_T.index
        aggregate_df = []
        for column in ['A','C','G','T']:
            subset_columns = [column,'Position']
            final = final_T[subset_columns]
            final['Base'] = column
            final.columns = ['Frequency','Position','Base']
            aggregate_df.append(final)
        df = pd.concat(aggregate_df)
        df = df.fillna(0)

        fig = px.line(df, x="Position", y="Frequency", color='Base')
        #fig = px.pie(match_non_match_df, values='Value', names='Match')
        fig.update_layout(
            title={
            'text' : 'Zipcode Cycles - Per Base Sequence',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig.update_layout(
            yaxis_range=[0,1],
            font=dict(
                family="Arial",
                size=17,
            ))
        
        return(fig)
        
    @dashapp.callback(
    dash.dependencies.Output('bead_sequence_table', 'children'),
    [dash.dependencies.Input('button_bead', 'n_clicks')],
    [State('bead-dropdown','value')])
    def return_zipcode_particlecall_table(n_clicks, value):
        
        df = pd.read_csv(value+'/IntegratorResults/PC/longest_barcode_beads.txt',sep='\t')
        
        
       # print(df.head())
        table = html.Div([dash_table.DataTable(
            columns=[{"name": i, "id": i} for i in df.columns],
            style_table={
                'overflowX': 'scroll',
                'width':'100%',
                'margin':'auto'},
            data=df.to_dict('records'),
            virtualization=True,
            page_size=40,
            sort_action="native",
            style_header={'backgroundColor': 'rgb(30, 30, 30)'},

            style_filter={
                'backgroundColor': '#A5A3A6',
                'color': 'black'
            },  
            style_cell={
                'backgroundColor': 'rgb(0, 25, 51)',
                'color': 'white',
                'textAlign':'center',
                'fontSize':10, 
                'font-family':'sans-serif'},
            
            export_format='xlsx',
            export_headers='display',
            merge_duplicate_headers=True)]
        )
        
        
        return([table])
        
        
    @dashapp.callback(
        dash.dependencies.Output('intermediate-value-bead-scatter', 'figure'),
        [dash.dependencies.Input('button_bead', 'n_clicks')],
        [dash.dependencies.State('bead-dropdown', 'value')])
    def return_zipcode_chip_scatterplot(n_clicks,value):
        
        bead_sequence = pd.read_csv(value+'/IntegratorResults/PC/longest_barcode_beads.txt',sep='\t')
        
        
        fig = px.scatter(bead_sequence,x=bead_sequence['col'], y=bead_sequence['row'],hover_data=['rc_barcode'])
        #fig = px.pie(match_non_match_df, values='Value', names='Match')
        fig.update_layout(
            title={
            'text' : 'Bead Chip Location',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig['layout']['yaxis']['autorange'] = "reversed"
        fig.update_layout(
            yaxis_range=[0,255],
            font=dict(
                family="Arial",
                size=17,
            ))

        return(fig)
    

        #return(html_file_obj)
        
    @dashapp.callback(
        dash.dependencies.Output('intermediate-value-bead-number', 'children'),
        [dash.dependencies.Input('button_bead', 'n_clicks')],
        [dash.dependencies.State('bead-dropdown', 'value')])
    def return_zipcode_number(n_clicks,value):
        
        bead_sequence = pd.read_csv(value+'/IntegratorResults/PC/longest_barcode_beads.txt',sep='\t')
        
        bead_number = len(bead_sequence)
        
        string_value = 'Number of Beads:'+ str(bead_number)

        return(string_value)
    
    
    
    
    @dashapp.callback(
    dash.dependencies.Output('matched_barplot', 'figure'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_levenshtein_frequency_barlot(n_clicks,value):
        
        
        df = pd.read_csv(value+'/12_cycle/12.Cycle_Matching.Analysis.csv',index_col=0)
        
        mismatch_table = pd.DataFrame((df['Levenshtein_distance'].value_counts()))
        fig = px.bar(mismatch_table,x=mismatch_table.index,y='Levenshtein_distance')
        fig.update_layout(
            title={
            'text' : 'Bar Plot - Frequency of Number of Mismatches',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig.update_layout(
            font=dict(
                family="Arial",
                size=11,
            ))
        

        return(fig)
        
    @dashapp.callback(
    dash.dependencies.Output('matched_frequency_barplot', 'children'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_levenshtein_mismatch_html(n_clicks,value):
        
        
        p= value+'/12_cycle/mismatch_frequency/Full.matching12.Cycle_Matching.Analysis.csv.html'
        
        #mismatch_table = pd.DataFrame((df['Levenshtein_distance'].value_counts()))
        html_file = open(p,"r").read()

        
        html_file_obj = html.Iframe(srcDoc=html_file,style={"height": 510, "width": 1500}) 
        #html_file_obj = dash_dangerously_set_inner_html.DangerouslySetInnerHTML(html_file)
        return(html_file_obj)

     
        
     
    
    @dashapp.callback(
    dash.dependencies.Output('matched_table', 'children'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_levenshtein_frequency_table(n_clicks,value):
        
        
        df = pd.read_csv(value+'/12_cycle/12.Cycle_Matching.Analysis.csv',index_col=0)
        
        table = html.Div([dash_table.DataTable(
            columns=[{"name": i, "id": i} for i in df.columns],
            style_table={
                'overflowX': 'scroll',
                'width':'100%',
                'margin':'auto'},
            data=df.to_dict('records'),
            virtualization=True,
            page_size=40,
            sort_action="native",
            style_header={'backgroundColor': 'rgb(30, 30, 30)'},

            style_filter={
                'backgroundColor': '#A5A3A6',
                'color': 'black'
            },  
            style_cell={
                'backgroundColor': 'rgb(0, 25, 51)',
                'color': 'white',
                'textAlign':'center',
                'fontSize':10, 
                'font-family':'sans-serif'},
            
            export_format='xlsx',
            export_headers='display',
            merge_duplicate_headers=True)]
        )
        
        
        
        
        return(table)
    
    @dashapp.callback(
    dash.dependencies.Output('matched_zipcode_protein', 'children'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_protein_bead_table(n_clicks,value):
        
        
        df = pd.read_csv(value+'/12_cycle/12_Cycle_Protein_Bead.Matching.Analysis.csv',index_col=0)
        
        table = html.Div([dash_table.DataTable(
            columns=[{"name": i, "id": i} for i in df.columns],
            style_table={
                'overflowX': 'scroll',
                'width':'100%',
                'margin':'auto'},
            data=df.to_dict('records'),
            virtualization=True,
            page_size=40,
            sort_action="native",
            style_header={'backgroundColor': 'rgb(30, 30, 30)'},

            style_filter={
                'backgroundColor': '#A5A3A6',
                'color': 'black'
            },  
            style_cell={
                'backgroundColor': 'rgb(0, 25, 51)',
                'color': 'white',
                'textAlign':'center',
                'fontSize':10, 
                'font-family':'sans-serif'},
            
            export_format='xlsx',
            export_headers='display',
            merge_duplicate_headers=True)]
        )
        return([table])
    
    @dashapp.callback(
    dash.dependencies.Output('protein-zipcode-number', 'children'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_proteinzipcode_number(n_clicks,value):
        
        df = pd.read_csv(value+'/12_cycle/12_Cycle_Protein_Bead.Matching.Analysis.csv',index_col=0)
        
        bead_number = len(df)
        
        string_value = 'Number of Protein Chambers with Zipcodes:'+ str(bead_number)

        return(string_value)
    

    
    
    @dashapp.callback(
    dash.dependencies.Output('matched_final_summary_table', 'children'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_levenshtein_match_table(n_clicks,value):
        
        
        df = pd.read_csv(value+'/12_cycle/12_Cycles.Matched.Cells.csv',index_col=0)
        def get_subset(df,mismatch,identifier):
            df_subset_0 = (df[df['Levenshtein_distance'] <= mismatch])
            count = df_subset_0['Cell_Identifier'].value_counts()
            new_0 = pd.DataFrame(count).T
            new_0.index = [identifier]
            new_0['Total'] = new_0.sum(axis=1)
            return(new_0)
        
        perfect = get_subset(df,0,'Perfect')
        mismatch_1 = get_subset(df,1,'< One Mismatch')
        mismatch_2 = get_subset(df,2,'< Two Mismatch')
        
        
        
        df = (pd.concat([perfect, mismatch_1,mismatch_2]))
        df.insert (0, "Threshold", df.index)
        
        total_column = df.pop('Total')
        df = pd.concat([df, total_column], 1)
        


        table = html.Div([dash_table.DataTable(
            columns=[{"name": i, "id": i} for i in df.columns],

            data=df.to_dict('records'),
            
            sort_action="native",
            style_header={'backgroundColor': 'rgb(30, 30, 30)'},

            style_filter={
                'backgroundColor': '#A5A3A6',
                'color': 'black'
            },  
            style_cell={
                'backgroundColor': 'rgb(0, 25, 51)',
                'color': 'white',
                'textAlign':'center',
                'fontSize':25, 
                'font-family':'sans-serif'},
            fill_width=False,
            export_format='xlsx',
            export_headers='display',
            merge_duplicate_headers=True)]
        )
        
        return([table])
    
    
    
    @dashapp.callback(
    dash.dependencies.Output('matched_final_table', 'children'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_zipcode_protein_illumina_table(n_clicks,value):
        
        
        df = pd.read_csv(value+'/12_cycle/12_Cycles.Matched.Cells.csv',index_col=0)
        
        table = html.Div([dash_table.DataTable(
            columns=[{"name": i, "id": i} for i in df.columns],
            style_table={
                'overflowX': 'scroll',
                'width':'100%',
                'margin':'auto'},
            data=df.to_dict('records'),
            virtualization=True,
            page_size=40,
            sort_action="native",
            style_header={'backgroundColor': 'rgb(30, 30, 30)'},

            style_filter={
                'backgroundColor': '#A5A3A6',
                'color': 'black'
            },  
            style_cell={
                'backgroundColor': 'rgb(0, 25, 51)',
                'color': 'white',
                'textAlign':'center',
                'fontSize':10, 
                'font-family':'sans-serif'},
            
            export_format='xlsx',
            export_headers='display',
            merge_duplicate_headers=True)]
        )
        
        
        
        
        
        
        return(table)
    
    @dashapp.callback(
    dash.dependencies.Output('matched_zipcode_protein_location', 'figure'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_protein_zipcode_table(n_clicks,value):
        
        
        bead_sequence = pd.read_csv(value+'/12_cycle/12_Cycle_Protein_Bead.Matching.Analysis.csv',index_col=0)
        
        fig = px.scatter(bead_sequence,x=bead_sequence['column'], y=bead_sequence['row'],hover_data=['rc_barcode'])
        #fig = px.pie(match_non_match_df, values='Value', names='Match')
        fig.update_layout(
            title={
            'text' : 'Chambers with Zipcode Sequences & Protein Data',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig['layout']['yaxis']['autorange'] = "reversed"
        fig.update_layout(
            yaxis_range=[0,255],
            font=dict(
                family="Arial",
                size=17,
            ))

        return(fig)
    
    @dashapp.callback(
    dash.dependencies.Output('matched_final_scatter', 'figure'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_final_chip_scatter(n_clicks,value):
        
        
        bead_sequence = pd.read_csv(value+'/12_cycle/12_Cycles.Matched.Cells.csv',index_col=0)
        
        fig = px.scatter(bead_sequence,x=bead_sequence['column'], y=bead_sequence['row'],hover_data=['rc_barcode'])
        #fig = px.pie(match_non_match_df, values='Value', names='Match')
        fig.update_layout(
            title={
            'text' : 'Chambers with Zipcodes/Protein/RNA',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig['layout']['yaxis']['autorange'] = "reversed"
        fig.update_layout(
            yaxis_range=[0,255],
            font=dict(
                family="Arial",
                size=17,
            ))

        return(fig)
        
        
        
    @dashapp.callback(
    dash.dependencies.Output('protein_matched_barplot', 'figure'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [State('matched-dropdown','value')])
    def return_matched_final_proteinExp_barplot(n_clicks,value):
        
        df =  pd.read_csv(value+'/12_cycle/12_Cycles.Matched.Cells.csv',index_col=0)


        zero = get_box_plot(df,'Zero Cell')
        single  = get_box_plot(df,'Single Cell')
        two  = get_box_plot(df,'Two Cell')
        three  = get_box_plot(df,'Multi Cell')
        
        combined = (pd.concat([zero,single,two,three]))
        
        fig = px.bar(combined, x="Mean Expression", y="Protein",
             color='classification', barmode='group',orientation='h',
             height=800,width=1100)
        
        fig.update_layout(
            title={
            'text' : 'Matched Cells - Proteomic Mean Expression vs. Chamber Classification',
            'x':0.5,
            'xanchor': 'center'
        })
        fig.update_yaxes(tickfont_size=14)
        fig.update_yaxes(autorange="reversed")
        fig.update_layout(
            font=dict(
                family="Arial",
                size=17,
            ))
        

        return(fig)
    
    
    
    
    
    
    @dashapp.callback(
    dash.dependencies.Output('matched_table_2', 'children'),
    [dash.dependencies.Input('matched-button', 'n_clicks')],
    [dash.dependencies.State('matched-dropdown', 'value')])
    def return_matched_summary_table(n_clicks,value):
        
        def get_length(df):
            df_subset_0 = len(df[df['Levenshtein_distance'] <= 0])
            df_subset_1 = len(df[df['Levenshtein_distance'] <= 1])
            df_subset_2 = len(df[df['Levenshtein_distance'] <= 2])
            return(df_subset_0,df_subset_1,df_subset_2)
        


        csv = pd.read_csv(value+'/12_cycle/12.Cycle_Matching.Analysis.csv',index_col=0)
        df = pd.DataFrame(get_length(csv))
        df['Threshold'] = ['Perfect','<= 1 Mismatch','<= 2 Mismatch']
        df.columns = ['Number of Matches','Threshold']
        df = df[['Threshold','Number of Matches']]
        #df.index = ['Perfect','<= 1 Mismatch','<= 2 Mismatch']
        
        table = html.Div([dash_table.DataTable(
            columns=[{"name": i, "id": i} for i in df.columns],

            data=df.to_dict('records'),
            
            sort_action="native",
            style_header={'backgroundColor': 'rgb(30, 30, 30)'},

            style_filter={
                'backgroundColor': '#A5A3A6',
                'color': 'black'
            },  
            style_cell={
                'backgroundColor': 'rgb(0, 25, 51)',
                'color': 'white',
                'textAlign':'center',
                'fontSize':25, 
                'font-family':'sans-serif'},
            
            export_format='xlsx',
            export_headers='display',
            merge_duplicate_headers=True)]
        )
        
        
        
        
        return(table)
    
    
    ## Interval Updates
    
    
    @dashapp.callback(Output('my-dropdown', 'options'),
    [Input('interval_component', 'n_intervals')])
    def update_interval_transcriptome_sample(num):
        
        
        qc_list = []



        for root, dirs, files in os.walk('/home/isodev/Processed_Data/Illumina_RNA/'):
            for file in files:
                if file.endswith('5000_output.html'):
                    p = os.path.join(root)
                    drop_dict = {}

                    sample = p.split('/')[5]
                    #file_split_1 = file.split('_R2')[0]
                    #sample_name = file_split_1.split('QC.')[1]
                    drop_dict['label'] = sample
                    drop_dict['value'] = p

                    qc_list.append(drop_dict)
                    
        return(qc_list)
    
    @dashapp.callback(Output('bead-dropdown', 'options'),
    [Input('interval_component', 'n_intervals')])
    def update_interval_protein_sample(num):
        
        
        bead_list = []
        d= '/mnt/duomicsystemsruns/zipcode_results/'
        for path in os.listdir(d):
            full_path = os.path.join(d, path)
            drop_dict = {}
            drop_dict['label'] = full_path.split('/')[4]
            drop_dict['value'] = full_path
            bead_list.append(drop_dict)
                
                    
        return(bead_list)
    
    
    
    @dashapp.callback(Output('bead-matching-dropdown', 'options'),
    [Input('interval_component', 'n_intervals')])
    def update_interval_zipcode_sample(num):
        
        
    
        bead_list = []
        d= '/mnt/duomicsystemsruns/zipcode_results/'
        for path in os.listdir(d):
            full_path = os.path.join(d, path)
            drop_dict = {}
            drop_dict['label'] = full_path.split('/')[4]
            drop_dict['value'] = full_path
            bead_list.append(drop_dict)
                
                    
        return(bead_list)
    
    
    
    
    @dashapp.callback(Output('matching-report-dropdown', 'options'),
    [Input('interval_component', 'n_intervals')])
    def update_interval_matching_log(num):
        
        

        matching_log_path = '/home/isodev/Processed_Data/logs/matching'

        matching_report_list = []

        for root, dirs, files in os.walk(matching_log_path):
            for file in files:
                if file.endswith('.html'):
                    p = os.path.join(root,file)
                    drop_dict = {}


                    drop_dict['label'] = p.split('/')[-1]
                    drop_dict['value'] = p

                    matching_report_list.append(drop_dict)
        return(matching_report_list)
    
    
    
        
    @dashapp.callback(Output('matched-dropdown', 'options'),
    [Input('interval_component', 'n_intervals')])
    def update_interval_matched_sample(num):
        
        
        
        matched_data_list = []
        matched_folders= '/home/isodev/Processed_Data/Matched_Data/'
        for path in os.listdir(matched_folders):
            full_path = os.path.join(matched_folders, path)
            drop_dict = {}
            drop_dict['label'] = full_path.split('/')[5]
            drop_dict['value'] = full_path
            matched_data_list.append(drop_dict)
            
        return(matched_data_list)

    
 
    
    @dashapp.callback(Output('tmux_table', 'children'),
    [Input('interval_component_2', 'n_intervals')])
    def update_interval_active_pipelines(num):
        
        

        proc = subprocess.Popen(["tmux", "ls"],stdout=subprocess.PIPE)

        if sys.version_info[0] < 3: 
            from StringIO import StringIO
        else:
            from io import StringIO

        b = StringIO(proc.communicate()[0].decode('utf-8'))

        df = pd.read_csv(b, sep=",",header=None)
        
        df.columns = ['Active Pipelines']
       
        table = html.Div([dash_table.DataTable(
                columns=[{"name": i, "id": i} for i in df.columns],
                data=df.to_dict('records'),
                row_selectable='multi',
                selected_rows=[],
                page_size=15,
                editable=True,
                sort_action="native",
                style_header={'backgroundColor': 'rgb(0,128,0)'},

                style_filter={
                    'backgroundColor': '#A5A3A6',
                    'color': 'black'
                },  
                style_cell={
                    'backgroundColor': 'rgb(0, 25, 51)',
                    'color': 'white',
                    'textAlign':'center',
                    'fontSize':30, 
                    'font-family':'sans-serif'},
                export_format='xlsx',
                export_headers='display',
                merge_duplicate_headers=True)]
            )

        return([table])
    
    
    @dashapp.callback([
     dash.dependencies.Output('matching_output', 'children')],
    [dash.dependencies.Input('matching-button', 'n_clicks')],
    [dash.dependencies.State('user_project_name', 'value'),
     dash.dependencies.State('alevin-matching-dropdown', 'value'),
    dash.dependencies.State('protein-matching-dropdown', 'value'),
    dash.dependencies.State('bead-matching-dropdown', 'value'),
    dash.dependencies.State('bead-selection-matching-dropdown','value')])
    def execute_nextflow_matching_pipeline(n_clicks,user_project_name,alevin_value,protein_value,bead_value,basecaller_selection):
        
        
        if n_clicks:
            
            #print(user_project_name)
            #print(alevin_value)
            
            bead_value = bead_value + str(basecaller_selection)
            
            
            alevin_input_split = alevin_value.split('/')
            alevin_input_split.pop()
            alevin_input = ("/".join(alevin_input_split))

            #project_directory = value.split('/')[-1]
#             #project_id = project_directory.split('.')[0]
#             project_id_report = '/home/isodev/Processed_Data/logs/'+user_project_name + '.html'
            print('\n')
            print('\n')
            print('\n')

    
            print(alevin_input)
            print(bead_value)
            print(protein_value)
            print('\n')
            print('\n')
            print('\n')
            print('\n')
            project_id_report = '/home/isodev/Processed_Data/logs/matching/'+user_project_name + '.html'

            subprocess.run(["tmux", "new-session","-d" ,"-s", user_project_name])
            #subprocess.run(["tmux", "send-keys","conda activate duomic"])
            subprocess.run(["tmux", "send-keys","conda activate duomic","C-m"])
           # subprocess.run(["tmux", "send-keys","exec nextflow run /home/isodev/development/dev/matching/match.dashboard.nf --alevin_input {} --zipcode_input {} --protein_input {} --project_id {} -with-report {}".format(alevin_input,bead_value,protein_value,user_project_name,project_id_report),"C-m"])
            subprocess.run(["tmux", "send-keys","exec nextflow run /home/isodev/development/dev/matching/match.dashboard.nf --alevin_input {} --zipcode_input {} --protein_input {} --project_id {} -with-report {}".format(alevin_input,bead_value,protein_value,user_project_name,project_id_report),"C-m"])

            
            

            return([])   
        
    @dashapp.callback(
    dash.dependencies.Output('matching-report-html_image', 'children'),
    [dash.dependencies.Input('matching-report-button', 'n_clicks')],
    [dash.dependencies.State('matching-report-dropdown', 'value')])
    
    def return_matched_nextflow_report(n_clicks,value):
        

        
        
        print(value)
        
        html_file = open(value,"r").read()
        
        html_file_obj = html.Iframe(srcDoc=html_file,style={"height": 1500, "width": 1600}) 
        #html_file_obj = dash_dangerously_set_inner_html.DangerouslySetInnerHTML(html_file)
        return(html_file_obj)
    
    
    
    
    @dashapp.callback(
    dash.dependencies.Output('rna-analysis-data-dropdown', 'style'),
    [dash.dependencies.Input('analysis-button', 'n_clicks')],
    [dash.dependencies.Input('analysis-type-dropdown', 'value')])

    def show_hide_element(n_clicks,value):
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        print(value)
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        
        if value == 'Transcriptome':
            return {'display': 'block'}
        if value == 'Matched':
            return {'display': 'none'}

            
    
    @dashapp.callback(
    dash.dependencies.Output('matched-analysis-data-dropdown', 'style'),
    [dash.dependencies.Input('analysis-button', 'n_clicks')],
    [dash.dependencies.Input('analysis-type-dropdown', 'value')])

    def show_hide_element(n_clicks,value):
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        print(value)
        print('\n')
        print('\n')
        print('\n')
        print('\n')
        
        if value == 'Matched':
            return {'display': 'block'}
        if value == 'Transcriptome':
            return {'display': 'none'}


    
    @dashapp.callback(
    Output('intermediate-value-bead','style'),
    [Input('button_bead','n_clicks')])

    def show_hide_element_button(value):

        if value == None:
            return {'display': 'none'}
        else:
            return {'display': 'inline-block'}
        
        
    
    @dashapp.callback(
    Output('card-content','style'),
    [Input('card-tabs','active_tab')])

    def active_tab_callback(active_tab):
        print(active_tab)
        if active_tab == 'tab-1':
            return {'display': 'inline'}
        else:
            return {'display': 'none'}

    @dashapp.callback(
    Output('card-content5','style'),
    [Input('card-tabs','active_tab')])

    def active_tab_callback(active_tab):
        print(active_tab)
        if active_tab == 'tab-5':
            return {'display': 'inline'}
        else:
            return {'display': 'none'}
        
    @dashapp.callback(
    Output('card-content6','style'),
    [Input('card-tabs','active_tab')])

    def active_tab_callback(active_tab):
        print(active_tab)
        if active_tab == 'tab-6':
            return {'display': 'inline'}
        else:
            return {'display': 'none'}

        
    @dashapp.callback(
    Output('card-content2','style'),
    [Input('card-tabs','active_tab')])

    def tab_content(active_tab):
        if active_tab == 'tab-2':
            return {'display': 'inline'}
        else:
            return {'display': 'none'}
        
    @dashapp.callback(
    Output('card-content3','style'),
    [Input('card-tabs','active_tab')])

    def tab_content(active_tab):
        if active_tab == 'tab-3':
            return {'display': 'inline'}
        else:
            return {'display': 'none'}

        
    @dashapp.callback(
    Output('card-content4','style'),
    [Input('card-tabs','active_tab')])

    def tab_content(active_tab):
        if active_tab == 'tab-4':
            return {'display': 'inline'}
        else:
            return {'display': 'none'}
        
    
    @dashapp.callback(
    Output('alevin-tab-body','style'),
    [Input('button','n_clicks')])

    def show_hide_element_button(value):

        if value == None:
            return {'display': 'none'}
        else:
            return {'display': 'inline'}
        
    @dashapp.callback(
    Output('analyze-tab-body','style'),
    [Input('analysis-button','n_clicks')])

    def show_hide_element_button(value):

        if value == None:
            return {'display': 'none'}
        else:
            return {'display': 'inline'}
        
        
    @dashapp.callback(
    Output('analysis-tab-body','style'),
    [Input('button','n_clicks')])

    def show_hide_element_button(value):

        if value == None:
            return {'display': 'none'}
        else:
            return {'display': 'inline'}
        
    @dashapp.callback(
    Output('protein-tab-body','style'),
    [Input('button_protein','n_clicks')])

    def show_hide_element_button(value):

        if value == None:
            return {'display': 'none'}
        else:
            return {'display': 'inline'}
        
    @dashapp.callback(
    Output('bead-tab-body','style'),
    [Input('button_bead','n_clicks')])

    def show_hide_element_button(value):

        if value == None:
            return {'display': 'none'}
        else:
            return {'display': 'inline'}
     
    @dashapp.callback(
    Output('match-tab-body','style'),
    [Input('matched-button','n_clicks')])

    def show_hide_element_button(value):

        if value == None:
            return {'display': 'none'}
        else:
            return {'display': 'inline'}