import dash as dcc
from dash import html
import dash_bootstrap_components as dbc
import pandas as pd
from dash.dependencies import Input, Output
from functools import reduce
from functools import reduce
import scipy
import scipy.stats
from functools import reduce
import numpy 
import os
import pandas as pd
import dash_core_components as dcc
from dash.dependencies import Output, Input, State
import dash_table
from flask_login import current_user
from app.models import User, Post, Ownership
import base64


def get_layout():


    qc_list = []



    for root, dirs, files in os.walk('/home/isodev/Processed_Data/Illumina_RNA/'):
        for file in files:
            if file.endswith('default_output.html'):
                p = os.path.join(root)
                drop_dict = {}

                sample = p.split('/')[5]
                #file_split_1 = file.split('_R2')[0]
                #sample_name = file_split_1.split('QC.')[1]
                drop_dict['label'] = sample
                drop_dict['value'] = p

                qc_list.append(drop_dict)
                


                
#     bead_list = []

#     for root, dirs, files in os.walk('/home/isodev/Processed_Data/Zipcode/'):
#         for file in files:

#             p = os.path.join(root,file)
#             drop_dict = {}


#             drop_dict['label'] = p.split('/')[5]
#             drop_dict['value'] = p

#             bead_list.append(drop_dict)
                
    
    bead_list = []
    d= '/mnt/duomicsystemsruns/zipcode_results/'
    for path in os.listdir(d):
        full_path = os.path.join(d, path)
        drop_dict = {}
        drop_dict['label'] = full_path.split('/')[4]
        drop_dict['value'] = full_path
        bead_list.append(drop_dict)


    protein_list = []
    d= '/home/isodev/Processed_Data/Protein/'
    for path in os.listdir(d):
        full_path = os.path.join(d, path)
        drop_dict = {}
        drop_dict['label'] = full_path.split('/')[5]
        drop_dict['value'] = full_path
        protein_list.append(drop_dict)



    matched_data_list = []
    matched_folders= '/home/isodev/Processed_Data/Matched_Data/'
    for path in os.listdir(matched_folders):
        full_path = os.path.join(matched_folders, path)
        drop_dict = {}
        drop_dict['label'] = full_path.split('/')[5]
        drop_dict['value'] = full_path
        matched_data_list.append(drop_dict)
        

    
    d= '/home/isodev/Processed_Data/matching/logs'

    report_html_list = []

    for root, dirs, files in os.walk(d):
        for file in files:
            if file.endswith('.html'):
                p = os.path.join(root,file)
                drop_dict = {}


                drop_dict['label'] = p.split('/')[-1]
                drop_dict['value'] = p

                report_html_list.append(drop_dict)

                
                
        
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
                
                
    
    experiment_list = []
    d= '/home/isodev/Processed_Data/Experiment/'
    for path in os.listdir(d):
        full_path = os.path.join(d, path)
        drop_dict = {}
        drop_dict['label'] = full_path.split('/')[5]
        drop_dict['value'] = full_path
        experiment_list.append(drop_dict)            
                

    logo = os.getcwd() + '/app/static/images/isoplexis-Logo-225px_w.png' # replace with your own image
    #Docker
    #logo = '/app/app/static/images/isoplexis-Logo-225px_w.png' # replace with your own image
    encoded_image = base64.b64encode(open(logo, 'rb').read())

    search_bar = dbc.Row(
        [
            dbc.Col(dbc.Input(type="search", placeholder="Search")),
            dbc.Col(
                dbc.Button("Search", color="primary", className="ml-2"),
                width="auto",
            ),
        ],
        #no_gutters=True,
        className="ml-auto flex-nowrap mt-3 mt-md-0",
        align="center",
    )

    navbar = dbc.Navbar(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()), height="100")),
                    ],
                    align="center",
                    #no_gutters=True,
                ),
                href="/login",

            ),
            dbc.NavItem(dbc.NavLink("Home", active=True, href="/",external_link=True)),
            dbc.NavItem(dbc.NavLink("Logout", active=True, href="/logout",external_link=True)),
            dbc.NavItem(dbc.NavLink("Contact Us", active=True, href="",external_link=True)),
            dbc.NavbarToggler(id="navbar-toggler"),

            dbc.Collapse(search_bar, id="navbar-collapse", navbar=True),
        ],

        color="dark",
        dark=False,
    )



    jumbotron = dbc.Jumbotron(
    [
        dbc.Container(
            [dbc.Row([
                html.H2("Duomic Pipeline QC", className="display-3"),
            ], justify="center", align="center", className="h-50"),
             html.Hr(className="my-2"),
             dbc.Row([
                html.P("Quality Control Dashboard interfacing with Duomic Pipelines for Transcriptomic | Proteomics | and Matched Data"),
            ], justify="center", align="center", className="h-50"),
            dbc.Row([html.Div(id='tmux_table',style= { 'display': 'flex'},className="justify-content-center")], justify="center", align="center", className="h-50")],
    fluid=True)])
    
    
    #tmux_table =  dbc.Row([html.Div(id='tmux_table',style= { 'display': 'flex'},className="justify-content-center")], justify="center", align="center", className="h-50")
    
    # Tab 1 Alevin 
    alevin_button = html.Div(
        [
            dbc.Button("Submit", id="button",color="primary")
        ],style={'padding': 10}
    )



    alevin_dropdown_analysis = html.Div(dcc.Dropdown(
            id='my-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=qc_list,
             placeholder="Select RNA Sample",
        ),style={'padding': 10,'width':'25%','text-align':'center'})
    
    alevin_matching_dropdown = html.Div(dcc.Dropdown(
            id='alevin-matching-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=qc_list,
            placeholder="Select RNA Sample",
        ),style={'padding': 10,'width':'18%','text-align':'center'})
    
    
    
    
    alevin_scanpy_processed = html.Div(id="scanpy-qc-processed",style= { 'display': 'none'})
    
    alevin_scanpy = html.Div(dcc.Graph(id="scanpy-qc"),style= { 'display': 'flex'},className="justify-content-center")
    alevin_umap_top_expressed = html.Div(dcc.Graph(id="scanpy-umap-topgenes"),style= { 'display': 'flex'},className="justify-content-center")
    alevin_umap = html.Div(dcc.Graph(id="scanpy-umap"),style= { 'display': 'flex'},className="justify-content-center")
    
    alevin_umap_bar = dbc.Row([alevin_umap,alevin_umap_top_expressed], justify="center", align="center")
        
    #alevin_marker_table= html.Div(id ='marker_rank_table',style= { 'display': 'inline-block','textAlign' : 'center' },className="twelve columns")
    #alevin_marker_table = html.Div(id='marker_rank_table',style= { 'display': 'inline'},className="justify-content-center")
    #alevin_scanpy_go_table = html.Div(id='scanpy-qc-go-table',style= { 'display': 'flex'},className="justify-content-center")
    alevin_marker_table = html.Div([html.Div(id='marker_rank_table', style={'width': '99%'})],
                    style={'display': 'flex', 'justify-content': 'center'})
    
    
    
    alevin_scanpy_go_table = html.Div(id='scanpy-qc-go-table',style= { 'display': 'inline'},className='twelve columns')
    
    
    cell_barcode_number = dbc.Row([
                html.H2(id='intermediate-value-illumina-number', className="display-3"),
            ], justify="center", align="center", className="h-50")
    
    
    #dcc.Store
    
    alevin_scanpy_intermediate = html.Div(dcc.Store(id='intermediate-value-scanpy'), style={'display': 'none'})
    alevin_row = html.Div([alevin_dropdown_analysis,alevin_button],style=dict(display='flex'),className='justify-content-center')
    
    
    
    matched_button = html.Div(
        [
            dbc.Button("Submit", id="matched-button",color="primary")
        ],style={'padding': 10}
    )


    matched_data_dropdown = html.Div(dcc.Dropdown(
            id='matched-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=matched_data_list,
            placeholder="Select Matched Data",
        ),style={'padding': 10,'width':'25%','text-align':'center'})
    
    matched_data_row = html.Div([matched_data_dropdown,matched_button],style=dict(display='flex'),className='justify-content-center')
    
    
    
    #  Tab 2
    
    
    protein_button = html.Div(
        [
            dbc.Button("Submit", id="button_protein",color="primary")
        ],style={'padding': 10}
    )



    protein_input_dropdown = html.Div(dcc.Dropdown(
            id='protein-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=protein_list,
            placeholder="Select Protein Sample",
            #value=["test"]
        ),style={'padding': 10,'width':'25%','text-align':'center'})
    
    
    protein_matching_dropdown = html.Div(dcc.Dropdown(
            id='protein-matching-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=protein_list,
            placeholder="Select Protein Sample"
        ),style={'padding': 10,'width':'18%','text-align':'center'})
    
    
    protein_table = html.Div(id='protein_combined_table',style= { 'display': 'flex'},className="justify-content-center")
    protein_barplot = html.Div(dcc.Graph(id="protein_barplot"),style= { 'display': 'flex'},className="justify-content-center")
    protein_scatter = html.Div(dcc.Graph(id="protein_scatterplot"),style= { 'display': 'flex'},className="justify-content-center")
    protein_boxplot = html.Div(dcc.Graph(id="protein_boxplot"),style= { 'display': 'flex'},className="justify-content-center")
    protein_barplot_mean = html.Div(dcc.Graph(id="protein_barplot_mean"),style= { 'display': 'flex'},className="justify-content-center")
    
    
    #protein_tsne = html.Div(dcc.Graph(id="protein_tsne"),style= { 'display': 'flex'},className="justify-content-center")
    protein_umap = html.Div(dcc.Graph(id="protein_umap"),style= { 'display': 'flex'},className="justify-content-center")
    

    
    
    protein_hidden = html.Div(id='protein_table_hidden', style={'display': 'none'})
    matched_hidden = html.Div(id='matched_hidden', style={'display': 'none'})
     # Tab 3
        
    bead_button = html.Div(
        [
            dbc.Button("Submit", id="button_bead",color="primary")
        ],style={'padding': 10}
    )



    bead_barcode_files = html.Div(dcc.Dropdown(
            id='bead-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=bead_list,
            placeholder="Select Zipcode Sample",
        ),style={'padding': 10,'width':'25%','text-align':'center'})
    
    bead_matching_dropdown = html.Div(dcc.Dropdown(
            id='bead-matching-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=bead_list,
            placeholder="Select Zipcode Sample",
        ),style={'padding': 10,'width':'18%','text-align':'center'})

    bead_basecaller_selection_dropdown = html.Div(dcc.Dropdown(
            id='bead-selection-matching-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=[{'label':'Particle Call', 'value': "/IntegratorResults/PC/longest_barcode_beads.txt"}, {'label': 'Local Call', 'value': "/IntegratorResults/LC/longest_barcode_beads.txt"}],
            placeholder="Select BaseCaller",
        ),style={'padding': 10,'width':'18%','text-align':'center'})
    
    
    
    protein_row = html.Div([protein_input_dropdown,protein_button,protein_hidden],style=dict(display='flex'),className='justify-content-center')
    protein_row_2 = html.Div([protein_barplot,protein_scatter,protein_boxplot],style={'display':'flex'},className='justify-content-center')
    protein_row_3 = html.Div([protein_umap],style={'display':'flex'},className='justify-content-center')
        
    bead_row = html.Div([bead_barcode_files,bead_button],style=dict(display='flex'),className='justify-content-center')
                           
                        
    

    table = html.Div(id='output-data-upload')


    alevin_output_div = dbc.Row([html.Div(id='intermediate-value',style= { 'display': 'flex'}),
                                    html.Div(id='intermediate-value-2',style= { 'display': 'flex'})], justify="center", align="center")
    
    
    fastqc_output_div = dbc.Row([html.Div(id='intermediate-value-fastqc-r1',style= { 'display': 'flex'}),
                                    html.Div(id='intermediate-value-fastqc-r2',style= { 'display': 'flex'})], justify="center", align="center")
    
    fastqscreen_html = dbc.Row([html.Div(id='fastqscreen-html',style= { 'display': 'flex'})], justify="center", align="center")
    

    bead_scatterplot_analysis = html.Div(dcc.Graph(id="intermediate-value-bead-line"),style= { 'display': 'flex','textAlign': 'center'},className="six columns")
    bead_scatterplot_analysis_2 = html.Div(dcc.Graph(id="intermediate-value-bead-line-2"),style= { 'display': 'flex','textAlign': 'center'},className="six columns")
    bead_scatter_location = html.Div(dcc.Graph(id="intermediate-value-bead-scatter"),style= { 'display': 'flex','textAlign': 'center'},className="twelve columns")
    bead_scatter_table = html.Div(id='bead_sequence_table',style= { 'display': 'flex'},className="justify-content-center")
    
    bead_row_plot = dbc.Row([bead_scatterplot_analysis,bead_scatterplot_analysis_2,bead_scatter_location],justify="center", align="center")

    
    bead_number = dbc.Row([
                html.H2(id='intermediate-value-bead-number', className="display-3"),
            ], justify="center", align="center", className="h-50")
    

    # Interval 
    interval_update = html.Div(dcc.Interval(id='interval_component_2',
                interval=100000,
               
    ))
    
    interval_update_2 = html.Div(dcc.Interval(id='interval_component',
                interval=350000,
               
    ))
    
    
    # Pipeline body
    
    matching_button = html.Div(
        [
            dbc.Button("Submit", id="matching-button",color="success")
        ],style={'padding': 10}
    )
    
    matching_report_dropdown = html.Div(dcc.Dropdown(
            id='matching-report-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=matching_report_list,
            placeholder="Choose Pipeline QC Report"
        ),style={'padding': 10,'width':'25%','text-align':'center'})
    
    
    matching_report_button = html.Div(
        [
            dbc.Button("Submit", id="matching-report-button",color="primary")
        ],style={'padding': 10}
    )
    
    select_matching_report = dbc.Row(html.Div(html.H3("Select Pipeline Report"),style={'padding': 5,'margin-right': "15px"}), justify="center", align="center")
     
    matching_report_html = dbc.Row([html.Div(id='matching-report-html_image',style= { 'display': 'flex'})], justify="center", align="center")
    
    
    
    matching_report_row = html.Div([select_matching_report,matching_report_dropdown, matching_report_button],style=dict(display='flex'),className='justify-content-center')
    
    user_input = html.Div(dcc.Input(id="user_project_name", type="text", placeholder="", style={'marginRight':'10px','verticalAlign':'middle'}))
    
    matching_div = dbc.Row([html.Div(id='matching_output',style= { 'display': 'flex'})], justify="center", align="center")
    
        
    
    matching_pipeline_row = html.Div([alevin_matching_dropdown,protein_matching_dropdown,bead_matching_dropdown,bead_basecaller_selection_dropdown],style=dict(display='flex'),className='justify-content-center')
    
    name_matching = dbc.Row(html.Div(html.H3("Name Matching Analysis"),style={'padding': 5,'margin-right': "15px"}), justify="center", align="center")
    
    
    matching_pipeline_row_2 = html.Div([name_matching,user_input, matching_button],style=dict(display='flex'),className='justify-content-center')
    
    matching_analysis_header = dbc.Row(html.Div(html.H2("Select Datasets to Match"),style={'padding': 5}), justify="center", align="center")
    
    
    matching_body = html.Div([matching_analysis_header,matching_pipeline_row,matching_div,matching_pipeline_row_2,matching_report_row,matching_report_html],className='justify-content-center')
    
    
    
    rna_show = html.Div([alevin_scanpy_processed,alevin_scanpy,alevin_umap_bar,alevin_marker_table,alevin_output_div,fastqscreen_html,fastqc_output_div,table],id='alevin-tab-body',className='row',style={'display': 'none'})
    rna_body = html.Div([alevin_row,cell_barcode_number,rna_show])
    
    protein_show = html.Div([protein_row_2,protein_row_3,protein_barplot_mean,protein_table],id='protein-tab-body',className='row',style={'display': 'none'})
    protein_body = html.Div([protein_row,protein_show],className='justify-content-center')
    
    
    
    
    bead_show = html.Div([bead_number,bead_row_plot,bead_scatter_table],id='bead-tab-body',className='row',style={'display': 'none'})
    bead_body = html.Div([bead_row,bead_show])
    #layout = html.Div([navbar,jumbotron,dropdown_analysis,button,hidden_div,table])

    
    matched_freq_barplot = dbc.Row([html.Div(id='matched_frequency_barplot',style= { 'display': 'flex'})], justify="center", align="center")
    matched_barplot = html.Div(dcc.Graph(id="matched_barplot"),style= { 'display': 'flex'},className="justify-content-center")
    matched_table = html.Div(id='matched_table',style= { 'display': 'flex'},className="justify-content-center")
    matched_table_2 = html.Div(id='matched_table_2',style= { 'display': 'flex'},className="justify-content-center")
    matched_row_2 = dbc.Row([html.Div([matched_table,matched_barplot,matched_table_2],style={'display':'flex'},className='justify-content-center')],justify="center", align="center")
    
    protein_barplot_mean = html.Div(dcc.Graph(id="protein_matched_barplot"),style= { 'display': 'flex'},className="justify-content-center")
    
    
    matched_final_table = html.Div(id='matched_final_table',className="justify-content-center")
    matched_final_table_summary = dbc.Row([html.Div(id='matched_final_summary_table',className="justify-content-center")],justify="center", align="center")
    #matched_final_row = dbc.Row([html.Div([matched_final_table,matched_final_table_summary],style={'display':'flex'},className='justify-content-center')],justify="center", align="center")
    #matched_final_row = dbc.Row([html.Div([matched_final_table,matched_final_table_summary],style={'display':'inline-block'},className='justify-content-center')],justify="center", align="center")
    matched_final_scatter = html.Div(dcc.Graph(id="matched_final_scatter"),style= { 'display': 'flex'},className="justify-content-center")
    
    matched_final_row= html.Div([matched_final_table,matched_final_scatter],style={'display':'in-line'},className='justify-content-center')
    
    
    
    zipcoder_vs_illumina_header = dbc.Row(html.Div(html.H2("Zipcoder Sequences vs. Illumina Sequences"),style={'padding': 5}), justify="center", align="center")
    zipcoder_vs_protein_header = dbc.Row(html.Div(html.H2("Zipcoder Chambers vs. Protein Chambers"),style={'padding': 5}), justify="center", align="center")
    zipcoder_vs_illumina_vs_protein = dbc.Row(html.Div(html.H2("Zipcoder vs. Illumina vs. Protein"),style={'padding': 5}), justify="center", align="center")
    
    # Zipcode vs. Protein
    matched_zipcode_vs_protein_table = html.Div(id='matched_zipcode_protein',style= { 'display': 'flex'},className="justify-content-center")
    matched_zipcode_vs_protein_scatter = html.Div(dcc.Graph(id="matched_zipcode_protein_location"),style= { 'display': 'flex'},className="justify-content-center")
    
    protein_zipcode_match_number = dbc.Row([
                html.H2(id='protein-zipcode-number', className="display-3"),
            ], justify="center", align="center", className="h-50")
    
    
    matched_row_3 = html.Div([protein_zipcode_match_number,matched_zipcode_vs_protein_table,matched_zipcode_vs_protein_scatter],style={'display':'in-line'},className='justify-content-center')
    
    matched_show = html.Div([zipcoder_vs_illumina_vs_protein,matched_final_table_summary,protein_barplot_mean,matched_final_row,zipcoder_vs_illumina_header,matched_row_2,matched_freq_barplot,zipcoder_vs_protein_header,matched_row_3],id='match-tab-body',className='row',style={'display': 'none'})
    matched_body = html.Div([matched_data_row,matched_show],className='justify-content-center')

    
    # Analysis Tab
    
    
    
    analysis_type_dropdown = html.Div(dcc.Dropdown(
            id='analysis-type-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=[{'label': 'Transcriptome', 'value': 'Transcriptome'},
            {'label': 'Matched', 'value': 'Matched'}],
            placeholder="Choose Datatype to Analyze"
        ),style={'padding': 10,'width':'25%','text-align':'center'})
    
    
    rna_analysis_dropdown = html.Div(dcc.Dropdown(
            id='rna-analysis-data-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=qc_list,
            placeholder="Choose Transcriptome Input Data",
            multi=True
        ),style={'padding': 10,'width':'25%','text-align':'center'})
    
       
    experiment_analysis_dropdown = html.Div(dcc.Dropdown(
            id='matched-analysis-data-dropdown',
            #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
            options=experiment_list,
            placeholder="Choose Matched Data",
            multi=True
        ),style={'padding': 10,'width':'25%','text-align':'center'})
    
    analyze_button = html.Div(
        [
            dbc.Button("Submit", id="analysis-button",color="success")
        ],style={'padding': 10}
    )
    
    rna_analysis_row = html.Div([analysis_type_dropdown,analyze_button],style=dict(display='flex'),className='justify-content-center')
    
    rna_analysis_row_2 = html.Div([rna_analysis_dropdown],style=dict(display='flex'),className='justify-content-center')
    rna_analysis_row_3= html.Div([experiment_analysis_dropdown],style=dict(display='flex'),className='justify-content-center')
    #rna_analysis_row_4 = html.Div([rna_analysis_row_2,rna_analysis_row_3],style=dict(display='flex'),className='justify-content-center')
    #na_analysis_row_2 = dbc.Row([html.Div([rna_analysis_dropdown],style={'display':'flex'},className='justify-content-center')],justify="center", align="center")
    
    
    rna_show = html.Div([rna_analysis_row_2,rna_analysis_row_3],id='analyze-tab-body',className='row',style={'display': 'none'})
    
    rna_analysis_body = html.Div([rna_analysis_row,rna_show,alevin_scanpy_intermediate],className='justify-content-center')
    
    
    

    tab_layouts = dbc.Card(
    [
        dbc.CardHeader(
            dbc.Row(dbc.Tabs(
                [
                    dbc.Tab(label="Illumina RNA-seq Summary", tab_id="tab-1",label_style={"color": "#00AEF9"}),
                    dbc.Tab(label="IsoSpeak Protein Summary", tab_id="tab-2",label_style={"color": "#00AEF9"}),
                    dbc.Tab(label="Zipcoder Sequence Summary", tab_id="tab-3",label_style={"color": "#00AEF9"}),
                    dbc.Tab(label="Matching Analysis Pipeline", tab_id="tab-5",label_style={"color": "#d15406"}),
                    dbc.Tab(label="Matched Data Summary", tab_id="tab-4",label_style={"color": "#00AEF9"}),
                    dbc.Tab(label="Analysis Pipelines", tab_id="tab-6",label_style={"color": "#d15406"})                  
                ],
                id="card-tabs",
                card=True,
                active_tab="tab-1",
            ),justify="center", align="center", className="h-50")
        ),
        dbc.CardBody(rna_body,id='card-content',style={'display': 'none'}),
        dbc.CardBody(protein_body,id='card-content2',style={'display': 'none'}),
        dbc.CardBody(bead_body,id='card-content3',style={'display': 'none'}),
        dbc.CardBody(matching_body,id='card-content5',style={'display': 'none'}),
        dbc.CardBody(matched_body,id='card-content4',style={'display': 'none'}),
        dbc.CardBody(rna_analysis_body,id='card-content6',style={'display': 'none'})
        
        
    ]
)




    layout = html.Div([navbar,jumbotron,tab_layouts,interval_update,interval_update_2])
    
    
    
    	# return(layout)
 
    return(layout)

layout = (get_layout())
