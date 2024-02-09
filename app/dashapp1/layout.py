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



SIDEBAR_STYLE = {
        "position": "fixed",
        "top": 62.5,
        "left": 0,
        "bottom": 0,
        "width": "16rem",
        "height": "100%",
        "z-index": 1,
        "overflow-x": "hidden",
        "transition": "all 0.5s",
        "padding": "0.5rem 1rem",
        "background-color": "#f8f9fa",
    }

SIDEBAR_HIDEN = {
    "position": "fixed",
    "top": 62.5,
    "left": "-16rem",
    "bottom": 0,
    "width": "16rem",
    "height": "100%",
    "z-index": 1,
    "overflow-x": "hidden",
    "transition": "all 0.5s",
    "padding": "0rem 0rem",
    "background-color": "#f8f9fa",
}

# the styles for the main content position it to the right of the sidebar and
# add some padding.
CONTENT_STYLE = {
    "transition": "margin-left .5s",
    "margin-left": "18rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}

CONTENT_STYLE1 = {
    "transition": "margin-left .5s",
    "margin-left": "2rem",
    "margin-right": "2rem",
    "padding": "2rem 1rem",
    "background-color": "#f8f9fa",
}


def get_layout():

    
    navbar = dbc.NavbarSimple(
    children=[
        dbc.Button("Toggle", outline=True, color="secondary", className="mr-1", id="btn_sidebar"),
    ],
    brand="single cell RNA-seq",
    brand_href="#",
    color="dark",
    dark=True,
    fluid=True,
)


    # the style arguments for the sidebar. We use position:fixed and a fixed width
    
    sidebar = html.Div(
        [
            html.H2("Sidebar", className="display-4"),
            html.Hr(),
            html.P(
                "A simple sidebar layout with navigation links", className="lead"
            ),
            dbc.Nav(
                [
                    dbc.NavLink("Page 1", href="/page-1", id="page-1-link"),
                    dbc.NavLink("Page 2", href="/page-2", id="page-2-link"),
                    dbc.NavLink("Page 3", href="/page-3", id="page-3-link"),
                ],
                vertical=True,
                pills=True,
            ),
        ],
        id="sidebar",
        style=SIDEBAR_STYLE,
    )

    content = html.Div(

        id="page-content",
        style=CONTENT_STYLE)





    layout = html.Div(
    [
        dcc.Store(id='side_click'),
        dcc.Location(id="url"),
        navbar,
        sidebar,
        content,
    ],
)


    #layout = html.Div([navbar,jumbotron,tab_layouts,interval_update,interval_update_2])
    
    
    
    	# return(layout)
 
    return(layout)

layout = (get_layout())


#     qc_list = []



#     for root, dirs, files in os.walk('/'):
#         for file in files:
#             if file.endswith('default_output.html'):
#                 p = os.path.join(root)
#                 drop_dict = {}

#                 sample = p.split('/')[5]
#                 #file_split_1 = file.split('_R2')[0]
#                 #sample_name = file_split_1.split('QC.')[1]
#                 drop_dict['label'] = sample
#                 drop_dict['value'] = p

#                 qc_list.append(drop_dict)
                

                
#     logo = os.getcwd() + '/app/static/images/isoplexis-Logo-225px_w.png' # replace with your own image
#     #Docker
#     #logo = '/app/app/static/images/isoplexis-Logo-225px_w.png' # replace with your own image
#     encoded_image = base64.b64encode(open(logo, 'rb').read())

#     search_bar = dbc.Row(
#         [
#             dbc.Col(dbc.Input(type="search", placeholder="Search")),
#             dbc.Col(
#                 dbc.Button("Search", color="primary", className="ml-2"),
#                 width="auto",
#             ),
#         ],
#         #no_gutters=True,
#         className="ml-auto flex-nowrap mt-3 mt-md-0",
#         align="center",
#     )

#     navbar = dbc.Navbar(
#         [
#             html.A(
#                 # Use row and col to control vertical alignment of logo / brand
#                 dbc.Row(
#                     [
#                         dbc.Col(html.Img(src='data:image/png;base64,{}'.format(encoded_image.decode()), height="100")),
#                     ],
#                     align="center",
#                     #no_gutters=True,
#                 ),
#                 href="/login",

#             ),
#             dbc.NavItem(dbc.NavLink("Home", active=True, href="/",external_link=True)),
#             dbc.NavItem(dbc.NavLink("Logout", active=True, href="/logout",external_link=True)),
#             dbc.NavItem(dbc.NavLink("Contact Us", active=True, href="",external_link=True)),
#             dbc.NavbarToggler(id="navbar-toggler"),

#             dbc.Collapse(search_bar, id="navbar-collapse", navbar=True),
#         ],

#         color="dark",
#         dark=False,
#     )



#     jumbotron = dbc.Jumbotron(
#     [
#         dbc.Container(
#             [dbc.Row([
#                 html.H2("BioDash Analysis Tools", className="display-3"),
#             ], justify="center", align="center", className="h-50"),
#              html.Hr(className="my-2"),
#              dbc.Row([
#                 html.P(""),
#             ], justify="center", align="center", className="h-50")],
#     fluid=True)])
    
    
#     #tmux_table =  dbc.Row([html.Div(id='tmux_table',style= { 'display': 'flex'},className="justify-content-center")], justify="center", align="center", className="h-50")
    
#     # Tab 1 Alevin 
#     alevin_button = html.Div(
#         [
#             dbc.Button("Submit", id="button",color="primary")
#         ],style={'padding': 10}
#     )



#     alevin_dropdown_analysis = html.Div(dcc.Dropdown(
#             id='my-dropdown',
#             #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
#             options=qc_list,
#              placeholder="Select RNA Sample",
#         ),style={'padding': 10,'width':'25%','text-align':'center'})
    
#     alevin_matching_dropdown = html.Div(dcc.Dropdown(
#             id='alevin-matching-dropdown',
#             #options=[{'label':'OLFR2 vs. GPI1', 'value': "'OLFR2','GPI1'"}, {'label': 'GPI1 vs. OLFR2', 'value': "'GPI1','OLFR2'"}],
#             options=qc_list,
#             placeholder="Select RNA Sample",
#         ),style={'padding': 10,'width':'18%','text-align':'center'})
    
    
    
    
#     alevin_scanpy_processed = html.Div(id="scanpy-qc-processed",style= { 'display': 'none'})
    
#     alevin_scanpy = html.Div(dcc.Graph(id="scanpy-qc"),style= { 'display': 'flex'},className="justify-content-center")
#     alevin_umap_top_expressed = html.Div(dcc.Graph(id="scanpy-umap-topgenes"),style= { 'display': 'flex'},className="justify-content-center")
#     alevin_umap = html.Div(dcc.Graph(id="scanpy-umap"),style= { 'display': 'flex'},className="justify-content-center")
    
#     alevin_umap_bar = dbc.Row([alevin_umap,alevin_umap_top_expressed], justify="center", align="center")

#     alevin_marker_table = html.Div([html.Div(id='marker_rank_table', style={'width': '99%'})],
#                     style={'display': 'flex', 'justify-content': 'center'})
    
    
    
#     alevin_scanpy_go_table = html.Div(id='scanpy-qc-go-table',style= { 'display': 'inline'},className='twelve columns')
    
    
    
#     #dcc.Store
    
#     alevin_scanpy_intermediate = html.Div(dcc.Store(id='intermediate-value-scanpy'), style={'display': 'none'})
#     alevin_row = html.Div([alevin_dropdown_analysis,alevin_button],style=dict(display='flex'),className='justify-content-center')
    
    
       
    

#     table = html.Div(id='output-data-upload')


#     alevin_output_div = dbc.Row([html.Div(id='intermediate-value',style= { 'display': 'flex'}),
#                                     html.Div(id='intermediate-value-2',style= { 'display': 'flex'})], justify="center", align="center")
    
    
#     fastqc_output_div = dbc.Row([html.Div(id='intermediate-value-fastqc-r1',style= { 'display': 'flex'}),
#                                     html.Div(id='intermediate-value-fastqc-r2',style= { 'display': 'flex'})], justify="center", align="center")
    
#     fastqscreen_html = dbc.Row([html.Div(id='fastqscreen-html',style= { 'display': 'flex'})], justify="center", align="center")
    

#     bead_scatterplot_analysis = html.Div(dcc.Graph(id="intermediate-value-bead-line"),style= { 'display': 'flex','textAlign': 'center'},className="six columns")
#     bead_scatterplot_analysis_2 = html.Div(dcc.Graph(id="intermediate-value-bead-line-2"),style= { 'display': 'flex','textAlign': 'center'},className="six columns")
#     bead_scatter_location = html.Div(dcc.Graph(id="intermediate-value-bead-scatter"),style= { 'display': 'flex','textAlign': 'center'},className="twelve columns")
#     bead_scatter_table = html.Div(id='bead_sequence_table',style= { 'display': 'flex'},className="justify-content-center")
    
#     bead_row_plot = dbc.Row([bead_scatterplot_analysis,bead_scatterplot_analysis_2,bead_scatter_location],justify="center", align="center")

    
#     bead_number = dbc.Row([
#                 html.H2(id='intermediate-value-bead-number', className="display-3"),
#             ], justify="center", align="center", className="h-50")
    

#     # Interval 
#     interval_update = html.Div(dcc.Interval(id='interval_component_2',
#                 interval=100000,
               
#     ))
    
#     interval_update_2 = html.Div(dcc.Interval(id='interval_component',
#                 interval=350000,
               
#     ))
    
    
    
    
#     rna_show = html.Div([alevin_scanpy_processed,alevin_scanpy,alevin_umap_bar,alevin_marker_table,alevin_output_div,fastqscreen_html,fastqc_output_div,table],id='alevin-tab-body',className='row',style={'display': 'none'})
#     rna_body = html.Div([alevin_row,rna_show])
    
    
    

#     tab_layouts = dbc.Card(
#     [
#         dbc.CardHeader(
#             dbc.Row(dbc.Tabs(
#                 [
#                     dbc.Tab(label="Illumina RNA-seq Summary", tab_id="tab-1",label_style={"color": "#00AEF9"})               
#                 ],
#                 id="card-tabs",
#                 card=True,
#                 active_tab="tab-1",
#             ),justify="center", align="center", className="h-50")
#         ),
#         dbc.CardBody(rna_body,id='card-content',style={'display': 'none'})
        
        
#     ]
# )




#     layout = html.Div([navbar,jumbotron,tab_layouts,interval_update,interval_update_2])
    
    
    
#     	# return(layout)
 
    #return(layout)

#layout = (get_layout())
