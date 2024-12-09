# -*- coding: utf-8 -*-
"""
Created on Thu Nov  9 13:35:49 2023

@author: esra.bardakci
"""

import os
from Bio import SeqIO
import numpy as np
from scipy.signal import savgol_filter
from scipy.signal import find_peaks
from scipy import stats
import plotly.express as px
import pandas as pd
from itertools import combinations
from scipy.stats import pearsonr
import dash
import dash_core_components as dcc
import dash_bootstrap_components as dbc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import plotly.graph_objs as go
import plotly.offline as pyo
import plotly.graph_objs as go
from dash import Dash, dash_table, dcc, html
from random import randint
from plotly.subplots import make_subplots
from sklearn.linear_model import Ridge
from sklearn.linear_model import RidgeCV
import base64
import os
import re
from urllib.parse import quote as urlquote
from flask import Flask, send_from_directory


HeliksSizeS2 = [39, 43, 50, 67, 77, 87, 96, 106]
HeliksSizeS = [39, 43, 50, 67, 77, 87, 96, 106, 136]
HeliksSizeL=[39,49,67,77,107,135,153,163,202,235,300,360,398,450,487,496]
Liz120=[15,20,25,35,50,62,80,110,120]
Liz500=[35,50,75,100,139,150,160,200,300,340,350,400,450,490,500]

size_dict={'HeliksSizeS':HeliksSizeS,
           'HeliksSizeS2':HeliksSizeS2,
           'HeliksSizeL':HeliksSizeL,
           'Liz120':Liz120,
           'Liz500':Liz500}

dye_dict = {"1": 'blue',
            "2": 'green',
            "3": 'black',
            "4": 'red',
            "5": 'orange'}


def read_sizes(liz_path):
    liz_list=[]
    with open(liz_path,'r') as f:
        for line in f:
            for word in line.strip().split(','):
                try:
                    liz_list.append(int(word))
                except:
                    continue
    return liz_list



def LizStd(filepath,sizes,threshold):
    # size ladder datasının okunması
    data = SeqIO.read(filepath, "abi")
    firstAnnotatedData = data.annotations["abif_raw"]
    data205 = firstAnnotatedData["DATA205"]



    data205_up = np.array(data205)
    # dar bir scale kullanılarak orijinal data üzerinde smoothing uygulandı ve peakler hesaplandı.
    data205_up2 = savgol_filter(data205_up, 7, 2)
    data205_up2[data205_up2 <= threshold] = 0
    peaks, _ = find_peaks(data205_up2, prominence=1, distance=30)

    def find_best_correlation(S, P, tol):
        n = len(S)
        result_p = [subset for i, subset in enumerate(
            combinations(P, n)) if pearsonr(subset, S)[0] > tol]
        result_p_rho = [pearsonr(subset, S)[0] for i, subset in enumerate( 
            combinations(P, n)) if pearsonr(subset, S)[0] > tol]
        while len(result_p)==0 and tol >=0.1:
            tol-=0.001
            print(tol)
            result_p = [subset for i, subset in enumerate(
            combinations(P, n)) if pearsonr(subset, S)[0] > tol] 
            result_p_rho = [pearsonr(subset, S)[0] for i, subset in enumerate( 
            combinations(P, n)) if pearsonr(subset, S)[0] > tol]
            if tol <0.1:
                print('tol value is smaller than 0.1')
                break
        
        return result_p[result_p_rho.index(max(result_p_rho))],tol

    # LIZ500 size standard
    #sizes = [39, 43, 50, 67, 77, 87, 96, 106, 136]
    S = np.array(sizes)
    P = np.array(peaks)
    # find true peaks of ladder
    selected_peaks = list(find_best_correlation(S, P, tol=1)[0])
    #ladder_bases = pd.DataFrame({'peak_index':matched_peaks,'base':sizes})

    remained_peaks = list(set(peaks) ^ set(selected_peaks))

    sizepair = {}
    count = 0
    for i in peaks:
        if i in selected_peaks:
            # print('True')
            sizepair[i] = sizes[count]
            count += 1
        if i in remained_peaks:
            # print('False')
            sizepair[i] = 'no match' 

    return sizepair, peaks, data205_up2,find_best_correlation(S, P, tol=1)[1]


########################################################################

def read_bins(bin_file,panel_name):

    file = open(bin_file, 'r')
    panel_lines = []
    marker_lines = []
    counter = -1
    for line in file:
        counter += 1
        if 'Panel' in line:
            panel_lines.append(counter)
        if 'Marker Name' in line:
            marker_lines.append(counter)

    panel_dict = {}
    for i in panel_lines:
        with open(bin_file, 'r') as fp:
            x = fp.readlines()[i].strip().split('\t')[1]
            panel_dict[x] = i

    updated_marker_lines = []
    for i, j in enumerate(panel_lines):
        if j == panel_dict[panel_name]:
            try:
                updated_marker_lines = list(filter(
                    lambda marker_lines: marker_lines > panel_lines[i] and marker_lines < panel_lines[i+1], marker_lines))
                #print('1. ',updated_marker_lines)
            except:
                updated_marker_lines = list(
                    filter(lambda marker_lines: marker_lines > panel_lines[i], marker_lines))
                #print('2. ',updated_marker_lines)

    marker_list = []
    all_dataframes = []

    # print('updated_marker_lines')
    # print(updated_marker_lines)
    #print('panel_lines: ',panel_lines)
    #print('panel_dict: ',panel_dict)

    for i in range(len(updated_marker_lines)):
        bin_infos = []
        if i != len(updated_marker_lines)-1:
            # print('i',i)
            with open(bin_file, 'r') as fp:
                x = fp.readlines()[updated_marker_lines[i]
                                 :updated_marker_lines[i+1]]
                a = [a.split('\t') for a in x]
                marker = a[0][1].strip('\n')
                marker_list.append(marker)
                for i in a[2:]:
                    # print('*******1*******')
                    #print('marker: ',marker)
                    if len(i) == 4:
                        i[3] = i[3].strip('\n')
                        # print(i)
                        bin_infos.append(i)
                locals()[marker] = pd.DataFrame(bin_infos, columns=[
                    'allele', 'start_size', 'end_size', 'color'])
                # print(locals()[marker])

        else:
            # print('ii',i)
            # print('updated_marker_lines[i]',updated_marker_lines[i])
            b = [x for x in panel_lines if x > updated_marker_lines[i]]
            #print('b: ',b)
            #print('b[0]: ',b[0])

            if len(b) != 0:
                with open(bin_file, 'r') as fp:
                    x = fp.readlines()[updated_marker_lines[i]:b[0]]
                    a = [a.split('\t') for a in x]
                    #print('a: ')
                    # print(a)
                    marker = a[0][1].strip('\n')
                    marker_list.append(marker)
                    for i in a[2:]:
                        #print('i: ')
                        # print(i)
                        # print('******2********')
                        #print('marker: ',marker)
                        if len(i) == 4:
                            i[3] = i[3].strip('\n')
                            # print(i)
                            bin_infos.append(i)
                    locals()[marker] = pd.DataFrame(bin_infos, columns=[
                        'allele', 'start_size', 'end_size', 'color'])
                    # print(locals()[marker])

            else:
                with open(bin_file, 'r') as fp:
                    x = fp.readlines()[updated_marker_lines[i]:]
                    a = [a.split('\t') for a in x]
                    marker = a[0][1].strip('\n')
                    marker_list.append(marker)
                    for i in a[2:]:
                        # print('******2********')
                        #print('marker: ',marker)
                        if len(i) == 4:
                            i[3] = i[3].strip('\n')
                            # print(i)
                            bin_infos.append(i)
                    locals()[marker] = pd.DataFrame(bin_infos, columns=[
                        'allele', 'start_size', 'end_size', 'color'])

        all_dataframes.append(locals()[marker])

    return marker_list, all_dataframes

########################################################################

def write_bin_list(binpath,panel_name):
    #bin file okunurken marker sayısını kaydeder ki app call yaparken input olarak verebilelim.
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project")
    if not os.path.exists(UPLOAD_DIRECTORY):
        os.makedirs(UPLOAD_DIRECTORY) 
    fpath=os.path.join(UPLOAD_DIRECTORY, 'bin_list.txt')  
    
    bins = read_bins(binpath,panel_name)    
    bin_list = []
    for i in range(len(bins[1])):
        locals()['df-{}'.format(i)] = bins[1][i]
        bin_list.append(locals()['df-{}'.format(i)])
        
    f = open(fpath,"w")
    f.write((str(len(bin_list))))
    f.close()
    
#--------------------------------------------------#

def read_bin_list():
    #read from file to get length of bin_list which indicates count of bin tables
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project")
    fpath=os.path.join(UPLOAD_DIRECTORY, 'bin_list.txt')  
    file = open(fpath, 'r')  
    len_bin_list=file.read()
    file.close()
    final=re.findall('\d+',len_bin_list)
    
    if len_bin_list.isdigit():
        #print('---f----lenbinlist',final[0])
        return int(len_bin_list)

    

########################################################################


def get_sizes(filepath, size_keys, size_values, threshold):

    def read_fsa(filepath, labels=['DATA9', 'DATA10', 'DATA11', 'DATA12']):
        """takes fsa file and label of each dataset in given file"""
        data = SeqIO.read(filepath, "abi")
        D00 = data.annotations["abif_raw"]
        D0 = {k: D00[k] for k in labels}
        return pd.DataFrame(D0)

    D = read_fsa(filepath)

    def get_peak_table(data, threshold):
        """get peak index, start, end and intensity"""
        np.set_printoptions(precision=2)
        variable = savgol_filter(data, 7, 2)
        # min height for track data
        variable[variable <= threshold] = 0
        peaks, coordinates = find_peaks(variable, prominence=1, height=threshold,distance=10)
        P = pd.DataFrame(coordinates)
        P.insert(0, 'peak_index', peaks)
        Dvar = pd.DataFrame(variable, columns=['peak_heights'])
        return P, Dvar

    # peaks for track 9 and 10
    P9 = get_peak_table(D.DATA9, threshold)
    P10 = get_peak_table(D.DATA10, threshold)
    P11 = get_peak_table(D.DATA11, threshold)
    P12 = get_peak_table(D.DATA12, threshold)

    x = np.array(size_values).reshape(-1, 1)
    y = np.array(size_keys).reshape(-1, 1)
    # define x and y parameters
  
    # build a regression model based on sizes and ladder peaks
    ridge = RidgeCV(alphas=np.arange(0.0005, 10, 0.001)).fit(x, y)
    r_alpha = ridge.alpha_
    model_ridge = Ridge(alpha=r_alpha).fit(x, y)

    track_list = []
    track_list.append(P9[0])
    track_list.append(P10[0])
    track_list.append(P11[0])
    track_list.append(P12[0])

    # predict corresponding bases of the track DATA peaks
    all_peak_data = []
    peak_dict = {}
    for i, r in zip(track_list, range(len(track_list))):
        if len(i) > 0:
            # predict corresponding bases of the DATA9 peaks
            values = np.array(i.peak_index).reshape(-1, 1)
            ridge_x_predicted = model_ridge.predict(values)
            # insert predicted bases to the peak table
            i['predicted_base'] = ridge_x_predicted.round(3)
            locals()['D'+str(r+9)
                     ] = i.query('predicted_base>0').reset_index(drop=True)
            locals()['D'+str(r+9)].index += 1
            peak_dict['D'+str(r+9)] = locals()['D'+str(r+9)]
        else:
            locals()['D'+str(r+9)] = i

    all_peak_data.append(peak_dict)
    

    return P9[1], P10[1], P11[1], P12[1], all_peak_data, model_ridge


########################################################################


def get_sizes_v2(filepath, size_keys, size_values, threshold):

    def read_fsa(filepath, labels=['DATA9', 'DATA10', 'DATA11', 'DATA12']):
        """takes fsa file and label of each dataset in given file"""
        data = SeqIO.read(filepath, "abi")
        D00 = data.annotations["abif_raw"]
        D0 = {k: D00[k] for k in labels}
        return pd.DataFrame(D0)

    D = read_fsa(filepath)

    def get_peak_table(data, threshold):
        """get peak index, start, end and intensity"""
        np.set_printoptions(precision=2)
        variable = savgol_filter(data, 7, 2)
        # min height for track data
        #variable[variable<=threshold] = 0
        # height=20 idi.
        peaks, coordinates = find_peaks(
            variable, prominence=1, height=threshold,distance=10)
        P = pd.DataFrame(coordinates)
        P.insert(0, 'peak_index', peaks)
        Dvar = pd.DataFrame(variable, columns=['peak_heights'])
        return P, Dvar

    # peaks for track 9 and 10
    P9 = get_peak_table(D.DATA9, threshold)
    P10 = get_peak_table(D.DATA10, threshold)
    P11 = get_peak_table(D.DATA11, threshold)
    P12 = get_peak_table(D.DATA12, threshold)

    x = np.array(size_values).reshape(-1, 1)
    y = np.array(size_keys).reshape(-1, 1)
    # define x and y parameters

    # build a regression model based on sizes and ladder peaks
    ridge = RidgeCV(alphas=np.arange(0.0005, 10, 0.001)).fit(x, y)
    r_alpha = ridge.alpha_
    model_ridge = Ridge(alpha=r_alpha).fit(x, y)

    track_list = []
    track_list.append(P9[0])
    track_list.append(P10[0])
    track_list.append(P11[0])
    track_list.append(P12[0])

    # predict corresponding bases of the track DATA peaks
    all_peak_data = []
    peak_dict = {}
    for i, r in zip(track_list, range(len(track_list))):
        if len(i) > 0:
            # predict corresponding bases of the DATA9 peaks
            values = np.array(i.peak_index).reshape(-1, 1)
            ridge_x_predicted = model_ridge.predict(values)
            # insert predicted bases to the peak table
            i['predicted_base'] = ridge_x_predicted.round(3)
            locals()['D'+str(r+9)
                     ] = i.query('predicted_base>0').reset_index(drop=True)
            locals()['D'+str(r+9)].index += 1
            peak_dict['D'+str(r+9)] = locals()['D'+str(r+9)]
        else:
            locals()['D'+str(r+9)] = i

    all_peak_data.append(peak_dict)
    

    return P9[1], P10[1], P11[1], P12[1], all_peak_data, model_ridge

########################################################################

def get_threshold(filepath,threshold):
    
    def read_fsa(filepath, labels=['DATA9', 'DATA10', 'DATA11', 'DATA12']):
        """takes fsa file and label of each dataset in given file"""
        data = SeqIO.read(filepath, "abi")
        D00 = data.annotations["abif_raw"]
        D0 = {k: D00[k] for k in labels}
        return pd.DataFrame(D0)

    D = read_fsa(filepath)
    
    
    def get_peak_table(data, threshold):
        """get peak index, start, end and intensity"""
        np.set_printoptions(precision=2)
        variable = savgol_filter(data, 7, 2)
        # min height for track data
        #variable[variable<=threshold] = 0
        # height=20 idi.
        peaks, coordinates = find_peaks(
            variable, prominence=1, height=threshold)
        P = pd.DataFrame(coordinates)
        P.insert(0, 'peak_index', peaks)
        Dvar = pd.DataFrame(variable, columns=['peak_heights'])
        return P, Dvar

    # peaks for track 9 and 10
    P9 = get_peak_table(D.DATA9, threshold)
    P10 = get_peak_table(D.DATA10, threshold)
    P11 = get_peak_table(D.DATA11, threshold)
    P12 = get_peak_table(D.DATA12, threshold)
    
    track_list = []
    track_list.append(P9[0])
    track_list.append(P10[0])
    track_list.append(P11[0])
    track_list.append(P12[0])
    
    track_heights=[] 
    for i in track_list: 
        track_heights.extend(list(i['peak_heights'])) 
        track_heights.sort(reverse=True) 
    threshold_mean=int(np.array(track_heights[:15]).mean().round())
    
    return threshold_mean     
    
########################################################################

server = Flask(__name__)
app = dash.Dash(external_stylesheets=[dbc.themes.BOOTSTRAP,dbc.icons.FONT_AWESOME, dbc.icons.BOOTSTRAP],server=server)

app.layout= html.Div([
    html.A([ 
        dbc.Row([ 
            dbc.Col([ 
                html.A([
                    html.Img(src=app.get_asset_url('helikslogo.jpg'),
                                style={'align':'center',
                                       'width':'30%'
                                       })
                    
                    
                    ],href='https://heliksdna.com/'),

                ]),
        
            dbc.Col([ 
                html.H2(children='GENE MAPPER', style={'text-align': 'center', 'color': 'dark-grey'})
                ]),
                
            dbc.Col()
            
             ], style={
                    'height':'auto',
                    'width':'auto',
                    'display':'flex',
                    'background-color':'white',
                    'align-items': 'center'
                }),
        
    ]),
                 
                 
    html.Div([
        
        dcc.Tabs(id='genemapper-tabs', value='Upload Data', children=[
            dcc.Tab(
                label='Step 1:Upload Data', 
                value='Upload Data',
                children=html.Div([ #start
                    dbc.Row([ 
                        dbc.Col([ 
                            html.Br(),
                            html.Br(),
                            dbc.Row([
                                html.Div('Upload your fsa file'),
                                dcc.Upload(id='fsa',
                                           children= html.Div([
                                               'Drag and Drop or ', 
                                               html.A('Select a File')], 
                                               style={  
                                                   'width': '100%',
                                                   'height': '65px',
                                                   'lineHeight': '60px',
                                                   'borderWidth': '1px',
                                                   'borderStyle': 'dashed',
                                                   'borderRadius': '5px',
                                                   'textAlign': 'center'}),multiple=True)
                                ]),
                            
                            html.Br(),
                            dbc.Row([
                                html.Div('Upload your bin file'),
                                dcc.Upload(id='bin', 
                                           children=html.Div([ 
                                               'Drag and Drop or ',
                                               html.A('Select a File')], 
                                               style={ 
                                                   'width': '100%',
                                                   'height': '65px',
                                                   'lineHeight': '60px',
                                                   'borderWidth': '1px',
                                                   'borderStyle': 'dashed',
                                                   'borderRadius': '5px',
                                                   'textAlign': 'center'}),multiple=True)
                                ]),
                            
                            html.Br(),
                            dbc.Row([ 
                                html.Div('Upload or Select your Liz Standard'),
                                html.Div([
                                    html.Div([
                                        dcc.Upload(id='liz', 
                                                   children=html.Div([ 
                                                       'Drag and Drop or ', 
                                                       html.A('Select a File')],  
                                                       style={ 
                                                           'width': '100%',
                                                           'height': '65px',
                                                           'lineHeight': '60px',
                                                           'borderWidth': '1px',
                                                           'borderStyle': 'dashed',
                                                           'borderRadius': '5px',
                                                           'textAlign': 'center'}),multiple=True)
                                        ]),
                                    
                                    html.Div([
                                        dcc.Dropdown(['HeliksSizeS','HeliksSizeS2','HeliksSizeL',
                                                      'Liz120','Liz500'], id='select-sizes', style=dict(width='100%'))
                                    ])


                                     
                                    ],style={'display':'inline-block'})

                             
                                ])

                                ],style={'margin-right':'30px'}
                                    
                                    ),
                        
                        dbc.Col([
                            html.Br(),
                            html.Br(),
                            dbc.Row([
                                html.Div('Enter the panel name'),
                                html.Div([
                                    html.Div([
                                        dcc.Input(id='panel-name', type='string',placeholder="Type your panel name here...",debounce=False,step=1,style=dict(heigh='27px',width='290px'))
                                        
                                        ]),
                                    html.Div([ 
                                        dcc.Dropdown(['FMF_HELIKS_V2','FMF-V2_HLKS','FMF'], id='select-panel', style=dict(width='290px'))  ])
                                    
                                    
                                    ])
            
                                ]),
                            html.Br(),
                            dbc.Row([
                                html.Div('Enter the project name'),
                                html.Div([ 
                                    html.Div([
                                        dcc.Input(id='project-name', type='string',placeholder="Type your project name here...",debounce=False,step=1,style=dict(heigh='27px',width='290px'))
                                        
                                        ])
                                    ])
                                
                                ]),

                            
                            html.Br(),    
                            dbc.Row([
                                html.Div('Give threshold value for Size Ladder'),
                                html.Div([
                                    html.Div([
                                        dcc.Input(id='sizeladder-threshold',type='number',value=250,placeholder="Type detected threshold value here...",debounce=False, step=1,style=dict(height='27px',width='290px'))
                                                
                                        ])
                                    
                                    ])
                                
                                ]),
                            html.Br(),
                            
                            
                            dbc.Row([
                                html.Div('Upload project file'),
                                dcc.Upload(id='project', 
                                           children=html.Div([ 
                                               'Drag and Drop or ', 
                                               html.A('Select a File')],  
                                               style={ 
                                                   'width': '100%',
                                                   'height': '60px',
                                                   'lineHeight': '60px',
                                                   'borderWidth': '1px',
                                                   'borderStyle': 'dashed',
                                                   'borderRadius': '5px',
                                                   'textAlign': 'center'}))
                                ])
                                
                                ])
                        
                    
                        ]),
                    html.Br(),
                    dbc.Row(
                        html.Div([
                            dbc.Button(
                                id='first-submit',
                                            n_clicks=0,
                                            children='Submit',
                                            style={'fontSize': 21,
                                                   'backgroundColor': '#7FB6C6',
                                                   'color': 'white',
                                                   'height': '50px',
                                                   'text-align': 'center'},size="lg", className="me-1"
                                
                                ),
                            
                            dbc.Button(id="return-size-matching", color="light",
                                       className="bi bi-arrow-right-square",
                                       style={
                                           'fontSize':21,
                                           'height':'50px',
                                           'text-align':'center'
                                           }
                                       )
                            
                            
                            ],style={'textAlign':'center'})

                        ),
                    html.Br(),
                    dbc.Row([ 
                        html.H2("Action List"),
                        html.Ul(id="action-list",children='no action available yet!')
                        
                        ])#finish
                    
                
                    
                    ])
                
                
                
                ),
            
            dcc.Tab(label='Step 2:Size Matching',
                    value='Size Matching',
                    children=html.Div([
                    html.Div( 
                    id='edit-size'
                    ),
                    dcc.Graph(id='graph'),
                    html.Div([
                        dbc.Button(id='submit-button',
                                   n_clicks=0,
                                   children='Submit',
                                   style={'fontSize': 21,
                                          'backgroundColor': '#7FB6C6',
                                          'color': 'white',
                                          'height': '50px',
                                          'text-align': 'center'}) 
                        ],style={'display': 'flex','textAlign':'center', 'width': '250px','margin':'auto'})
               
                                ])), #finish-2 
            
            dcc.Tab( 
                label='Step 3:Size Graph', 
                id='Size Graph',
                value='Size Graph',
                children=html.Div([
                    dcc.Tabs(id="subtabs1", value='Subtab1', children=[

                        dcc.Tab(label='Size Estimating', id='subtab1', value='Subtab1',
                                children=html.Div([
                                    html.Div([
                                        html.Div('Set threshold:',
                                                  style=dict(padding='10px')),
                                        dcc.Input(id='threshold', type='number', debounce=False, step=1, style=dict(
                                            height='28px')),
                                        html.Abbr("\u2754", title="This section is for setting the threshold value of the minimum peak heights to be detected. ", style=dict(
                                            color='white'))

                                    ], style={'display': 'flex', 'margin-left': '10px'}),
                                    html.Br(),

                                    html.Div(children=[
                                        html.Div('Attitude for peaks:',
                                                 style=dict(padding='10px')),
                                        dcc.Dropdown(
                                            ['Remove', 'Ignore'], 'Remove', id='select-mode', style=dict(width='145px')),
                                        html.Abbr("\u2754", title="This section is about whether you want to delete or just ignore the peaks with peak height below the threshold value.", style=dict(
                                            color='white'))
                                    ], style={'display': 'flex', 'width': '70%', 'margin-left': '10px'}),
                                    dcc.Graph(id='size-graph')
                                ])

                                ),



                        dcc.Tab(label='Show Bins', id='subtab2', value='Subtab2',
                                children=html.Div([
                                    dcc.Graph(id='bin-graph'),
                                    html.Div(id='bin-tables'),

                                    html.Div([
                                        dbc.Button(id='call',
                                                    n_clicks=0,
                                                    children='Call Alleles',
                                                    style={'fontSize': 21,
                                                           'backgroundColor': '#7FB6C6',
                                                           'color': 'white',
                                                           'height': '50px',
                                                           'text-align': 'center'}
                                                    ),
                                        dbc.Button(id='save',
                                                    n_clicks=0,
                                                    children='Save State',
                                                    style={'fontSize': 21,
                                                           'backgroundColor': '#7FB6C6',
                                                           'color': 'white',
                                                           'height': '50px',
                                                           'text-align': 'center'})
                                        ], style={'display': 'flex','textAlign':'center', 'width': '250px','margin':'auto'}),
                                    dbc.Modal(
                                        [
                                            dbc.ModalHeader(dbc.ModalTitle("Save Project State")),
                                            dbc.ModalBody("Project state is saved! You can find the project state file under the project folder."),
                                            dbc.ModalFooter( dbc.Button("Close", id="close", n_clicks=0))
                                        
                                        ],
                                        id="modal-save",
                                        centered=True,
                                        is_open=False,
                                    )
    
    
                                ]))

                    ])

                ])
            ), #finish-step3
            
            dcc.Tab(
                label='Step 4:Allele Calling',
                id='allele-call',
                value='Allele Calling',
                children=html.Div([
                    html.Div('Captured alleles will be seen here...',id='alleles-header',style={'padding':'5px'}),
                    html.Div(
                        id='alleles'                       
                        )
           
                    ]

                )
            ) #finish step4
            
            
            
            
            
            
            
            
            ]) #tabs finish
        
        
    
        
        
        ])
    

    
    
    ])

                 
        
##-----------------------------------------------------------------------##                 
@app.callback(
    Output("action-list", "children"),
    State("project-name", "value"),
    State("sizeladder-threshold","value"),
    [State('fsa','filename'),
    State('fsa','contents') ],
    State('select-sizes','value'),
    [State('bin','filename'),
     State('bin','contents')],
    [State('liz','filename'),
     State('liz','contents')],
    State('panel-name','value'),
    State('select-panel','value'),
    Input('first-submit','n_clicks'),
    prevent_initial_call=True
)

def update_output(projectname, sizethreshold, uploaded_filenames, uploaded_file_contents,selectedsize,bin_filename,bin_filecontent,liz_filename,liz_filecontent,panel_name,panelselect,n_clicks): 
    
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project") 
    UPLOAD_DIRECTORY =os.path.join(UPLOAD_DIRECTORY,str(projectname))
    UPLOAD_DIRECTORY = os.path.join(UPLOAD_DIRECTORY ,"app_uploaded_files")
    if not os.path.exists(UPLOAD_DIRECTORY):
        os.makedirs(UPLOAD_DIRECTORY) 
        
    #######panel kısmı
    def return_panel():
        if panel_name is not None:
            panel=panel_name
        if panel_name is None and panelselect is not None:
            panel=panelselect
        return panel
    
    panel_name=return_panel()
        
    print('PANEL NAMEEEEE',panel_name)
        
        
    
    def save_files(uploaded_filename,uploaded_filecontents,UPLOAD_DIRECTORY): 
        if uploaded_filename is not None and uploaded_filecontents is not None:
            for name,content in zip(uploaded_filename,uploaded_filecontents):
                data = content.encode("utf8").split(b";base64,")[1]
                with open(os.path.join(UPLOAD_DIRECTORY, name), "wb") as fp:
                    fp.write(base64.decodebytes(data))
    
    save_files(uploaded_filenames,uploaded_file_contents,UPLOAD_DIRECTORY)
    save_files(bin_filename,bin_filecontent,UPLOAD_DIRECTORY) 
    save_files(liz_filename,liz_filecontent,UPLOAD_DIRECTORY)



    input_dict={}
   
    
   
    ################################################################################# 

    def uploaded_files():
        """List the files in the upload directory."""
        files = [] 
        for filename in os.listdir(UPLOAD_DIRECTORY): 
            path = os.path.join(UPLOAD_DIRECTORY, filename)
            print(path)
            if os.path.isfile(path): 
                files.append(filename) 
        return files
    
    #################################################################################
    binpath=os.path.join(UPLOAD_DIRECTORY,bin_filename[0])
    write_bin_list(binpath,panel_name)
    #################################################################################    
    
    
    files = uploaded_files()

    if projectname is not None and len(projectname) !=0:
        input_dict['projectname']='Project name is entered: {} '.format(projectname)
    
    if sizethreshold is not None:
        input_dict['sizethreshold']='Threshold value is set: {}'.format(sizethreshold)
    
    
    if uploaded_filenames is not None and len(uploaded_filenames) !=0: 
        #print('uploaded file names',uploaded_filenames)
        input_dict['fsa']='fsa file is uploaded: {}'.format(uploaded_filenames[0])
        
        
        
    if bin_filename is not None and len(bin_filename) !=0:
        input_dict['bin']='bin file is uploaded: {}'.format(bin_filename[0])
        
        
    if liz_filename is not None and len(liz_filename) !=0:
        input_dict['liz']='liz standard is uploaded: {}'.format(liz_filename[0])
    
    
    if selectedsize is not None:
        input_dict['liz']='liz standard is selected: {}'.format(selectedsize)
        
        
    if panel_name is not None and len(panel_name) !=0:
        input_dict['panel_name']='Panel name is entered: {} '.format(panel_name)
        
    if panelselect is not None:
        input_dict['panel_name']='Panel name is selected: {} '.format(panelselect)
  
           
    if len(input_dict)!= 0 : 
        return [html.Li(i) for i in input_dict.values() ]
    else:
        return 'no action available yet!'
    
    
##---------------------------------------------------------------------------##    

@app.callback(Output('edit-size','children'),
              State('fsa','filename'),
              State('liz','filename'),
              State('select-sizes','value'),
              State('project-name','value'),
              Input('return-size-matching','n_clicks'),
              State('sizeladder-threshold','value'),
              prevent_initial_call=True
    )
def size_screen(fsa,liz,selectedsize,projectname,n_clicks,threshold):
    
    #print('selected siiiiiizeeee')
    #print(selectedsize)
    #print(type(selectedsize))
    
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project") 
    UPLOAD_DIRECTORY =os.path.join(UPLOAD_DIRECTORY,str(projectname))
    UPLOAD_DIRECTORY = os.path.join(UPLOAD_DIRECTORY ,"app_uploaded_files")
    if not os.path.exists(UPLOAD_DIRECTORY):
        os.makedirs(UPLOAD_DIRECTORY) 
        
    
    path=os.path.join(UPLOAD_DIRECTORY, fsa[0])
     
    #liz part-eğer liz seçilmişse onu kullan, yüklendiyse dosyadan oku
    
    def return_sizes():
        if liz is not None: 
            liz_path=os.path.join(UPLOAD_DIRECTORY,liz[0]) 
            sizes=read_sizes(liz_path)
            
        if liz is None and selectedsize is not None:
            sizes = size_dict[selectedsize]
            #print('sizeeees0',sizes)
        return sizes

    sizes=return_sizes()
    liz = LizStd(path.encode('unicode_escape').decode(),sizes,threshold)
    #print('lizzzz')
    size_table=dash_table.DataTable( 
        id='edit-size-table',
        columns=[{'id': 'status', 'name': 'status'}],
        data=[
            {'status': i} for i in list(liz[0].values())
            ],
        editable=True
        )
    
    #print('size tableee',size_table)
    
    
    return size_table
    

##---------------------------------------------------------------------------##  
 
@app.callback(Output('graph','figure'),
              Input('edit-size-table', 'data'),
              State('fsa','filename'),
              State('liz','filename'),
              State('select-sizes','value'),
              State('project-name','value'),
              State('sizeladder-threshold','value'),
              prevent_initial_call=True
    )
def size_graph(data,fsa,liz,selectedsize,projectname,threshold):
    #grafiği ayrı return etmemizin amacı, size_screen'de return edilen datatable'ı burda input olarak vermek
    #ve her yapılan değişiklikte call back'i aktive edip grafikte yansıtmak.
    
    #print('dataaa----------aa',data)
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project") 
    UPLOAD_DIRECTORY =os.path.join(UPLOAD_DIRECTORY,str(projectname))
    UPLOAD_DIRECTORY = os.path.join(UPLOAD_DIRECTORY ,"app_uploaded_files")
    if not os.path.exists(UPLOAD_DIRECTORY):
        os.makedirs(UPLOAD_DIRECTORY) 
        
        
    path=os.path.join(UPLOAD_DIRECTORY, fsa[0])
    
    
    #liz part-eğer liz seçilmişse onu kullan, yüklendiyse dosyadan oku
    def return_sizes():
        if liz is not None: 
            liz_path=os.path.join(UPLOAD_DIRECTORY,liz[0]) 
            sizes=read_sizes(liz_path)
            
        if liz is None and selectedsize is not None:
            sizes = size_dict[selectedsize]
    
        return sizes
    
    sizes=return_sizes() 
    
    liz = LizStd(path.encode('unicode_escape').decode(),sizes,threshold)
    
    #graph part
    df = pd.DataFrame(data, columns=['status'])
    
    #print('df statuuu-----uuus',df)
    #print('dataaa666aa',data)
    
    DetectedPeak = go.Scatter(x=liz[1],
                              y=liz[2][liz[1]],
                              mode='markers+text',
                              text=df['status'],
                              textposition='top center',
                              name='peaks',
                              marker=dict(
        size=10,
        color='black')
    )
    SizeLadder = go.Scatter(y=liz[2],
                            mode='lines',
                            name='size ladder',
                            line=dict(color='orange'))

    data = [DetectedPeak, SizeLadder]

    layout = go.Layout(
        title='Size Matching (tol={})'.format(liz[3]),
        hovermode='closest'
    )
    
    #print('text  df['status]')
    #print(df['status'])

    fig = go.Figure(data=data, layout=layout, layout_xaxis_range=[
                    liz[1].min()-200, liz[1].max()+200])
    fig.update_traces(
        hovertemplate="<b><i>height</i>: %{y}</b> <br><b><i>datapoint</i>: %{x}</b>")
    
    return fig

##---------------------------------------------------------------------------##  

@app.callback(Output('size-graph', 'figure'),
              Output('threshold','value'),
              Input('threshold', 'value'),
              State('edit-size-table', 'data'),
              State('fsa','filename'),
              State('liz','filename'),
              State('select-sizes','value'),
              State('sizeladder-threshold','value'),
              State('project-name','value'),
              Input('select-mode', 'value'),
              [Input('submit-button', 'n_clicks')])

def SizeEstimatingTab(threshold,data,fsa,liz,selectedsize,ladderthreshold,projectname,drop,n_clicks):
    
    
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project") 
    UPLOAD_DIRECTORY =os.path.join(UPLOAD_DIRECTORY,str(projectname))
    UPLOAD_DIRECTORY = os.path.join(UPLOAD_DIRECTORY ,"app_uploaded_files")
    if not os.path.exists(UPLOAD_DIRECTORY):
        os.makedirs(UPLOAD_DIRECTORY) 
        
    print('THRESHOLDDDDDDDD')
    print(threshold)
    path=os.path.join(UPLOAD_DIRECTORY, fsa[0])
    
    if threshold is None: 
        avg_threshold=int((get_threshold(path, 150)*30)/100)
    else:
        avg_threshold=threshold
    
    
    ####-------------------------sizeestimatingpart------------------------####
    
    #liz kısmı
    def return_sizes():
        if liz is not None: 
            liz_path=os.path.join(UPLOAD_DIRECTORY,liz[0]) 
            sizes=read_sizes(liz_path)
            
        if liz is None and selectedsize is not None:
            sizes = size_dict[selectedsize]
    
        return sizes
    
    sizes=return_sizes() 
    
    liz = LizStd(path.encode('unicode_escape').decode(),sizes,ladderthreshold)
    

    
    
    
    
    df = pd.DataFrame(data, columns=['status'])
    
    
    size_std = list(df['status'])
    liz_peaks = list(liz[1])


    sizepair = {}
    for i, j in zip(size_std, liz_peaks):
        if i == 'no match' or i=="":
            pass
        else:
            sizepair[int(i)] = j
            
    #print('sizepaiiiiiiiiiiir',sizepair)
            
    get = get_sizes_v2(path, list(sizepair.keys()),
                       list(sizepair.values()), avg_threshold)
    prediction_model = get[5]
    plot_index = {'D9': get[0], 'D10': get[1], 'D11': get[2], 'D12': get[3]}
    

    if drop == 'Remove' and n_clicks != 0:
        get = get_sizes(path, list(sizepair.keys()),
                        list(sizepair.values()), avg_threshold)
        prediction_model = get[5]

        plot_index = {'D9': get[0], 'D10': get[1],
                      'D11': get[2], 'D12': get[3]}
        fig = make_subplots(rows=len(get[4][0].keys()), cols=1)
        
        
        print('getttttttttttttt',get[4][0]['D9'].columns)
        print('plotttt index',plot_index)

        for j, i in enumerate(list(get[4][0].keys())):
            # convert index to size
            plot_index[i]['index'] = plot_index[i].index
            track_predict = prediction_model.predict(
                np.array(plot_index[i]['index']).reshape(-1, 1))
            plot_index[i]['track_size'] = track_predict
            

            # plot only one track and put bases as labels
            plot1 = go.Scatter(x=get[4][0][i]['predicted_base'],
                               y=get[4][0][i]['peak_heights'],
                               mode='markers+text',
                               text=get[4][0][i]['predicted_base'],
                               textposition='top center',
                               name=dye_dict[str(j+1)]+'-peaks',
                               marker=dict(
                                   color=dye_dict[str(j+1)], size=8)
                               )
            plot2 = go.Scatter(x=plot_index[i][plot_index[i]['track_size'] <= int(list(get[4][0][i]['predicted_base'][-1:])[0])+50].query('track_size >0')['track_size'],
                               y=plot_index[i][plot_index[i]['track_size'] <= int(list(
                                   get[4][0][i]['predicted_base'][-1:])[0])+50].query('track_size >0')['peak_heights'],
                               mode='lines', name=dye_dict[str(j+1)]+'-lines',
                               line=dict(color=dye_dict[str(j+1)]))
            data = [plot1, plot2]
            for a in data:
                fig.append_trace(a, row=j+1, col=1)
                fig.update_traces(
                    hovertemplate="<b><i>height</i>: %{y}</b> <br><b><i>size</i>: %{x}</b>")
            fig.update_xaxes(showspikes=True, spikemode="across")
            fig.update_layout(height=1500, title="Patient Data Plots")
            
    if drop == 'Ignore' and n_clicks != 0:

        fig = make_subplots(rows=len(get[4][0].keys()), cols=1)

        for j, i in enumerate(list(get[4][0].keys())):
            # convert index to size
            plot_index[i]['index'] = plot_index[i].index
            track_predict = prediction_model.predict(
                np.array(plot_index[i]['index']).reshape(-1, 1))
            plot_index[i]['track_size'] = track_predict

            # plot only one track and put bases as labels
            plot1 = go.Scatter(x=get[4][0][i]['predicted_base'],
                               y=get[4][0][i]['peak_heights'],
                               mode='markers+text',
                               text=get[4][0][i]['predicted_base'],
                               textposition='top center',

                               name=dye_dict[str(j+1)]+'-peaks',
                               marker=dict(
                                   color=dye_dict[str(j+1)], size=8)
                               )
            plot2 = go.Scatter(x=plot_index[i][plot_index[i]['track_size'] <= int(list(get[4][0][i]['predicted_base'][-1:])[0])+50].query('track_size >0')['track_size'],
                               y=plot_index[i][plot_index[i]['track_size'] <= int(list(
                                   get[4][0][i]['predicted_base'][-1:])[0])+50].query('track_size >0')['peak_heights'],
                               mode='lines', name=dye_dict[str(j+1)]+'-lines',
                               line=dict(color=dye_dict[str(j+1)]))
            data = [plot1, plot2]
            for a in data:
                fig.append_trace(a, row=j+1, col=1)
                fig.update_traces(
                    hovertemplate="<b><i>height</i>: %{y}</b> <br><b><i>size</i>: %{x}</b>")

            fig.update_xaxes(showspikes=True, spikemode="across")
            fig.update_layout(height=1500, title="Patient Data Plots")

            
    
    return fig,avg_threshold

##-----------------------------------------------------------------------##  


@app.callback(Output('bin-graph', 'figure'),
              Input('threshold', 'value'),
              State('edit-size-table', 'data'),
              State('bin','filename'),
              State('fsa','filename'),
              State('liz','filename'),
              State('select-sizes','value'),
              State('sizeladder-threshold','value'),
              State('project-name','value'),
              State('panel-name','value'),
              State('select-panel','value'),
              Input('select-mode', 'value'),
              [Input('submit-button', 'n_clicks')],
              [Input('bin-table{}'.format(i), 'derived_virtual_data')
               for i in range(read_bin_list())],
              [Input('bin-table{}'.format(i), 'derived_virtual_selected_rows')
               for i in range(read_bin_list())])

def ShowBinsTab(threshold,data,binfile,fsa,liz,selectedsize,ladderthreshold,projectname,panel_name,panelselect,drop,n_clicks,*bin_data): 
   
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project") 
    UPLOAD_DIRECTORY =os.path.join(UPLOAD_DIRECTORY,str(projectname))
    UPLOAD_DIRECTORY = os.path.join(UPLOAD_DIRECTORY ,"app_uploaded_files")
    if not os.path.exists(UPLOAD_DIRECTORY):
        os.makedirs(UPLOAD_DIRECTORY) 
        
    print('bin dataaaaaaaaaaaaaaaaa',bin_data)   
    l=int(len(bin_data)/2)
    marked_rows=bin_data[l::1]
    rows=bin_data[:l:1]
    
    selected_rows=[]
    for c,i in enumerate(marked_rows):
        try:
            for j in i:
                #print('j',j,'C',c)
                #print('related rows',rows[c][j])
                selected_rows.append(rows[c][j])
        except:
            selected_rows=[]
    
    print('SELECTED ROWS')
    
    print(selected_rows)
            
        
    #print('THRESHOLDDDDDDDD')
    #print(threshold)
    path=os.path.join(UPLOAD_DIRECTORY, fsa[0])
    binpath=os.path.join(UPLOAD_DIRECTORY,binfile[0])
    
    if threshold is None: 
        avg_threshold=int((get_threshold(path, 150)*30)/100)
    else:
        avg_threshold=threshold
    
        
    ####-------------------------sizeestimatingpart------------------------####
    
    #liz kısmı
    def return_sizes():
        if liz is not None: 
            liz_path=os.path.join(UPLOAD_DIRECTORY,liz[0]) 
            sizes=read_sizes(liz_path)
            
        if liz is None and selectedsize is not None:
            sizes = size_dict[selectedsize]
    
        return sizes
    
    sizes=return_sizes() 
    
    
    liz = LizStd(path.encode('unicode_escape').decode(),sizes,ladderthreshold)
    
     
    #######panel kısmı
    def return_panel():
        if panel_name is not None:
            panel=panel_name
        if panel_name is None and panelselect is not None:
            panel=panelselect
        return panel
    
    panel_name=return_panel()
        
    print('PANEL NAMEEEEE',panel_name)
        
    
    df = pd.DataFrame(data, columns=['status'])
    
    
    size_std = list(df['status'])
    liz_peaks = list(liz[1])


    sizepair = {}
    for i, j in zip(size_std, liz_peaks):
        if i == 'no match' or i=="":
            pass
        else:
            sizepair[int(i)] = j
    
    bins = read_bins(binpath,panel_name)
    
    ############get sizes for ignore option#############
    get = get_sizes_v2(path, list(sizepair.keys()),
                       list(sizepair.values()), avg_threshold)
    prediction_model = get[5]
    plot_index = {'D9': get[0], 'D10': get[1], 'D11': get[2], 'D12': get[3]}
    #####################################################
    
    if drop == 'Remove' and n_clicks != 0: 
        print('removeeeee')
        get = get_sizes(path, list(sizepair.keys()),
                        list(sizepair.values()), threshold)
        prediction_model = get[5]

        plot_index = {'D9': get[0], 'D10': get[1],
                      'D11': get[2], 'D12': get[3]}
        fig_binr = go.Figure()
        pick = []
        for j, i in enumerate(list(get[4][0].keys())):
            # convert index to size
            plot_index[i]['index'] = plot_index[i].index
            track_predict = prediction_model.predict(
                np.array(plot_index[i]['index']).reshape(-1, 1))
            plot_index[i]['track_size'] = track_predict

            # plot only one track and put bases as labels
            plot = go.Scatter(x=plot_index[i][plot_index[i]['track_size'] <= int(list(get[4][0][i]['predicted_base'][-1:])[0])+50].query('track_size >0')['track_size'],
                              y=plot_index[i][plot_index[i]['track_size'] <= int(list(
                              get[4][0][i]['predicted_base'][-1:])[0])+50].query('track_size >0')['peak_heights'],
                              mode='lines', name=dye_dict[str(j+1)],
                              line=dict(color=dye_dict[str(j+1)]))

            pick.append(plot)

        # panel kısmı

        colors = []
        for i in range(len(bins[1])):
            bins[1][i]['color'] = bins[1][i]['color'].str.lower()
            colors.append(bins[1][i].color[0])
            colors.append(bins[1][i].color[1].lower())

        colors = list(set(colors))
        
     

        color_dict = {}
        for color in colors: 
            color_dict[color] = []

        for color in colors:
             for i in range(len(bins[1])):
                color_dict[color].append(
                    list(bins[1][i]['allele'][bins[1][i]['color'] == color]))

        for i in color_dict.keys():
             color_dict[i] = list(
                set([a for sublists in color_dict[i] for a in sublists]))

        # print('color_dict:')
        # print(color_dict)

        #######

        fig_binr.add_traces(pick)
    
    
    
        bin_list = []
        for i in range(len(bins[1])):
            locals()['df-{}'.format(i)] = bins[1][i]
            bin_list.append(locals()['df-{}'.format(i)])
    
    
        #add bins
        if len(selected_rows) != 0: 
            for i in selected_rows:
                    fig_binr.add_shape(x0=i['start_size'], x1=i['end_size'], y0=20000, xref='x1',
                                   fillcolor=i['color'], opacity=0.2, name=i['allele'], editable=True)

            fig_binr.update_traces(
                hovertemplate="<b><i>height</i>: %{y}</b> <br><b><i>size</i>: %{x}</b>")
            fig_binr.update_xaxes(showspikes=True, spikemode="across")
            fig_binr.update_layout(height=500, title="Patient Data Plots")
        
    
    
    if drop == 'Ignore' and n_clicks != 0:
        print('ignoreeee')
        fig_binr = go.Figure()
        pick = []
        for j, i in enumerate(list(get[4][0].keys())):
            # convert index to size
            plot_index[i]['index'] = plot_index[i].index
            track_predict = prediction_model.predict(
                np.array(plot_index[i]['index']).reshape(-1, 1))
            plot_index[i]['track_size'] = track_predict

            # plot only one track and put bases as labels
            plot = go.Scatter(x=plot_index[i][plot_index[i]['track_size'] <= int(list(get[4][0][i]['predicted_base'][-1:])[0])+50].query('track_size >0')['track_size'],
                              y=plot_index[i][plot_index[i]['track_size'] <= int(list(
                                  get[4][0][i]['predicted_base'][-1:])[0])+50].query('track_size >0')['peak_heights'],
                              mode='lines', name=dye_dict[str(j+1)],
                              line=dict(color=dye_dict[str(j+1)]))

            pick.append(plot)
            
      

        # panel kısmı
        colors = []
        for i in range(len(bins[1])):
            bins[1][i]['color'] = bins[1][i]['color'].str.lower()
            colors.append(bins[1][i].color[0])
            colors.append(bins[1][i].color[1].lower())

        colors = list(set(colors))

        color_dict = {}
        for color in colors:
            color_dict[color] = []

        for color in colors:
            for i in range(len(bins[1])):
                color_dict[color].append(
                    list(bins[1][i]['allele'][bins[1][i]['color'] == color]))

        for i in color_dict.keys():
            color_dict[i] = list(
                set([a for sublists in color_dict[i] for a in sublists]))

        # print('color_dict:')
        # print(color_dict)

        #######

        fig_binr.add_traces(pick)

        if len(selected_rows) != 0: 
            for i in selected_rows:
                fig_binr.add_shape(x0=i['start_size'], x1=i['end_size'], y0=20000, xref='x1',
                                   fillcolor=i['color'], opacity=0.2, name=i['allele'], editable=True)

        fig_binr.update_traces(
            hovertemplate="<b><i>height</i>: %{y}</b> <br><b><i>size</i>: %{x}</b>")
        fig_binr.update_xaxes(showspikes=True, spikemode="across")
        fig_binr.update_layout(height=500, title="Patient Data Plots")
        
    return fig_binr




@app.callback(Output('bin-tables','children'),
              Input('threshold', 'value'),
              State('edit-size-table', 'data'),
              State('bin','filename'),
              State('fsa','filename'),
              State('liz','filename'),
              State('select-sizes','value'),
              State('sizeladder-threshold','value'),
              State('project-name','value'),
              State('panel-name','value'),
              State('select-panel','value'),
              Input('select-mode', 'value'),
              [Input('submit-button', 'n_clicks')])

def AddBins(threshold,data,binfile,fsa,liz,selectedsize,ladderthreshold,projectname,panel_name,panelselect,drop,n_clicks): 
   
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project") 
    UPLOAD_DIRECTORY =os.path.join(UPLOAD_DIRECTORY,str(projectname))
    UPLOAD_DIRECTORY = os.path.join(UPLOAD_DIRECTORY ,"app_uploaded_files")
    if not os.path.exists(UPLOAD_DIRECTORY):
        os.makedirs(UPLOAD_DIRECTORY) 
        
         
    #######panel kısmı
    def return_panel():
        if panel_name is not None:
            panel=panel_name
        if panel_name is None and panelselect is not None:
            panel=panelselect
        return panel
    
    panel_name=return_panel()
        
    print('PANEL NAMEEEEEee',panel_name)
        
    print('THRESHOLDDDDDDDD')
    #print(threshold)
    path=os.path.join(UPLOAD_DIRECTORY, fsa[0])
    binpath=os.path.join(UPLOAD_DIRECTORY,binfile[0])
    
    if threshold is None: 
        avg_threshold=int((get_threshold(path, 150)*30)/100)
    else:
        avg_threshold=threshold
    

    
    bins = read_bins(binpath,panel_name)
   
    
    bin_list = []
    for i in range(len(bins[1])):
        locals()['df-{}'.format(i)] = bins[1][i]
        bin_list.append(locals()['df-{}'.format(i)])
    
    
    bintable=[
    dash_table.DataTable(id='bin-table{}'.format(j),
                         data=i.to_dict('records'),
                         columns=[{'id': 'allele', 'name': [str(bins[0][j]), 'allele'],'editable':False},
                                  {'id': 'start_size', 'name': [
                                      str(bins[0][j]), 'start_size']},
                                  {'id': 'end_size', 'name': [
                                      str(bins[0][j]), 'end_size']},
                                  {'id': 'color', 'name': [str(bins[0][j]), 'color']}],
                         editable=True,
                         merge_duplicate_headers=True,
                         row_selectable="multi",
                         selected_rows=[],
                         page_action="native",
                         style_cell={
                             'maxWidth': '5px',
                             'whiteSpace': 'normal'},
                         style_table={
                             'overflowX': 'scroll'},
                         style_header={
                             'backgroundColor': 'rgb(49, 51, 53)',
                             'color': 'white',
                             'fontWeight': 'bold'},
                         style_data_conditional=[{'if': {'filter_query':'{color}="Blue"',
                                                         'column_id': col},
                                                  'border': '1px solid blue',
                                                  'color': 'blue'} for col in ['allele','start_size','end_size','color']]+
                                                 [{'if': {'filter_query':'{color}="Green"',
                                                                                 'column_id': col},
                                                   'border': '1px solid green',
                                                   'color': 'green'} for col in ['allele','start_size','end_size','color']]+

                                                  [{'if': {'filter_query':'{color}="Red"',
                                                           'column_id': col},
                                                    'border': '1px solid red',
                                                    'color': 'red'} for col in ['allele','start_size','end_size','color']]+
                                                   [{'if': {'filter_query':'{color}="Purple"',
                                                            'column_id': col},
                                                     'border': '1px solid purple',
                                                     'color': 'purple'} for col in ['allele','start_size','end_size','color']]+
                                                   [{'if': {'filter_query':'{color}="Yellow"',
                                                            'column_id': col},
                                                     'border': '1px solid yellow'} for col in ['allele','start_size','end_size','color']]+
                                                   [{'if': {'column_editable': False},
                                                     'cursor': 'not-allowed'}]                                                                                       
                         ) for j, i in enumerate(bin_list)]
    
    return bintable

##---------------------------save project modal call------------------------------##  


@app.callback(  
    Output("modal-save", "is_open"),
    Input("save", "n_clicks"),
    Input('close','n_clicks'),
    State('project-name','value'),
    State('panel-name','value'),
    State('select-panel','value'),
    State('fsa','filename'),
    State('bin','filename'),
    State('liz','filename'),
    State('select-sizes','value'),
    State('sizeladder-threshold','value'),
    State('edit-size-table', 'data'),
    State("modal-save", "is_open"),
    [State('bin-table{}'.format(i), 'derived_virtual_data')
     for i in range(read_bin_list())])
def save_project(click1,click2,projectname,panel_name,panelselect,fsa,binfile,liz,selectedsize,threshold,data,is_open,*bin_data):
    
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project") 
    UPLOAD_DIRECTORY_ =os.path.join(UPLOAD_DIRECTORY,str(projectname))
    UPLOAD_DIRECTORY = os.path.join(UPLOAD_DIRECTORY_ ,"app_uploaded_files")
    path=os.path.join(UPLOAD_DIRECTORY, fsa[0])
    
    #######panel kısmı
    def return_panel():
        if panel_name is not None:
            panel=panel_name
        if panel_name is None and panelselect is not None:
            panel=panelselect
        return panel
    
    panel_name=return_panel()
        
    print('PANEL NAMEEEEE',panel_name)
    
    binpath=os.path.join(UPLOAD_DIRECTORY,binfile[0])
    bins = read_bins(binpath,panel_name)
    
         
    
    #liz part-eğer liz seçilmişse onu kullan, yüklendiyse dosyadan oku
    
    def return_sizes():
        if liz is not None: 
            liz_path=os.path.join(UPLOAD_DIRECTORY,liz[0]) 
            sizes=read_sizes(liz_path)
            
        if liz is None and selectedsize is not None:
            sizes = size_dict[selectedsize]
            #print('sizeeees0',sizes)
        return sizes

    sizes=return_sizes()
    liz = LizStd(path.encode('unicode_escape').decode(),sizes,threshold)

    #currentdir = os.getcwd()
    #final_dir=os.path.join(currentdir, 'project')   
    #final_dir=os.path.join(final_dir,projectname)
    size_ladder=[]
    if click1 or click2:
        if not os.path.exists(UPLOAD_DIRECTORY_):
            os.makedirs(UPLOAD_DIRECTORY_, exist_ok=True) 
            #print(UPLOAD_DIRECTORY_)
            #print('done')
        file=os.path.join(UPLOAD_DIRECTORY_,'{}-state-{}.txt'.format(projectname,click1))
        f = open(file,"w")
        f.write('##MAIN FILE PATH :\n')
        f.write(path)
        f.write('\n\n##PANEL NAME :\n')
        f.write(panel_name) 
        #print('dataaaaaa',data)
        #print(list(liz[1]))
        for i in data:
            #print('iiiiii')
            #print(i['status'])
            size_ladder.append(i['status'])
            
        f.write('\n\n##SIZE LADDER :\n')
        f.write(str(size_ladder))
        f.write('\n\n##SIZE PAIR :\n')
        f.write(str(list(liz[1])))
        f.write('\n\n##MARKER LIST :\n')
        f.write(str(bins[0]))
        f.write('\n\n##BIN DATA :')
        for j,i in enumerate(bin_data):
            bin_tables = pd.DataFrame(i).to_string(header=True, index=False)
            f.write('\n::MarkerName::')
            f.write(str(bins[0][j]+'::'))
            f.write('\n')
            f.write(bin_tables)
        
        f.close()
        return not is_open

    return is_open
    
##-----------------------------allele calling-----------------------------##  




@app.callback(Output('alleles', 'children'),
              Output('alleles-header','children'),
              State('threshold', 'value'),
              Input('call','n_clicks'),
              State('edit-size-table', 'data'),
              State('bin','filename'),
              State('fsa','filename'),
              State('liz','filename'),
              State('select-sizes','value'),
              State('sizeladder-threshold','value'),
             
              State('project-name','value'),
              State('panel-name','value'),
              State('select-panel','value'),
              [State('bin-table{}'.format(i), 'derived_virtual_data')
               for i in range(read_bin_list())],
              [State('bin-table{}'.format(i), 'derived_virtual_selected_rows')
               for i in range(read_bin_list())],
              prevent_initial_call=True)


def call_alleles(threshold,call_clicks,data,binfile,fsa,liz,selectedsize,ladderthreshold,projectname,panel_name,panelselect,*bin_data):
    #diğer call-backte değiştirilen dataları bu call-backte kullanıyoruz. 
    #yeni call-back açılmasındaki amaç hem show-bins üzerinde yapılan değişiklikleri eş zamanlı yansıtıp
    #hem de allele calling kısmına butona basmadan önce değişiklikleri yansıtmamak, allele calling kısmında state, diğer call-backte input kullandım.
    header='Alleles Captured:'
    UPLOAD_DIRECTORY = os.path.join(os.getcwd() ,"project") 
    UPLOAD_DIRECTORY =os.path.join(UPLOAD_DIRECTORY,str(projectname))
    UPLOAD_DIRECTORY = os.path.join(UPLOAD_DIRECTORY ,"app_uploaded_files")
    if not os.path.exists(UPLOAD_DIRECTORY):
        os.makedirs(UPLOAD_DIRECTORY) 
    
    binpath=os.path.join(UPLOAD_DIRECTORY,binfile[0])
    path=os.path.join(UPLOAD_DIRECTORY, fsa[0])
    
    #######panel kısmı
    def return_panel():
        if panel_name is not None:
            panel=panel_name
        if panel_name is None and panelselect is not None:
            panel=panelselect
        return panel
    
    panel_name=return_panel()
        
    print('PANEL NAMEEEEE',panel_name)
    
    l=int(len(bin_data)/2)
    rows=bin_data[:l:1]
    
    #liz kısmı
    def return_sizes():
        if liz is not None: 
            liz_path=os.path.join(UPLOAD_DIRECTORY,liz[0]) 
            sizes=read_sizes(liz_path)
            
        if liz is None and selectedsize is not None:
            sizes = size_dict[selectedsize]
    
        return sizes
    
    sizes=return_sizes() 
    
    liz = LizStd(path.encode('unicode_escape').decode(),sizes,ladderthreshold)
 
    
 
    df = pd.DataFrame(data, columns=['status'])

    bins = read_bins(binpath,panel_name)
    size_std = list(df['status'])
    liz_peaks = list(liz[1])

    bin_data_list = []
    for i in bin_data: 
        bin_data_list.append(pd.DataFrame(i))
            
        #print('bin data list')
        #print(bin_data_list)

    sizepair = {}
    for i, j in zip(size_std, liz_peaks):
      
        if i == 'no match' or i=="":
            pass
        else:
            sizepair[int(i)] = j
   

    get = get_sizes_v2(path, list(sizepair.keys()),
                       list(sizepair.values()), threshold)
    
    
    
    #--------------------ALLELE CALLING PART--------------------#
    peak_heights=[] 
    peak_index=[] 
    peak_size=[] 
    alleles=[] 
    marker_infos=[]
    colors=[]
    tracks=[]
    for z,i in enumerate(get[4][0].keys()):
        for j in range(len(get[4][0][i])):
            size=get[4][0][i].iloc[j,5] 
            for d,a in enumerate(rows):
                try:
                    for b in a:
                        #print('size:',size)
                        #print("float(b['start_size'])")
                        #print(float(b['start_size']))
                        
                        #print("float(b['end_size'])")
                        #print(float(b['end_size']))
                        

                        if float(size) >= float(b['start_size']) and float(size) <= float(b['end_size']) and str.lower(b['color']) == dye_dict[str(z+1)]:
                            peak_size.append(size)
                            peak_index.append(get[4][0][i].iloc[j,0])
                            peak_heights.append(get[4][0][i].iloc[j,1])
                            alleles.append(b['allele'])
                            marker_infos.append(bins[0][d])
                            colors.append(dye_dict[str(z+1)])
                            tracks.append(i)
                            #print(size,b['allele'],dye_dict[str(dyes[z+1])],bins[0][d]) 
                except:
                    pass
    output={}
    output['Track_Data']=tracks
    output['Index']=peak_index
    output['Size']=peak_size
    output['Height']=peak_heights
    output['Allele']=alleles
    output['Marker']=marker_infos
    output['Color']=colors

    allele_output=pd.DataFrame(output)
    allele_output.index+=1
    
    
    
    if call_clicks != 0:
        final=dash_table.DataTable(id='allele-calling',
                         data=allele_output.to_dict('records'),
                         columns=[{'id': 'Track_Data', 'name':'Track_Data'},
                                  {'id': 'Index', 'name':'Index'},
                                  {'id': 'Size', 'name':'Size'},
                                  {'id': 'Height', 'name':'Height'},
                                  {'id': 'Allele', 'name':'Allele'},
                                  {'id': 'Marker', 'name':'Marker'},
                                  {'id': 'Color', 'name':'Color'}
                                  ],
                         style_header={
                            'backgroundColor': 'rgb(33,46,82)',
                            'color': 'white',
                            'fontWeight': 'bold'},
                         style_cell={
                             'maxWidth': '5px',
                             'whiteSpace': 'normal'}
        )
    else:
        final=[]
        
    return final,header











##-----------------------------------------------------------------------##  

@app.callback(Output('genemapper-tabs', 'value'),
              [Input('return-size-matching', 'n_clicks'),
              Input('submit-button','n_clicks'),
              Input('call','n_clicks')])
def open_home_tab(n_clicks1,n_clicks2,n_clicks3):
    #ok butonuna basıldığında step2'ye yönlendirsin.
    btn=dash.callback_context.triggered[0]["prop_id"].split(".")[0]
    #ok butonuna basıldığında step2'ye yönlendirsin.
    
    if btn=='return-size-matching':
        return 'Size Matching'
    
    if btn=='submit-button':
        return 'Size Graph'
    
    if btn=='call':
        return 'Allele Calling'
    
    else: 
        return 'Upload Data'





app.run(debug=False) 