from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans


##=============== Plot function
def plot_scatter_3D(x,df,color_col,color_discrete_sequence_set=px.colors.qualitative.Dark24):
    fig = px.scatter_3d(df, x=x[:,0], y=x[:,1],z=x[:,2],
        color=df[color_col].astype(str),
        hover_data=['S','E_c_1','E_c_2','LIG_TYPE','LIG_UNAME'],
        size=df['Count']+2,
        size_max=30,
        color_continuous_scale='Plasma',
        color_discrete_sequence=color_discrete_sequence_set,
        #text=df['LIG_TYPE']
        title='<b>'+color_col+'</b>'
        )
    fig.update_layout(scene = dict(
                        xaxis_title='PCA 1',
                        yaxis_title='PCA 2',
                        zaxis_title='PCA 3',
                        xaxis=dict(backgroundcolor='white',gridcolor="#A19882"),
                        yaxis=dict(backgroundcolor='white',gridcolor="#A19882"),
                        zaxis=dict(gridcolor="#A19882"),
                        ))
    fig.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1
    )
    # projection
    
    
    
    
    fig.write_html('plot/[3D]'+color_col+"_3DPCA_plot.html")


def plot_scatter_2D(x,df,color_col,color_discrete_sequence_set=px.colors.qualitative.Dark24):
    fig = px.scatter(df, x=x[:,0], y=x[:,1],
        #color=df['E_c_1'].astype(str),
        color=df[color_col].astype(str),
        hover_data=['S','E_c_1','E_c_2','LIG_TYPE'],
        size=df['Count']+5,
        size_max=10,
        color_discrete_sequence=color_discrete_sequence_set,
        title='<b>'+color_col+'</b>'

        )
    fig.update_layout(scene = dict(
                        xaxis_title='PCA 1',
                        yaxis_title='PCA 2'))
    fig.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1
    )
    fig.write_html('plot/[2D]'+color_col+"PCA_plot.html")

df=pd.read_excel('CC_data.xlsx',sheet_name='data_check')
df=df[(df['type']=='E-Nu')]
#df=df[(df['LIG_TYPE']=='P;P;P')|(df['LIG_TYPE']=='C;C')]
x_col='E_BDE'
z_col='S_BDE'
x_label_col='ESmi1'
y_label_col='LIG_TYPE'
z_label_col='form_bond'
df['Count'] = 0
#df=df.sort_values('LIG_POS')
df=df.groupby(['E_BDE','Nu_BDE','S_BDE', # Bond dissociation energy
                'PA','He8_des', # Ligand parameters
                'LIG_TYPE','S','LIG_UNAME', # Ligand name
                'E_c_1','E_c_2', # Class of leaving group
                'type', # Nu-Nu, E-Nu, E-E
                'LIG','form_bond',
                'EHOMO','ELUMO',
                'S_1','S_2',
                'S_1_CNMR','S_2_CNMR',
                'LIG_VB1','LIG_VB2',
                ]).Count.count().reset_index()
df.sort_values(['LIG_TYPE'])

print('df before________\n',df)

#total_set=['E_BDE','Nu_BDE','S_BDE','EHOMO','ELUMO','PA','He8_des']
total_set=['E_BDE','Nu_BDE','S_BDE','EHOMO','ELUMO','PA','LIG_VB2','S_1_CNMR','S_2_CNMR']



data_df=df[total_set]#,'LIG_DEN']]
data_df=(data_df-data_df.min())/(data_df.max()-data_df.min())
#pca = PCA(n_components=3)
kmeans = KMeans(n_clusters=3, random_state=0).fit(data_df)
label_list=kmeans.labels_
print('len of label',len(label_list))
print('len of data_df',len(data_df))
print('len of df',len(df))

df['label']=label_list
#data_df['label']=label_list
data_df.to_excel('plot_data.xlsx')



# 2D ploting

pca = PCA(n_components=2)
X_transformed=pca.fit_transform(data_df)


plot_scatter_2D(X_transformed, df,'label')

#pca_score=pca.components_
#score_df=pd.DataFrame({'comp':total_set,'pca1':pca_score[0],'pca2':pca_score[1]})
#score_df.to_excel('plot/score2D.xlsx')



pca = PCA(n_components=3)
X_transformed=pca.fit_transform(data_df)


plot_scatter_3D(X_transformed, df,'label')

#pca_score=pca.components_
#score_df=pd.DataFrame({'comp':total_set,'pca1':pca_score[0],'pca2':pca_score[1],'pca3':pca_score[2]})
#score_df.to_excel('plot/score3D.xlsx')




