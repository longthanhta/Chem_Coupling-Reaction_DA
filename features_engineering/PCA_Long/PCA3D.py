from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn import preprocessing


##=============== Plot function
def plot_scatter(x,df):
    with open('p_graph.html', 'a') as f:
        fig = px.scatter_3d(df, x=x[:,0], y=x[:,1], z=x[:,2],
            #color=df['class'].astype(str),
            text=df['LIG'],
            color=df['yield'],
            hover_data=['ID','EGC','SCFA','SCF','LIG','yield','class'],
            symbol=df['LIG_GROUP'],
            )
        #fig.update_traces(marker=dict(size=10,
        #                      line=dict(width=2,
        #                                color='DarkSlateGrey')),
        #          selector=dict(mode='markers'))
        fig.update_traces(
            textposition='top center',
            textfont=dict(
            family="sans serif",
            #size=8,
            ),
            line=dict(width=2,
                    color='black'),
            )
        fig.update_layout(
            title = "fixed-ratio axes with compressed axes",
            xaxis = dict(
              scaleanchor = "y",
              scaleratio = 1,
            ),
            yaxis = dict(
              scaleanchor = "x",
              scaleratio = 1,
            ),
            legend= dict(
                traceorder="grouped",
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.01
                ),
            scene = dict(
                    xaxis_title='PCA1',
                    yaxis_title='PCA2',
                    zaxis_title='PCA3'),
        )
        fig.write_html("PCA3D_plot_yield.html")
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
        fig = px.scatter_3d(df, x=x[:,0], y=x[:,1], z=x[:,2],
            text=df['EGC'],
            color=df['class'].astype(str),
            #color=df['yield'],
            hover_data=['ID','EGC','SCFA','SCF','LIG','yield','class'],
            symbol=df['LIG_GROUP']
            #size=df['size']
            )
        fig.update_traces(
            textposition='top center',
            textfont=dict(
            family="sans serif",
            #size=8,
            ),
            line=dict(width=2,
                    color='black'),
            )
        fig.update_layout(
            title = "fixed-ratio axes with compressed axes",
            xaxis = dict(
              scaleanchor = "y",
              scaleratio = 1,
            ),
            yaxis = dict(
              scaleanchor = "x",
              scaleratio = 1,
            ),
            legend= dict(
                traceorder="grouped",
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.01
                ),
            scene = dict(
                    xaxis_title='PCA1',
                    yaxis_title='PCA2',
                    zaxis_title='PCA3'),
        )
        fig.write_html("PCA3D_plot.html")
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
        fig = px.scatter_3d(df, x=x[:,0], y=x[:,1], z=x[:,2],
            color=df['class'].astype(str),
            hover_data=['EGC','SCFA','SCF','LIG','yield'],
            symbol=df['LIG_GROUP']
            #size=df['size']
            )
        fig.update_traces(
            textposition='top center',
            textfont=dict(
            family="sans serif",
            #size=8,
            ),
            line=dict(width=2,
                    color='black'),
            )
        fig.update_layout(
            title = "fixed-ratio axes with compressed axes",
            xaxis = dict(
              scaleanchor = "y",
              scaleratio = 1,
            ),
            yaxis = dict(
              scaleanchor = "x",
              scaleratio = 1,
            ),
            legend= dict(
                traceorder="grouped",
                yanchor="top",
                y=0.99,
                xanchor="right",
                x=0.01
                ),
            scene = dict(
                    xaxis_title='PCA1',
                    yaxis_title='PCA2',
                    zaxis_title='PCA3'),
        )
        fig.write_html("PCA3DLigand.html")
        f.write(fig.to_html(full_html=False, include_plotlyjs='cdn'))
#print('dists',dists)
dists = pickle.load(open("bits.pkl", "rb"))
#dists = pickle.load(open("bitsMG.pkl", "rb"))
#'''
x = np.array(dists) #returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
dists = pd.DataFrame(x_scaled)
''

df = pd.read_excel('out.xlsx')
pca = PCA(n_components=3)
print('shape of df',df.shape)

X_transformed=pca.fit_transform(dists)
plot_scatter(X_transformed, df)
pcac=pca.components_
PCA_coef=pd.DataFrame(pcac)
PCA_coef=PCA_coef.abs()
header_df=pd.read_csv('header_lst.csv')
header_list=header_df['Item'].values.tolist()
PCA_coef.to_csv('PCA_coef.csv',index=None,header=header_list)
