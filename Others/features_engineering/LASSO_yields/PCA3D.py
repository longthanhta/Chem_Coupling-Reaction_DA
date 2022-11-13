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
    fig = px.scatter_3d(df, x=x[:,0], y=x[:,1], z=x[:,2],
        #color=df['class'].astype(str),
        #text=df['EG_class'],
        color=df['EG_class'],
        hover_data=['OID','EG_class','SCF_gen_uni'],
        #symbol=df['LIG_GROUP'],
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
    fig.write_html("PCA3D_plot.html")

dists_df = pd.read_csv("ft.csv")
#dists = pickle.load(open("bitsMG.pkl", "rb"))
#'''
x = np.array(dists_df) #returns a numpy array
min_max_scaler = preprocessing.MinMaxScaler()
x_scaled = min_max_scaler.fit_transform(x)
dists = pd.DataFrame(x_scaled)
''

df = pd.read_excel('CC_CC_data_Ph_Ph_yield_sub.xlsx')
pca = PCA(n_components=3)
print('shape of df',df.shape)

X_transformed=pca.fit_transform(dists)
plot_scatter(X_transformed, df)
pcac=pca.components_
PCA_coef=pd.DataFrame(pcac)
PCA_coef=PCA_coef.abs()
PCA_coef.to_csv('PCA_coef.csv',index=None,header=dists_df.columns.values)
