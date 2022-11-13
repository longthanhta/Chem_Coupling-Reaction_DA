from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import plotly.express as px
from sklearn.decomposition import PCA


##=============== Plot function
def plot_scatter(x,df):
    fig = px.scatter(df, x=x[:,0], y=x[:,1],
        color=df['class'].astype(str),
        hover_data=['EGC','SCFA','SCF'],
        )
    fig.update_traces(
        textposition='top center',
        textfont=dict(
        family="sans serif",
        size=12,
        )
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
            traceorder="grouped"
            )
    )
    fig.write_html("PCA_plot.html")
    fig2 = px.scatter(df, x=x[:,0], y=x[:,1],
        color=df['class'].astype(str),
        hover_data=['EGC','SCFA','SCF'],
        )
    fig2.update_traces(
        textposition='top center',
        textfont=dict(
        family="sans serif",
        size=12,
        )
        )
    fig2.update_layout(
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
            traceorder="grouped"
            )
    )
    fig2.write_html("PCA_plot_wo_text.html")

#print('dists',dists)
dists = pickle.load(open("bits.pkl", "rb"))
#dists = pickle.load(open("bitsMG.pkl", "rb"))

df = pd.read_excel('out.xlsx')
pca = PCA(n_components=2)
print('shape of df',df.shape)

X_transformed=pca.fit_transform(dists)
plot_scatter(X_transformed, df)
pcac=pca.components_
#element_list=dists.columns.values
#PCA_coef=pd.DataFrame(pcac,columns=element_list[1:])
#PCA_coef.to_csv('PCA_coef.csv')
