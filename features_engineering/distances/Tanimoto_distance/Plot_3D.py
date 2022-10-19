from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import plotly.express as px



##=============== Plot function
def plot_scatter(x,df):
    fig = px.scatter_3d(df, x=x[:,0], y=x[:,1], z=x[:,2],
    	color_discrete_sequence=px.colors.qualitative.Light24,
        color=df['FORM(upper)'].astype(str),hover_data=['REF','AUT','ID','BREAK','LEAVE','CHANGE','FORM','Count'])
    fig.write_html("plot_manual_clusters.html")


#print('dists',dists)
dists = pickle.load(open("distance_matrix.pkl", "rb"))
df = pd.read_csv('data_for_plot.csv')
X_tsne = TSNE(random_state=1,n_components=3,metric='precomputed',perplexity=30,n_iter=2000).fit_transform(dists)
print('X_tsne',X_tsne,len(X_tsne))
plot_scatter(X_tsne, df)
