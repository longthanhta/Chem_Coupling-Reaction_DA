from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import plotly.express as px
from threading import Thread
from multiprocessing import Process, Manager
from math import sqrt, floor


##=============== Plot function
def plot_scatter(x,df):
    fig = px.scatter(df, x=x[:,0], y=x[:,1],
        color=df['FORM(upper)'].astype(str),hover_data=['REF','AUT','ID','BREAK','LEAVE','CHANGE','FORM','Count'])
    fig.write_html("plot_manual_clusters.html")


#print('dists',dists)
dists = pickle.load(open("distance_matrix.pkl", "rb"))
df = pd.read_csv('data_for_plot.csv')
X_tsne = TSNE(random_state=1,metric='precomputed').fit_transform(dists)
#print('X_tsne',X_tsne,len(X_tsne))
plot_scatter(X_tsne, df)
