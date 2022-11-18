from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
import os


def main():
    global group_set
    plot_folder_name='plot_TSNE'
    try:
        os.mkdir(plot_folder_name)
    except:
        print('plot folder already exist')
    df=pd.read_excel('CC_data.xlsx',sheet_name='data_check')
    df=df[(df['E_c_1']=='O')]
    #df=df[(df['LIG_TYPE']=='P;P;P')|(df['LIG_TYPE']=='C;C')]
    df['Count'] = 0
    #df=df.sort_values('LIG_POS')
    group_set=['E_BDE','Nu_BDE','S_BDE', # Bond dissociation energy
                    'PA', # Ligand parameters
                    'LIG_TYPE','S','LIG_UNAME', # Ligand name
                    'E_c_1','E_c_2', # Class of leaving group
                    'type', # Nu-Nu, E-Nu, E-E
                    'LIG','formed_bond',
                    'S_1','S_2',
                    'S_1_CNMR','S_2_CNMR',
                    'LIG_VB']
    df=df.groupby(group_set).Count.count().reset_index()
    df.sort_values(['LIG_TYPE'])

    print('df before________\n',df)

    #total_set=['E_BDE','Nu_BDE','S_BDE','EHOMO','ELUMO','PA','He8_des']
    total_set=['E_BDE','Nu_BDE','S_BDE','PA','LIG_VB','S_1_CNMR','S_2_CNMR']


        # Performing kmean
    data_df=df[total_set]#,'LIG_DEN']]
    data_df=(data_df-data_df.min())/(data_df.max()-data_df.min())
    #pca = PCA(n_components=3)
    kmeans = KMeans(n_clusters=5, random_state=0).fit(data_df)
    label_list=kmeans.labels_
    
        # check data size
    print('len of label',len(label_list))
    print('len of data_df',len(data_df))
    print('len of df',len(df))

    



        # Configue TSNE

    learning_rate_set=500
    perplexity_set=50
    random_state_set=50
    n_iter_set=1000
        # Reduction
    X_transformed = TSNE(n_components=2,
    learning_rate=learning_rate_set,
    perplexity=perplexity_set,
    random_state=random_state_set,
    n_iter=n_iter_set).fit_transform(data_df)
    
    # Get tag
    tag='2D'+'_l'+str(learning_rate_set)+'_p'+str(perplexity_set)+'_r'+str(random_state_set)+'_i'+str(n_iter_set)
    plot_path=plot_folder_name+'/'+tag

    try:
        os.mkdir(plot_folder_name+'/'+str(tag))
    except:
        pass

        # assigned label
    label_list=['class ' + str(i+1) for i in label_list]
    df['label']=label_list
    df.to_excel(plot_path+'/data_w_label.xlsx')

        
    print('X_embedded____\n',X_transformed)
    #plot_scatter_2D(X_embedded, df,'label')
    for item in total_set:
        plot_scatter_2D(X_transformed, df,item,plot_path)


    


    plot_scatter_2D(X_transformed, df,'formed_bond',plot_path)
    plot_scatter_2D(X_transformed, df,'LIG_TYPE',plot_path)
    plot_scatter_2D(X_transformed, df,'LIG_UNAME',plot_path)
    plot_scatter_2D(X_transformed, df,'E_c_1',plot_path)
    plot_scatter_2D(X_transformed, df,'E_c_2',plot_path)
    plot_scatter_2D(X_transformed, df,'S_1',plot_path)
    plot_scatter_2D(X_transformed, df,'S_2',plot_path)
    plot_scatter_2D(X_transformed, df,'type',plot_path)
    plot_scatter_2D(X_transformed, df,'label',plot_path)

    # 3D ploting
'''
    # 3D ploting
    X_transformed = TSNE(n_components=3,
    learning_rate=learning_rate_set,
    perplexity=perplexity_set,
    random_state=random_state_set).fit_transform(data_df)
    
    # Configue TSNE
    learning_rate_set=1000
    perplexity_set=60
    random_state_set=0
    tag='3D_p60_l1000'
    plot_path=plot_folder_name+'/'+tag
 
    try:
        os.mkdir(plot_folder_name+'/'+str(tag))
    except:
        pass

    
    #print('X_embedded____\n',X_transformed.shapy,X_transformed)

    
    for item in total_set:
        plot_scatter_3D(X_transformed, df,item,plot_path)


    plot_scatter_3D(X_transformed, df,'formed_bond',plot_path)
    plot_scatter_3D(X_transformed, df,'LIG_TYPE',plot_path)
    plot_scatter_3D(X_transformed, df,'LIG_UNAME',plot_path)
    plot_scatter_3D(X_transformed, df,'E_c_1',plot_path)
    plot_scatter_3D(X_transformed, df,'E_c_2',plot_path)
    plot_scatter_3D(X_transformed, df,'S_1',plot_path)
    plot_scatter_3D(X_transformed, df,'S_2',plot_path)
    plot_scatter_3D(X_transformed, df,'type',plot_path)
'''



##=============== Plot function
def plot_scatter_2D(x,df,color_col,plot_path,color_discrete_sequence_set=px.colors.qualitative.Dark24):
    global group_set
    fig = px.scatter(df, x=x[:,0], y=x[:,1],
        #color=df['E_c_1'].astype(str),
        color=df[color_col],
        hover_data=group_set,
        size=df['Count']+5,
        size_max=10,
        color_discrete_sequence=color_discrete_sequence_set,
        title='<b>'+color_col+'</b>'

        )
    fig.update_layout(scene = dict(
                        xaxis_title='TSNE 1',
                        yaxis_title='TSNE 2'))
    fig.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1
    )
    fig.write_html(plot_path+'/'+'[2D]'+color_col+".html")
    

def plot_scatter_3D(x,df,color_col,plot_path,color_discrete_sequence_set=px.colors.qualitative.Dark24):
    fig = px.scatter_3d(df, x=x[:,0], y=x[:,1], z=x[:,2],
        #color=df['E_c_1'].astype(str),
        color=df[color_col],
        hover_data=['S','E_c_1','E_c_2','LIG_TYPE'],
        size=df['Count']+5,
        size_max=10,
        color_discrete_sequence=color_discrete_sequence_set,
        title='<b>'+color_col+'</b>'

        )
    fig.update_layout(scene = dict(
                        xaxis_title='TSNE 1',
                        yaxis_title='TSNE 2',
                        zaxis_title='TSNE 3'))
    fig.update_yaxes(
        scaleanchor = "x",
        scaleratio = 1
    )
    fig.write_html(plot_path+'/'+'[3D]'+color_col+".html")
   
#=======main=======
if __name__ == '__main__':
    main()
