import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np; np.random.seed(42)
import os
import pandas as pd
import matplotlib.patches as mpatches
from colorhash import ColorHash

def add_identity(axes, *line_args, **line_kwargs):
    identity, = axes.plot([], [], *line_args, **line_kwargs)
    def callback(axes):
        low_x, high_x = axes.get_xlim()
        low_y, high_y = axes.get_ylim()
        low = max(low_x, low_y)
        high = min(high_x, high_y)
        identity.set_data([low, high], [low, high])
    callback(axes)
    axes.callbacks.connect('xlim_changed', callback)
    axes.callbacks.connect('ylim_changed', callback)
    return axes

def get_color(item):
    return ColorHash(item).hex




if __name__ == '__main__':
    param_vec = pd.read_csv('params.csv',header=None).values.flatten()
    predtrain_vec = pd.read_csv('predtrain.csv',header=None).values.flatten()
    df_data = pd.read_excel('CC_data.xlsx',sheet_name='data_check')
    df_data=df_data[(df_data['E_c_1']=='O')]
    features_list=['S_2_CNMR','Nu_BDE','E_HOMO','E_LUMO','S_1_CNMR','PA','S_BDE','LIG_VB']
    bit_df=data_df[features_list]
    bit_np=bit_df.values
    yield_vec=data_df['E_BDE']
    header_lst=list(bit_df.columns.values)
    #print('bit_np\n',bit_np)
    #print('param_vec\n',len(param_vec))
    #print('yield_vec\n',len(yield_vec))
    # rescale param
    combine_lst_scaled=[]
    for index, row in enumerate(bit_np):
        sum=0
        for id2 in range(0,len(row)-1):
            sum+=row[id2]*param_vec[id2]
        combine_lst_scaled.append(sum)
    lig_list=[]
    for index, row in data_df.iterrows():
        #print(row)
        lig=row['label']
        lig_list.append(lig)
    min_v=min(combine_lst_scaled)
    min_i=combine_lst_scaled.index(min_v)
    max_v=max(combine_lst_scaled)
    max_i=combine_lst_scaled.index(max_v)
    print('min',min_v,'at',min_i)
    print('predicted yield at',min_i,'is',predtrain_vec[min_i])
    print('maxn',max_v,'at',max_i)
    print('predicted yield at',max_i,'is',predtrain_vec[max_i])
    range_estimated=max_v-min_v
    range_predicted=predtrain_vec[max_i]-predtrain_vec[min_i]
    unscale_fract=range_predicted/range_estimated
    unscale_fract=unscale_fract*0.8
    low_limit=predtrain_vec[min_i]
    high_limit=predtrain_vec[max_i]
    print('low_limit',low_limit)
    
    
    
    offset=0
    combine_lst=combine_lst_scaled
    combine_lst=[]
    combine_label=[]
    for index, row in enumerate(bit_np):
        sum=0
        for id2 in range(0,len(row)-1):
            sum+=row[id2]*param_vec[id2]*unscale_fract
        combine_lst.append(sum+offset)
    count=0
    scaled_coef=[]
    max_val_lst=[]
    min_val_lst=[]
    true_coef=[]
#'''


    new_header_lst=[]
    for ids,header in enumerate(header_lst):
        sc=param_vec[ids]*unscale_fract
        max_v=bit_df[header_lst[ids]].max()
        min_v=bit_df[header_lst[ids]].min()
        if np.abs(sc*(max_v-min_v)) == 0: continue
        new_header_lst.append(header)
        combine_label.append(str(param_vec[ids]*unscale_fract)[0:7])
        scaled_coef.append(sc)
        max_val_lst.append(max_v)
        min_val_lst.append(min_v)
        combine_label.append('.')
        combine_label.append(header_lst[ids])
        combine_label.append(' + ')
        true_coef.append(np.abs(sc*(max_v-min_v)))
        count+=1
        if count>4:
            combine_label.append('\n')
            count=0
    new_combine_lst=[]
    for index, row in enumerate(bit_np):
        sum=0
        for id2 in range(0,len(row)-1):
            if header_lst[id2] not in new_header_lst: continue
            sum+=row[id2]*param_vec[id2]*unscale_fract
        new_combine_lst.append(sum+offset)
    combine_label.append(str(offset)[0:7])
    combine_label=''.join(combine_label)
    scaled_df=pd.DataFrame({'features':new_header_lst,'coef':scaled_coef,'max':max_val_lst,'min':min_val_lst,'coef*(max-min)':true_coef})
    scaled_df.to_excel('scaled_df.xlsx',index=None)
    within_range=True
    for index,item in enumerate(new_combine_lst):
        if not(low_limit<item<high_limit):
        #if not(low_limit<item):
            within_range=False
    if within_range:
        condition=True

    # plot part:


    x=new_combine_lst
    y=yield_vec
    fig = plt.figure()
    ax = fig.add_subplot(111)
    unique_label=[]
    for index in range(0,len(x)):
        lig=lig_list[index]
        print('lig',lig)
        if lig not in unique_label:
            plt.scatter(x[index],y[index], s=100,c=get_color(lig),label=lig)
            unique_label.append(lig)
        else:
            plt.scatter(x[index],y[index], s=100,c=get_color(lig))




    # Draw the black line
    xl=np.arange(-1,1)
    #plt.plot(xl,xl,'k-')
    #plt.xlim(30,100)
    #plt.ylim(30,100)


    fig.set_size_inches(10, 10)

    #plt.title('log(RR) prediction META')
    plt.ylabel('actual', fontsize=16)



    plt.xlabel('predicted  \n'+combine_label, fontsize=10)
    #for i in range(len(x)):
    #    label=data_df.iloc[i]['entry']
    #    ax.annotate(label, (x[i], y[i]),color='black',fontsize=9,size=13)
    #os.chdir('..')
    #plt.legends()
    add_identity(ax, color='k', ls='--')
    ax.legend()
    #plt.show()
    plt.savefig('plot.png')
#'''#
