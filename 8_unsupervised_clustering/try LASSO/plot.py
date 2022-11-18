import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np; np.random.seed(42)
import os
import pandas as pd
import matplotlib.patches as mpatches
from colorhash import ColorHash
def get_color(egc):
    #return ColorHash(egc).hex
    if egc=='B.C': return 'green'
    if egc=='C.Si': return 'blue'



def hover(event):
    # if the mouse is over the scatter points
    if line.contains(event)[0]:
        # find out the index within the array from the event
        ind, = line.contains(event)[1]["ind"]
        # get the figure size
        w,h = fig.get_size_inches()*fig.dpi
        ws = (event.x > w/2.)*-1 + (event.x <= w/2.)
        hs = (event.y > h/2.)*-1 + (event.y <= h/2.)
        # if event occurs in the top or right quadrant of the figure,
        # change the annotation box position relative to mouse.
        ab.xybox = (xybox[0]*ws, xybox[1]*hs)
        # make annotation box visible
        ab.set_visible(True)
        # place it at the position of the hovered scatter point
        ab.xy =(x[ind], y[ind])
        # set the image corresponding to that point
        im.set_data(plt.imread(image_path[ind]))
    else:
        #if the mouse is not over a scatter point
        ab.set_visible(False)
    fig.canvas.draw_idle()





if __name__ == '__main__':
    bit_df=pd.read_csv('ft.csv')
    bit_np = pd.read_csv('ft.csv').values #input file
    param_vec = pd.read_csv('params.csv',header=None).values.flatten()
    predtrain_vec = pd.read_csv('predtrain.csv',header=None).values.flatten()
    data_df = pd.read_excel('data_w_label_2D_l5000_p30_r200_i1000.xlsx',sheet_name='4 Vs 6')
    df_features=df_data[['E_BDE','Nu_BDE','S_1_CNMR','S_2_CNMR','LIG_VB','PA']]
    yield_vec=data_df.Yield
    header_lst=list(bit_df.columns.values)
    #print('bit_np\n',bit_np)
    #print('param_vec\n',len(param_vec))
    #print('yield_vec\n',len(yield_vec))
    # rescale param
    combine_lst_scaled=[]
    color_lst=[]
    for index, row in enumerate(bit_np):
        sum=0
        for id2 in range(0,len(row)-1):
            sum+=row[id2]*param_vec[id2]
        combine_lst_scaled.append(sum)
    for index, row in data_df.iterrows():
        #print(row)
        EGC=row['EG_class_display']
        OID=row['OID']
        color_lst.append('green')
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
    offset=400
    print('offset',offset)
    offset=-1000
    low_limit=predtrain_vec[min_i]-10
    print('low_limit',low_limit)
    condition=False
    while not condition:
        offset+=20
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
        scaled_df.to_csv('scaled_df.csv',index=None)
        within_range=True
        for index,item in enumerate(new_combine_lst):
            if not(low_limit<=item<=150):
                within_range=False
        if within_range:
            condition=True


    # plot part:


    os.chdir('1_reactions')
    fig_name_lst=range(len(data_df))
    fig_name_lst= [str(v)+'.png' for v in fig_name_lst]
    fig_name_lst = np.array(fig_name_lst).astype('<U12')
    x=new_combine_lst
    y=yield_vec
    fig = plt.figure()
    ax = fig.add_subplot(111)
    line = plt.scatter(x,y, s=100,c=color_lst)
    image_path = np.asarray(fig_name_lst)



    # Draw the black line
    xl=np.arange(0,100)
    plt.plot(xl,xl,'k-')
    #plt.xlim(30,100)
    #plt.ylim(30,100)

    # create the annotations box
    image = plt.imread(image_path[0])
    im = OffsetImage(image, zoom=0.6)
    xybox=(90., 90.)
    ab = AnnotationBbox(im, (0,0), xybox=xybox, xycoords='data',
            boxcoords="offset points",  pad=0.3,  arrowprops=dict(arrowstyle="->"))
    # add it to the axes and make it invisible
    ax.add_artist(ab)
    ab.set_visible(False)
    # add callback for mouse moves
    fig.canvas.mpl_connect('motion_notify_event', hover)

    fig = plt.gcf()
    fig.set_size_inches(10.5, 9.5)

    plt.title('PPh2Me')
    plt.ylabel('yield', fontsize=16)

    label_lg_lst=list(set(data_df['EG_class_display'].values.tolist()))
    recs = []
    #os.chdir('..')
    #for i in range(0,len(label_lg_lst)):
    #    recs.append(mpatches.Rectangle((0,0),1,1,fc='get_color(label_lg_lst[i])'))
    #plt.legend(recs,label_lg_lst,loc=2)
    print(len(header_lst),len(scaled_coef),len(max_val_lst),len(min_val_lst))


    plt.xlabel('predicted yield \n'+combine_label, fontsize=10)
    for i in range(len(x)):
        label=data_df.iloc[i]['OID']
        ax.annotate(label, (x[i]-1, y[i]-2.6),color='red',fontsize=9,size=13)
    #os.chdir('..')
    plt.show()
    #plt.savefig('plot.png')
#'''#
