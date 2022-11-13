import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np; np.random.seed(42)
import os
import pandas as pd
import matplotlib.patches as mpatches

def get_color(egc):
    if egc == 'Br.B' : return '#F4D03F'
    if egc == 'Br.Br' : return '#78281F'
    if egc == 'Br.H' : return '#C39BD3'
    if egc == 'C.B' : return '#F39C12'
    if egc == 'Cl.B' : return '#E67E22'
    if egc == 'Cl.In' : return '#273746'
    if egc == 'F.Al' : return '#C39BD3'
    if egc =='F.B' : return '#D35400'
    if egc == 'I.B' : return '#C0392B'
    if egc == 'I.H' : return '#A3E4D7'
    if egc == 'N.B' : return '#E74C3C'
    if egc == 'O.Al' : return '#4D5656'
    if egc == 'O.B' : return '#F6DDCC'
    if egc == 'O.H' : return '#A9CCE3'
    if egc == 'O.Mg' : return '#5499C7'
    if egc == 'O.Sn' : return '#AED6F1'
    if egc == 'S.B' : return '#78281F'



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
    data_df = pd.read_excel('CC_CC_data_Ph_Ph_yield_sub.xlsx')
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
        print(row)
        EGC=row['EG_class_display']
        color_lst.append(get_color(EGC))
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
    offset=predtrain_vec[min_i]-min_v
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
    for ids in range(len(header_lst)):
        scaled_coef.append(param_vec[ids]*unscale_fract)
        if param_vec[ids] == 0: continue
        combine_label.append(str(param_vec[ids]*unscale_fract)[0:7])
        combine_label.append('.')
        combine_label.append(header_lst[ids])
        combine_label.append(' + ')
        count+=1
        if count>4:
            combine_label.append('\n')
            count=0
    combine_label.append(str(offset)[0:7])
    combine_label=''.join(combine_label)
#'''
    # plot part:
    os.chdir('1_reactions')
    fig_name_lst=range(len(data_df))
    fig_name_lst= [str(v)+'.png' for v in fig_name_lst]
    fig_name_lst = np.array(fig_name_lst).astype('<U12')
    x=combine_lst
    y=yield_vec
    fig = plt.figure()
    ax = fig.add_subplot(111)
    line = plt.scatter(x,y, s=100,c=color_lst)
    image_path = np.asarray(fig_name_lst)


    # Draw the black line
    xl=np.arange(0,100)
    plt.plot(xl,xl,'k-')
    plt.xlim(0,100)
    plt.ylim(0,100)

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

    plt.title('PCy3')
    plt.xlabel('predicted yield \n'+combine_label, fontsize=10)
    plt.ylabel('yield', fontsize=16)
    label_lg_lst=list(set(data_df['EG_class_display'].values.tolist()))
    recs = []
    for i in range(0,len(label_lg_lst)):
        recs.append(mpatches.Rectangle((0,0),1,1,fc=get_color(label_lg_lst[i])))
    plt.legend(recs,label_lg_lst,loc=2)
    scaled_df=pd.DataFrame({'features':header_lst,'coef':scaled_coef})
    scaled_df.to_csv('scaled_df.csv',index=None)
    plt.show()
#'''
