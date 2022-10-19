import pandas as pd
import scipy.stats
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np


input_path = 'result_16_3_2021_Ar6nonHet-Ar6nonHet_SubINFO_HetINFO_SubDescriptor_LigDescriptor_LGSmiles_NuEINFO_NuESubDescriptor_AAM_REF.xlsx'
df = pd.read_excel(input_path)


column_list = list(df.columns)

features1=list(['Nu_p_electronic',
                'Nu_m_electronic',
                'Nu_o_electronic',
                'Nu_p_steric',
                'Nu_m_steric',
                'Nu_o_steric',
                'Dentate_Num',
                'Rotatable_Boolean', 
                'donor_EN_AVG',
                'donor_hybrid_AVG',
                'SigmaDonatingDescriptor',
                'StericWholeDescriptor',
                'Steric3DDescriptor'])

x_label = 'Dentate_Num'
x = np.array([ i[0] for i in df[[x_label]].values.tolist()])
y_label = 'Nu_m_electronic'
y = np.array([ i[0] for i in df[[y_label]].values.tolist()])
OID_list = np.array([ i[0] for i in df[['OID']].values.tolist()])



def hover(event):
    # if the mouse is over the scatter points
    if line.contains(event)[0]:
        # find out the index within the array from the event
        data_pts = line.contains(event)[1]["ind"]
        ind = data_pts[0]
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
        plt.annotate(len(data_pts), (x[ind], y[ind]),color='red')
        #im.set_data(plt.text(x[ind], y[ind],str(len(data_pts))))
    else:
        #if the mouse is not over a scatter point
        ab.set_visible(False)
    fig.canvas.draw_idle()


slope, intercept, r, p, stderr = scipy.stats.linregress(x, y)
spearman_corr , uncorr_p_value = scipy.stats.spearmanr(x,y)
line_ = f'Regression line: y={intercept:.2f}+{slope:.2f}x, r_square={r**2:.6f}, spearman corr={spearman_corr:.2f}'
fig, ax = plt.subplots()


# get image

fig_name_lst= [str(v)+'.png' for v in OID_list]
fig_name_lst = np.array(fig_name_lst).astype('<U12')
image_path = np.asarray(fig_name_lst)
image = plt.imread(image_path[0])

im = OffsetImage(image, zoom=0.6)
xybox=(90., 90.)
ab = AnnotationBbox(im, (0,0), xybox=xybox, xycoords='data',
        boxcoords="offset points",  pad=0.3,  arrowprops=dict(arrowstyle="->"))

ax.add_artist(ab)
ab.set_visible(False)
fig.canvas.mpl_connect('motion_notify_event', hover)
fig = plt.gcf()
fig.set_size_inches(10.5, 9.5)



#ax.plot(x, y, linewidth=0, marker='s', label='Data points')
line = plt.scatter(x,y, s=30)
ax.plot(x, intercept + slope * x, label=line_)
ax.set_xlabel(x_label)
ax.set_ylabel(y_label)
ax.legend(facecolor='white')
plt.show()
