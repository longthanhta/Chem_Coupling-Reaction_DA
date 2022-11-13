import matplotlib.pyplot as plt
import matplotlib.image as mpimg
# import numpy as np

def imshow(img, cfboxes = [], tboxes = []):
    plt.imshow(img)
    
    def draw_box(box, color='red'):
        x1, y1, x2, y2 = box
        w = x2 - x1
        h = y2 - y1
        plt.gca().add_patch(
            plt.Rectangle((x1, y1), w, h,
                          fill=False, edgecolor=color, linewidth=2, alpha=0.5)
        )
       
    for box in cfboxes:
        draw_box(box, 'red')
    for box in tboxes:
        draw_box(box, 'blue')
    plt.plot()
    plt.show()


########### Load hand-extracted data to run algorithm ###########
data = {'catal-4': { 
			'cfboxes' : [[11,185,115,252],[125,180,276,243],[287,173,387,246],[394,155,531,257],[223,635,507,704]], 
			'tboxes' : [[20,258,101,275],[175,260,254,277],[279,256,353,273],[411,259,472,274],[429,680,511,697],[210,589,279,605]]
			},
		'anie-1': { 
			'cfboxes' : [[77,217,149,285],[402,211,536,275],[403,300,536,369],[403,823,541,889]], 
			'tboxes' : [[425,279,518,294],[424,372,522,387],[426,892,531,905]]
			}
		}


for img in data.keys():
	table = mpimg.imread(f'../database/{img}.png')
	# print(table, data[img]['cfboxes'], data[img]['tboxes'])
	imshow(table, data[img]['cfboxes'], data[img]['tboxes'])


import math
def distance_b(box1, box2):
    dx = max(0, box2[0]-box1[2], box1[0]-box2[2])
    dy = max(0, box2[1]-box1[3], box1[1]-box2[3])
    return math.sqrt(dx*dx + dy**2)
def distance_p(point1, point2):
    dx = point2[0] - point1[0]
    dy = point2[1] - point1[1]
    return math.sqrt(dx*dx + dy**2)
def angle_score(origin, sidepoint):
    x = sidepoint[0] - origin[0]
    y = sidepoint[1] - origin[1]
    angle = math.asin(abs(y)/math.sqrt(x**2+y**2))/math.pi # in [0,0.5]
    if x < 0: angle = 1 - angle
    if y < 0: angle = 2 - angle

    ## scoring by gradient
    if angle > 1.6: return 1 - (2-angle)/(2-1.75)
    elif angle > 0.8: return 1 - ((angle-0.8)/0.8)**2
    else: return 1

import numpy as np
catal_match = []

for img in ['anie-1']:
    table = mpimg.imread(f'../database/{img}.png')
    ncf = len(data[img]['cfboxes'])
    ntx = len(data[img]['tboxes'])

    for i in range(ncf):
        cbox = data[img]['cfboxes'][i]
        ccenter = [(cbox[0]+cbox[2])/2, (cbox[1]+cbox[3])/2]
        half_max_edge = max((cbox[2]-cbox[0])/2, (cbox[3]-cbox[1])/2)
        perimeter = cbox[2]-cbox[0]+cbox[3]-cbox[1]

        dist = []
        '''
        j = np.argmin(dist)
        catal_match.append([i,j])
        table = mpimg.imread(f'../database/{img}.png')
        imshow(table, data[img]['cfboxes'][i:i+1], data[img]['tboxes'][j:j+1])
        '''

        for j in range(ntx):
            tbox = data[img]['tboxes'][j]
            tcenter = [(tbox[0]+tbox[2])/2, (tbox[1]+tbox[3])/2] # can use numpy arithmetic for reuse
            d = distance_b(cbox, tbox) + distance_p(ccenter, tcenter)
            if d < perimeter+10: # 1. constraint on distance
                dist.append([j, angle_score(ccenter, tcenter)- d/perimeter]) # 2. scoring: quite complex @@

        print(dist)
        if len(dist) == 0:
            imshow(table, data[img]['cfboxes'][i:i+1], [])
        else:
            dist = np.array(dist)
            j = int(dist[np.argmax(dist[:,1]),0]) # max score
            print(j)
            imshow(table, data[img]['cfboxes'][i:i+1], data[img]['tboxes'][j:j+1])
        




