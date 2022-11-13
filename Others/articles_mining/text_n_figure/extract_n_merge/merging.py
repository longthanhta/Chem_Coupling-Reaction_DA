import math
import numpy as np

def distance_b(box1, box2):
    dx = max(0, box2[0]-box1[2], box1[0]-box2[2])
    dy = max(0, box2[1]-box1[3], box1[1]-box2[3])
    return math.sqrt(dx*dx + dy**2)
def distance_p(point1, point2):
    dx = point2[0] - point1[0]
    dy = point2[1] - point1[1]
    return math.sqrt(dx*dx + dy**2)

def merge_textboxes(tboxes):
    # return boxes as clusters
    # eg. [[b00,b01,...],[b10,b11,..],...]
    # how close two boxes is to be put into
    # the same cluster? 10% of their total semi perimeter

    clusters = []
    while len(tboxes) > 0:
        # starting with a new one
        print('=== new group ===')
        group = [0]
        left = list(range(1,len(tboxes)))

        # boxes is < 100 pixels from the box above
        # can be in the same group
        candidates = []
        for i in range(1,len(tboxes)):
            if distance_b(tboxes[i], tboxes[0]) < 100:
                candidates.append(i)

        has_thing_added = True
        while has_thing_added:
            has_thing_added = False
            ctemp = candidates.copy()
            for i in candidates:
                box = tboxes[i]

                for j in group:
                    b = tboxes[j]
                    distance = distance_b(box,b)
                    semi_perimeter = box[3]+box[2]-box[1]-box[0] + b[3]+b[2]-b[1]-b[0]
                    relative_distance = distance/semi_perimeter
                    if relative_distance < 0.1:
                        has_thing_added = True
                        group.append(i)
                        left.remove(i)
                        ctemp.remove(i)
                        break
            
            candidates = ctemp.copy()

        clusters.append(tboxes[group])
        tboxes = tboxes[left]

    return clusters


if __name__ == '__main__':
    data = {'catal-4': { 
            'cfboxes' : [[11,185,115,252],[125,180,276,243],[287,173,387,246],[394,155,531,257],[223,635,507,704]], 
            'tboxes' : [[20,258,101,275],[175,260,254,277],[279,256,353,273],[411,259,472,274],[429,680,511,697],[210,589,279,605],[20,258,101,270]]
            },
        'anie-1': { 
            'cfboxes' : [[77,217,149,285],[402,211,536,275],[403,300,536,369],[403,823,541,889]], 
            'tboxes' : [[425,279,518,294],[424,372,522,387],[426,892,531,905]]
            }
        }

    for img in data.keys():
        print(img)
        c = merge_textboxes(np.array(data[img]['tboxes']))
        print(data[img]['tboxes'])
        print(c)

