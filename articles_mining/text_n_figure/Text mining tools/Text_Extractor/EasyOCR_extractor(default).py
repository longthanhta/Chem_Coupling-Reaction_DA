# import library
import easyocr
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PIL import Image
import matplotlib.image as mpimg

def reading():
    # Setting
    reader = easyocr.Reader(['en'])
    #input file
    #reading('/media/gordon/Backup Plus/Documents/Research-Prof.Su/DataMiningProject/EsayOCR/Table_Extractor/A.png')
    img = cv2.imread("/media/gordon/Backup Plus/Documents/Research-Prof.Su/DataMiningProject/EsayOCR/Table_Extractor/A.png")

    # Image Target 
    global result
    result = reader.readtext("/media/gordon/Backup Plus/Documents/Research-Prof.Su/DataMiningProject/EsayOCR/Table_Extractor/A.png")    #list
    # output result
    print(result)


    # Create a Rectangle 
    for i in result:
        box_corr = i[0]
        x1 = int(box_corr[0][0])
        y1 = int(box_corr[0][1])
        x2= int(box_corr[2][0])
        y2 = int(box_corr[2][1])



        cv2.rectangle(img, (x1,y1), (x2,y2), (0, 255, 0), 1)

    cv2.imshow('My Image', img)
    cv2.waitKey(0)
    #cv2.destroyAllWindows()

    return result

 
reading()
