# import library
import easyocr
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PIL import Image
import matplotlib.image as mpimg

file = "0329test.png"

def box_plotting(file,informationLIST):
    img = cv2.imread(file)
    for i in informationLIST:
        box_corr = i[0]
        x1 = int(box_corr[0][0])
        y1 = int(box_corr[0][1])
        x2= int(box_corr[2][0])
        y2 = int(box_corr[2][1])



        cv2.rectangle(img, (x1,y1), (x2,y2), (0, 255, 0), 1)

    cv2.imshow('My Image', img)
    cv2.waitKey(0)
    cv2.destroyAllWindows()    
    
def reading():
    # Setting
    reader = easyocr.Reader(['en'])
    #input file
    #reading('/media/gordon/Backup Plus/Documents/Research-Prof.Su/DataMiningProject/EsayOCR/Table_Extractor/A.png')
    img = cv2.imread("0329test.png")

    # Image Target 
    global result
    result = reader.readtext("0329test.png")    #list
    # output result
    #print(result)


    # Create a Rectangle 
    for i in result:
        box_corr = i[0]
        x1 = int(box_corr[0][0])
        y1 = int(box_corr[0][1])
        x2= int(box_corr[2][0])
        y2 = int(box_corr[2][1])



        cv2.rectangle(img, (x1,y1), (x2,y2), (0, 255, 0), 1)

    #cv2.imshow('My Image', img)
    #cv2.waitKey(0)
    #cv2.destroyAllWindows()

    return result

 
Information_extracted = reading()
print("--------------Extracted Information----------")
print(Information_extracted)
box_plotting(file,Information_extracted)
print("--------------Restoring the yield related information----------")
Selected_Information=[]
for i in Information_extracted:
    if i[1].find("%") != -1:
        Selected_Information.append(i)
    else:
        pass
print(Selected_Information)
#box_plotting(file,Selected_Information)
print("--------------Detect the nonuseful paragraph information and delete it with charater > 30----------")
para_information = []

for i in Selected_Information:
    if len(i[1]) >= 30:
        print(i[1])
        para_information.append(i)
    else:
        pass    

for i in para_information:
    Selected_Information.remove(i)
print("--------------Finish para information del.---------")



box_plotting(file,Selected_Information)




