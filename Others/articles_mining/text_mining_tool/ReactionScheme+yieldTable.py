

#####################################################################
### This tool belong to Haibin Su's group                       #####
### Please ask for permission for any distribution of this tool.#####
#####################################################################



import numpy as np
import easyocr
import cv2
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from PIL import Image
import matplotlib.image as mpimg
import pandas as pd
import subprocess
import os
import shutil
from shutil import copyfile

FileInPut1_ReactionScheme = "Article1_P1.png"
FileInPut2_StructureOrientedTable = "Article1_T1.png"
ReactionSchemeSMILESandCOR = []
ReactionTableSMILESandCOR =[]

def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)

createFolder('ReactionScheme')
createFolder('StructureOriendedTable')
#createFolder('DataProcessingRecord')


createFolder('./ReactionScheme/Structure/')
createFolder('./ReactionScheme/Structure/Picture/')
createFolder('./ReactionScheme/Structure/SMILES/')
createFolder('./ReactionScheme/Reaction Condition/')

createFolder('./StructureOriendedTable/Structure')
createFolder('./StructureOriendedTable/Structure/Picture')
createFolder('./StructureOriendedTable/Structure/SMILES')
createFolder('./StructureOriendedTable/Yield')



def OSRAscheme(file,scalling_factor):
    print("OSRA processing")
    # Define command as string and then split() into list format
    cmd = 'osra -c '+file
    # Use shell to execute the command
    sp = subprocess.Popen(cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True)
    
    # Separate the output and error.
    # This is similar to Tuple where we store two values to two different variables
    out,err=sp.communicate()
    
    # Store the return code in rc variable
    rc=sp.wait()
    
    print('Return Code:',rc,'\n')
    print('output is: \n', out)
    print('error is: \n', err)  
    
    print("--------------Finished First SMILES and coordinate Extraction------------------")
    text_file = open(file+"SMILES&corr_record.txt", "w")
    text_file.write(out)
    text_file.close()
    text_file = open(file+"SMILES&corr_record.txt", "r")
    Lines = text_file.readlines()
    print("Number of SMILES detected:"+ str(len(Lines)))
    #print(Lines)
    print("--------------Restoring coordinate record----------------------------------------")
    box=[]
    for x in Lines:
        empty_position=x.find(" ")
        #print(x[empty_position+1:])
        box.append(x[empty_position+1:])
    print(box)
    XY1_corr = []
    XY2_corr = []
    for i in box:
        splitposition=i.find("-")
        XY1_corr.append(i[:splitposition])
        XY2_corr.append(i[splitposition+1:-1])   
    print(XY1_corr)
    print(XY2_corr)
    X1_corr = []
    X2_corr = []
    Y1_corr = []
    Y2_corr = []
    for i in XY1_corr:
        xposition=i.find("x")
        X1_corr.append(i[:xposition])
        Y1_corr.append(i[xposition+1:])   
    for i in XY2_corr:
        xposition=i.find("x")
        X2_corr.append(i[:xposition])
        Y2_corr.append(i[xposition+1:])   
    print(X1_corr)
    print(X2_corr)
    print(Y1_corr)
    print(Y2_corr)
    Scalling_Factor= scalling_factor
    print("--------------Rescalling coordinate by " + str(Scalling_Factor)+"----------------------------------------")
    X1_ncorr = []
    X2_ncorr = []
    Y1_ncorr = []
    Y2_ncorr = []
    for i in X1_corr:
        i = int(int(i)-int(i)*(Scalling_Factor))
        X1_ncorr.append(i)
    for i in Y1_corr:
        i = int(int(i)-int(i)*(Scalling_Factor))
        Y1_ncorr.append(i)
    for i in X2_corr:
        i = int(int(i)+int(i)*(Scalling_Factor))
        X2_ncorr.append(i)
    for i in Y2_corr:
        i = int(int(i)+int(i)*(Scalling_Factor))
        Y2_ncorr.append(i)
    print(X1_ncorr)
    print(X2_ncorr)
    print(Y1_ncorr)
    print(Y2_ncorr)
    print("--------------Finished Rescalling coordinate by " + str(Scalling_Factor)+"----------------------------------------")
    print("--------------Croping individual Structure----------------------------------------")
    img = cv2.imread(file)
    COUNT=0
    for i in X1_ncorr:
        x1=X1_ncorr[COUNT]
        x2=X2_ncorr[COUNT]
        y1=Y1_ncorr[COUNT]
        y2=Y2_ncorr[COUNT]
        Y =y1
        h = y2-y1
        X=x1
        w=x2-x1
        crop_img = img[Y:Y+h, X:X+w]
        cv2.imwrite('crop_img'+str(COUNT)+".jpg", crop_img)
        COUNT+=1
    print("--------------Finished Croping individual Structure----------------------------------------")
    print("--------------Individual SMILES Extraction----------------------------------------")
    COUNT=0
    global SMILES_LIST_SCHEME
    SMILES_LIST_SCHEME = []
    
    for i in X1_ncorr:
        individual_picture = "crop_img"+str(COUNT)+".jpg"
        cmd = 'osra '+individual_picture
        # Use shell to execute the command
        sp = subprocess.Popen(cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True)
    
    # Separate the output and error.
    # This is similar to Tuple where we store two values to two different variables
        out,err=sp.communicate()
    
    # Store the return code in rc variable
        rc=sp.wait()
    
    
    #print('Return Code:',rc,'\n')
        print(out)
        SMILES_LIST_SCHEME.append(out[:-1])
        text_file2 = open(file+"SMILES_second_extraction.txt", "a")
        corrlist = [X1_ncorr[COUNT],Y1_ncorr[COUNT],X2_ncorr[COUNT],Y2_ncorr[COUNT]]
        
        corrlist = str(corrlist)
        
        ReactionSchemeSMILESandCOR.append((out[:-1],corrlist))
        
        
        text_file2.write(corrlist)
        text_file2.write(out)
        text_file2.close()
        try:
            os.rename("crop_img"+str(COUNT)+'.jpg',out[:-1]+'.jpg')
            copyfile(out[:-1]+'.jpg', './ReactionScheme/Structure/Picture/'+out[:-1]+'.jpg')
        except:
            copyfile("crop_img"+str(COUNT)+'.jpg', './ReactionScheme/Structure/Picture/'+"crop_img"+str(COUNT)+'.jpg')
        """WILL BE FURTHER IMPROVED for different structures in 1 crop image"""
    #print('error is: \n', err) 
        COUNT+=1
    COUNT=0

    print("OSRA processing fininsed")
def RemoveWHITE(file): 
    img = Image.open(file) 
    img = img.convert("RGBA") 
  
    datas = img.getdata() 
  
    newData = [] 
  
    for items in datas: 
        if items[0] == 255 and items[1] == 255 and items[2] == 255: 
            newData.append((255, 255, 255, 0)) 
        else: 
            newData.append(items) 
  
    img.putdata(newData) 
    img.save(file+"Removed_White.png", "PNG") 
    print("White Removal: Successful") 
def  maintain_aspect_ratio_resize(image,  width=None,  height=None,  inter=cv2.INTER_AREA):
     dim =  None
     (h, w)  = image.shape[:2]
     if width is  None  and height is  None:
        return image

     if width is  None:
        r = height /  float(h)
        dim =  (int(w * r), height)
     else:
          r = width /  float(w)
          dim =  (width,  int(h * r))

     return cv2.resize(image, dim,  interpolation=inter)
def RemoveStructure(TEMPLATE,previous_image):
    template = cv2.imread(TEMPLATE)
    template = cv2.cvtColor(template, cv2.COLOR_BGR2GRAY)
    template = cv2.Canny(template,  50,  200)
    (tH, tW)  = template.shape[:2]
    cv2.imshow("template", template)
    
    original_image = cv2.imread(previous_image)
    final = original_image.copy()
    gray = cv2.cvtColor(original_image, cv2.COLOR_BGR2GRAY)
    found =  None
    
    for scale in np.linspace(0.2,  1.0,  20)[::-1]:
         resized = maintain_aspect_ratio_resize(gray,  width=int(gray.shape[1]  * scale))
         r = gray.shape[1]  /  float(resized.shape[1])
    
         if resized.shape[0]  < tH or resized.shape[1]  < tW:
            break
         canny = cv2.Canny(resized,  50,  200)
         detected = cv2.matchTemplate(canny, template, cv2.TM_CCOEFF)
         (_, max_val, _, max_loc)  = cv2.minMaxLoc(detected)
    
         if found is  None  or max_val > found[0]:
            found =  (max_val, max_loc, r)
    
    (_, max_loc, r)  = found
    (start_x, start_y)  =  (int(max_loc[0]  * r),  int(max_loc[1]  * r))
    (end_x, end_y)  =  (int((max_loc[0]  + tW)  * r),  int((max_loc[1]  + tH)  * r))
    
    
    cv2.rectangle(original_image,  (start_x, start_y),  (end_x, end_y),  (0,255,0),  2)
    cv2.imshow('detected', original_image)
    
    cv2.rectangle(final,  (start_x, start_y),  (end_x, end_y),  (255,255,255),  -1)
    cv2.imwrite('ReactionSchemeTextRemaing.png', final)
    #cv2.waitKey(0)
    cv2.destroyAllWindows()
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
def reading(file):
    # Setting
    reader = easyocr.Reader(['en'])
    #input file
    #reading('/media/gordon/Backup Plus/Documents/Research-Prof.Su/DataMiningProject/EsayOCR/Table_Extractor/A.png')
    img = cv2.imread(file)

    # Image Target 
    global result
    result = reader.readtext(file)    #list
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
def OSRAtable(file,scalling_factor):
    print("OSRA processing")
    # Define command as string and then split() into list format
    cmd = 'osra -c '+file
    # Use shell to execute the command
    sp = subprocess.Popen(cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True)
    
    # Separate the output and error.
    # This is similar to Tuple where we store two values to two different variables
    out,err=sp.communicate()
    
    # Store the return code in rc variable
    rc=sp.wait()
    
    print('Return Code:',rc,'\n')
    print('output is: \n', out)
    print('error is: \n', err)  
    
    print("--------------Finished First SMILES and coordinate Extraction------------------")
    text_file = open(file+"SMILES&corr_record.txt", "w")
    text_file.write(out)
    text_file.close()
    text_file = open(file+"SMILES&corr_record.txt", "r")
    Lines = text_file.readlines()
    print("Number of SMILES detected:"+ str(len(Lines)))
    #print(Lines)
    print("--------------Restoring coordinate record----------------------------------------")
    box=[]
    for x in Lines:
        empty_position=x.find(" ")
        #print(x[empty_position+1:])
        box.append(x[empty_position+1:])
    print(box)
    XY1_corr = []
    XY2_corr = []
    for i in box:
        splitposition=i.find("-")
        XY1_corr.append(i[:splitposition])
        XY2_corr.append(i[splitposition+1:-1])   
    print(XY1_corr)
    print(XY2_corr)
    X1_corr = []
    X2_corr = []
    Y1_corr = []
    Y2_corr = []
    for i in XY1_corr:
        xposition=i.find("x")
        X1_corr.append(i[:xposition])
        Y1_corr.append(i[xposition+1:])   
    for i in XY2_corr:
        xposition=i.find("x")
        X2_corr.append(i[:xposition])
        Y2_corr.append(i[xposition+1:])   
    print(X1_corr)
    print(X2_corr)
    print(Y1_corr)
    print(Y2_corr)
    Scalling_Factor= scalling_factor
    print("--------------Rescalling coordinate by " + str(Scalling_Factor)+"----------------------------------------")
    X1_ncorr = []
    X2_ncorr = []
    Y1_ncorr = []
    Y2_ncorr = []
    for i in X1_corr:
        i = int(int(i)-int(i)*(Scalling_Factor))
        X1_ncorr.append(i)
    for i in Y1_corr:
        i = int(int(i)-int(i)*(Scalling_Factor))
        Y1_ncorr.append(i)
    for i in X2_corr:
        i = int(int(i)+int(i)*(Scalling_Factor))
        X2_ncorr.append(i)
    for i in Y2_corr:
        i = int(int(i)+int(i)*(Scalling_Factor))
        Y2_ncorr.append(i)
    print(X1_ncorr)
    print(X2_ncorr)
    print(Y1_ncorr)
    print(Y2_ncorr)
    print("--------------Finished Rescalling coordinate by " + str(Scalling_Factor)+"----------------------------------------")
    print("--------------Croping individual Structure----------------------------------------")
    img = cv2.imread(file)
    COUNT=0
    for i in X1_ncorr:
        x1=X1_ncorr[COUNT]
        x2=X2_ncorr[COUNT]
        y1=Y1_ncorr[COUNT]
        y2=Y2_ncorr[COUNT]
        Y =y1
        h = y2-y1
        X=x1
        w=x2-x1
        crop_img = img[Y:Y+h, X:X+w]
        cv2.imwrite('crop_img'+str(COUNT)+".jpg", crop_img)
        COUNT+=1
    print("--------------Finished Croping individual Structure----------------------------------------")
    print("--------------Individual SMILES Extraction----------------------------------------")
    COUNT=0
    global SMILES_LIST_TABLE
    SMILES_LIST_TABLE = []  
    for i in X1_ncorr:
        individual_picture = "crop_img"+str(COUNT)+".jpg"
        cmd = 'osra '+individual_picture
        # Use shell to execute the command
        sp = subprocess.Popen(cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True)
    
    # Separate the output and error.
    # This is similar to Tuple where we store two values to two different variables
        out,err=sp.communicate()
    
    # Store the return code in rc variable
        rc=sp.wait()
    

    #print('Return Code:',rc,'\n')
        print(out)
        SMILES_LIST_TABLE.append(out[:-1])
        text_file2 = open(file+"SMILES_second_extraction.txt", "a")
        corrlist = [X1_ncorr[COUNT],Y1_ncorr[COUNT],X2_ncorr[COUNT],Y2_ncorr[COUNT]]
        
        corrlist = str(corrlist)
        
        ReactionTableSMILESandCOR.append((out[:-1],corrlist))
        
        
        text_file2.write(corrlist)
        text_file2.write(out)
        text_file2.close()
        try:
            os.rename("crop_img"+str(COUNT)+'.jpg',out[:-1]+'.jpg')
            copyfile(out[:-1]+'.jpg', './StructureOriendedTable/Structure/Picture/'+out[:-1]+'.jpg')
        except:
            copyfile("crop_img"+str(COUNT)+'.jpg', './StructureOriendedTable/Structure/Picture/'+"crop_img"+str(COUNT)+'.jpg')
        """WILL BE FURTHER IMPROVED for different structures in 1 crop image"""
    #print('error is: \n', err) 
        COUNT+=1
    COUNT=0
    
    
"""Reaction Scheme Extraction"""
OSRAscheme(FileInPut1_ReactionScheme,0.03)
RemoveWHITE(FileInPut1_ReactionScheme)
print("SMILES and corresponding coordinate extracted in REACTION SCHEME:")
print(ReactionSchemeSMILESandCOR)
x = FileInPut1_ReactionScheme
for i in ReactionSchemeSMILESandCOR:
    try:
        RemoveStructure(i[0]+".jpg",x)
        x = 'ReactionSchemeTextRemaing.png'
    except:
        print("Error: no this SMILES named picture    (NOT YET FIXED) skipped replacing this structure")
ReactionScheme_Information_extracted = reading('ReactionSchemeTextRemaing.png')
print("--------------Reaction Scheme: Extracted Information----------")
ReactionSchemeText = []
with open('ReactionScheme_Information.txt', 'w') as f:
    for i in ReactionScheme_Information_extracted :
        print(i[1])
        f.write(i[1])
        f.write('\n')
copyfile('ReactionScheme_Information.txt', 'ReactionScheme/Reaction Condition/ReactionScheme_Information.txt')

"""Structure-Oriented Table Extraction"""
OSRAtable(FileInPut2_StructureOrientedTable,0.03)
RemoveWHITE(FileInPut2_StructureOrientedTable)
print("SMILES and corresponding coordinate extracted in REACTION SCHEME:")
print(ReactionTableSMILESandCOR)

Information_extracted = reading(FileInPut2_StructureOrientedTable)
print("--------------Extracted Information----------")
print(Information_extracted)
#box_plotting(file,Information_extracted)
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
#box_plotting(file,Selected_Information)
print("--------------Finding the nearest boxes and combine the boxes.---------")
Information_extracted_with_centering= [] 
for i in Information_extracted:
    centerX = (i[0][0][0]+i[0][1][0])/2
    centerY = (i[0][0][1]+i[0][2][1])/2
    Information_extracted_with_centering.append([i[1],centerX,centerY,i[0]])
    
Selected_Information_with_centering = []
for i in Selected_Information:
    centerX = (i[0][0][0]+i[0][1][0])/2
    centerY = (i[0][0][1]+i[0][2][1])/2
    Selected_Information_with_centering.append([i[1],centerX,centerY,i[0]])
"Checking on the same horiziontal line or not AND nearest boxes"
percentage_Y=0.01
nearest_X_value = 15
for i in Selected_Information_with_centering:
    for y in Information_extracted_with_centering:
        if y == i:
            #print(i)
            pass
        elif abs((y[2]-i[2])/i[2])<=percentage_Y and (abs(i[3][0][0]-y[3][1][0]) <= nearest_X_value or abs(i[3][1][0]-y[3][0][0]) <= nearest_X_value):            
            combined_text = y[0] +" "+ i[0]
            #Re-scale the box
            allX1 = [i[3][0][0],i[3][2][0],y[3][0][0],y[3][2][0]]
            allX2 = [i[3][1][0],i[3][3][0],y[3][1][0],y[3][3][0]]
            allY1 = [i[3][0][1],i[3][2][1],y[3][0][1],y[3][2][1]]
            allY2 = [i[3][1][1],i[3][3][1],y[3][1][1],y[3][3][1]]
            allX1.sort()
            allX2.sort()
            allY1.sort()
            allY2.sort()          
            RescaledBox = [allX1[0],allY1[0]],[allX2[-1],allY1[0]],[allX1[0],allY2[-1]],[allX2[-1],allY2[-1]]
            for z in Selected_Information:
                if i[0] == z[1]:
                    ZZ=list(z)
                    Selected_Information.remove(z)
                    ZZ[1] = combined_text
                    ZZ[0] = RescaledBox
                    ZZ=tuple(ZZ)
                    print(ZZ)
                    Selected_Information.append(ZZ)
                else:
                    pass                  
           
        else:
            pass
"""CAUTION: The following plotting have some problem, BUT THE DATA SHOULD BE FINE (NEED RE-CHECK LATER)"""
#box_plotting(file,Selected_Information)
print("--------------Finish finding nearest horizontal boxes---------")
print("--------------Finish Combining horizontal boxes---------")
print("--------------Finish Updateing boxes information (Text and corrdinate)---------")
print("--------------Finding sub grouping boxes---------")


coordinateLIST=[]
TEXTlist=[]
CIlist=[]
#SMILES=[]
#IMAGE=[]
Yieldlist=[]
for i in Selected_Information:
     coordinateLIST.append(i[0])
     TEXTlist.append(i[1])
     CIlist.append(i[2])
     percentdire = i[1].find('%')
     Yieldlist.append((i[1][percentdire-2:percentdire+1]))
     
Selected_Information_with_centeringX = []
COUNT = 0
Selected_Information_with_centeringY = []
for i in coordinateLIST:
    centerX = (int(i[0][0])+int(i[1][0]))/2
    centerY = (int(i[0][1])+int(i[2][1]))/2
    Selected_Information_with_centeringX.append(centerX)
    Selected_Information_with_centeringY.append(centerY)
    
print(Selected_Information_with_centering)
d = {'coordinate': coordinateLIST, 'text': TEXTlist,'Confident Interval': CIlist,'yield':Yieldlist,'centerX':Selected_Information_with_centeringX,'centerY':Selected_Information_with_centeringY}
df = pd.DataFrame(data=d)
df = df.sort_values(by=['text'])
df.to_csv (FileInPut2_StructureOrientedTable+'yield_extraction.csv', index = False, header=True)
print("csv exported")
copyfile(FileInPut2_StructureOrientedTable+'yield_extraction.csv', './StructureOriendedTable/Yield/'+FileInPut2_StructureOrientedTable+'yield_extraction.csv')
""""Restore the SMILES list"""
print("SMILES_LIST_SCHEME:")
print(SMILES_LIST_SCHEME)
dd= {'SMILES': SMILES_LIST_SCHEME}
ddf = pd.DataFrame(data=dd)
ddf.to_csv (FileInPut1_ReactionScheme+'SMILES.csv', index = False, header=True)
copyfile(FileInPut1_ReactionScheme+'SMILES.csv', './ReactionScheme/Structure/SMILES/'+FileInPut1_ReactionScheme+'SMILES.csv')


print("SMILES_LIST_TABLE:")
print(SMILES_LIST_TABLE)
ddd= {'SMILES': SMILES_LIST_TABLE}
dddf = pd.DataFrame(data=ddd)
dddf.to_csv (FileInPut2_StructureOrientedTable+'SMILES.csv', index = False, header=True)
copyfile(FileInPut2_StructureOrientedTable+'SMILES.csv', './StructureOriendedTable/Structure/SMILES/'+FileInPut2_StructureOrientedTable+'SMILES.csv')

'''remove files'''
arr = os.listdir()
arr.remove("StructureOriendedTable")
arr.remove("ReactionScheme")
arr.remove("ReactionScheme+yieldTable.py")
arr.remove(FileInPut1_ReactionScheme)
arr.remove(FileInPut2_StructureOrientedTable)

for i in arr:
    os.remove(i) 

print("FINISHED")
