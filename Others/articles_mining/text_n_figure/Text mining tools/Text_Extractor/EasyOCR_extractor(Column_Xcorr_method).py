#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  1 18:24:08 2021

@author: gordon
"""
import torch
#Test is the CUDA connected to GPU
print("torch.cuda state:", torch.cuda.is_available())
import easyocr
from numpy import array

''' EasyOCR will give you a list. Here is the information of their construction:
    Example:
    result[0] : [([[11, 3], [75, 3], [75, 21], [11, 21]], 'Table 1', 0.3272569179534912)]
    result[0][0] is coordination of box
    result[0][1] is extracted informations
    result[0][2] is Confidence Interval  
    '''

def reading(file):
    # Setting
    reader = easyocr.Reader(['en']) 
    # Image Target 
    global result
    result = reader.readtext(file)    #list
    # output result
    #print(result)
    return result
    
def CheckingType():
    global all_extracted_text
    global all_coordinate
    all_coordinate = []
    all_extracted_text = []
    for i in result:
        #print(i)
        all_extracted_text.append(i[1])
        all_coordinate.append(i[0])
    #print(all_extracted_text)
    if "Entry" in all_extracted_text:
        print("Article Tpe: Entry")
        TakeAction_Entry_Type()
    else:
        print("Article Tpe: Others")

def TakeAction_Entry_Type():
    #Check Back Column_X_Corr (First Step of Searching for Yield in Entry type article)
    Back_Column_X_Corr = []
    for i in all_coordinate:
        Back_Column_X_Corr.append(i[1][0])
    
    print(Back_Column_X_Corr)
    unique_Back_Column_X_Corr = unique(Back_Column_X_Corr)
    print(unique_Back_Column_X_Corr)
    #Extract those X coordinate if they appera more than once
    X_That_Appeared_MoreThanOnce = []
    for i in unique_Back_Column_X_Corr:
        if countX(Back_Column_X_Corr,i) >1:
            X_That_Appeared_MoreThanOnce.append(i)
            X_That_Appeared_MoreThanOnce.sort()
    print(X_That_Appeared_MoreThanOnce)
            
    for i in X_That_Appeared_MoreThanOnce:
        for X in all_coordinate:
            if X == i[0][1][0]:
                print(i[1])
                
 
        
                    
def countX(lst, x): 
    count = 0
    for ele in lst: 
        if (ele == x): 
            count = count + 1
    return count

def unique(list1): 
  
    # intilize a null list 
    unique_list = [] 
      
    # traverse for all elements 
    for x in list1: 
        # check if exists in unique_list or not 
        if x not in unique_list: 
            unique_list.append(x) 
    # print list 
    #for x in unique_list: 
    #    print (x)
    return  unique_list
def main():
    #input file
    reading('/media/gordon/Backup Plus/Documents/Research-Prof.Su/DataMiningProject/EsayOCR/test05_Table.png')
    CheckingType()


main()




