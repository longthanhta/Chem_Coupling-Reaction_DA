import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from shutil import copy
pd.options.mode.chained_assignment = None #<- avoid pandas's annoucment
#------------------------------------------------------------------------------

def get_uni_list(df,prop):
    unique_element_lst=[]
    for i, row in df.iterrows():
        element_str=str(row[prop])
        if '.' in element_str: # A.B -> [A,B] ; A->[A]
            elements=element_str.split('.')
        else:
            elements=[element_str]

        # Get the unique element lst
        for e in elements:
            if len(e)>0 and e not in unique_element_lst: # None will have string len is 0
                unique_element_lst.append(e)
    return unique_element_lst


def get_color(string):
    if string=='H'  : return '#1F77B4'
    if string=='O'  : return '#FF7F0E'
    if string=='C'  : return '#1B9E77'
    if string=='Mg' : return '#D62728'
    if string=='B'  : return '#9467BD'
    if string=='N'  : return '#8C564B'
    if string=='S'  : return '#E377C2'
    if string=='Al' : return '#7F7F7F'
    if string=='Zn' : return '#BCBD22'
    if string=='Si' : return '#17BECF'
    if string=='X'  : return 'moccasin'
    if string=='In' : return '#C6DBEF'
    if string=='NA' : return 'gray'
    return 'black'




def get_uni_list(df,prop):
    unique_element_lst=[]
    for i, row in df.iterrows():
        element_str=str(row[prop])
        if '.' in element_str: # A.B -> [A,B] ; A->[A]
            elements=element_str.split('.')
        else:
            elements=[element_str]
        # Get the unique element lst
        for e in elements:
            if len(e)>0 and e != 'nan' and e not in unique_element_lst: # None will have string len is 0
                unique_element_lst.append(e)
    return unique_element_lst

def get_uni_pair_list(df,prop,same_scaffold=True):
    unique_pair_lst=['placeholder']
    for i, row in df.iterrows():
        pair_str=str(row[prop])
        if '.' in pair_str: # A.B -> [A,B] ; A->[A]
            A1,A2=pair_str.split('.')
        else:
            A1=pair_str
            A2='NA'
        # Get the unique pair lst
        #for e in pairs:
        pair_check=A1+'.'+A2
        pair_rev_check=A2+'.'+A1
        #print('debug pair,pair_rev',pair_check,pair_rev_check)
        for up in unique_pair_lst:
            if same_scaffold:
                if pair_check not in unique_pair_lst and pair_rev_check not in unique_pair_lst:
                    unique_pair_lst.append(pair_check)
            else:
                if pair_check not in unique_pair_lst:
                    unique_pair_lst.append(pair_check)
    unique_pair_lst.remove('placeholder')
    return unique_pair_lst

def filtering(df_in,prop,ue):
    #print('ue to filtering',ue)
    #print('df before filtering',df)
    check=[]
    for index,row in df_in.iterrows():
        element=row[prop]
        element=str(element) #<- 1,2 will be considered as number
        c=False
        if '.' in element:
        # This for the set with more than one Leaving group/ more than one ligands:
        # Like A.B,A will be counted as 2, A.A,B will be counted as only 1
            element_lst=element.split(".")
            for item in element_lst:
                if item == ue:
                    c=True
        else:
            if element == ue:
                c=True
        check.append(c)
    header_lst=list(df_in.columns.values)
    if 'check' in header_lst:
        df_in.drop(columns='check')
    df_in['check']=check
    df_filtered=df_in[(df_in['check']==True)]
    df_filtered=df_filtered.drop('check', 1)
    #print('df after filtering',df_filtered)
    return df_filtered

def filtering_pair(df_in,prop,up,same_scaffold=True): # LEAVE_A but unique will be a pair
    # up is unique_pair
    check=[]
    for index,row in df_in.iterrows():
        pair=row[prop]
        pair=str(pair) #<- 1,2 will be considered as number
        c=False
        if '.' in pair:
        # This for the set with more than one Leaving group/ more than one ligands:
        # Like A.B,A will be counted as 2, A.A,B will be counted as only 1
            A1,A2=pair.split(".")
            pair_check=A1+'.'+A2
            pair_rev_check=A2+'.'+A1
            if same_scaffold:
                if pair_check == up or pair_rev_check == up:
                    c=True
            else:
                if pair_check == up:
                    c=True
        else:
            if pair == up:
                c=True
        check.append(c)
    header_lst=list(df_in.columns.values)
    if 'check' in header_lst:
        df_in.drop(columns='check')
    df_in['check']=check
    #print('debug df, check for',up, df_in[['LEAVE_A','check']].to_string())
    df_filtered=df_in[(df_in['check']==True)]
    df_filtered=df_filtered.drop('check', 1)
    df_the_rest=df_in[(df_in['check']==False)]
    df_the_rest=df_the_rest.drop('check', 1)
    #print('df after filtering',df_filtered)
    return df_filtered,df_the_rest

def counting(df,unique_element_lst,prop,pair_mode=False,return_the_rest=False,same_scaffold=True):
    if not pair_mode:
        nor_lst=[]
        norp_lst=[]


        #print('unique_element_lst',unique_element_lst)
        for ue in unique_element_lst:
            if ue =='nan' : continue
            # Main one
            #print('ue',ue)
            # ue: unique element
            df_filtered=filtering(df,prop,ue)

            nor_lst.append(len(df_filtered))
        return nor_lst
    else:
        nor_lst=[]
        norp_lst=[]


        #print('unique_element_lst',unique_element_lst)
        for ue in unique_element_lst:
            if ue =='nan' : continue
            # Main one
            #print('ue',ue)
            # ue: unique element
            df_filtered,df_the_rest=filtering_pair(df,prop,ue,same_scaffold)

            nor_lst.append(len(df_filtered))
        if return_the_rest:
            return nor_lst,df_the_rest
        else:
            return nor_lst

def library_generating(cluster_df,EG,EGi_path,reactions_folder,mech_lst):
    cluster_df_EGi=filtering(cluster_df,'LEAVE_A',EG)
    EGisub_unique_lst=get_uni_list(cluster_df_EGi,'LEAVE_A')
    for EGisub in EGisub_unique_lst:
        EGisub_path=EGi_path+'/'+EGisub
        #oks.mkdir(EGisub_path)
        cluster_df_EGisub=filtering(cluster_df_EGi,'LEAVE_A',EGisub)
        prev_REF=''
        for index,row in cluster_df_EGisub.iterrows():
            OID=str(row['OID']) # OID = Our ID, to distinguish with Reaxys ID
            REF=str(row['FILE NAME'])
            if REF == prev_REF: continue #<- each article should only have 1 example
            prev_REF = REF
            scr=reactions_folder+'/'+OID+'.png'
            drt=EGisub_path+'/'+OID+'_'+REF+'.png'
            copy(scr,drt)
# Get mechanism picture
            for picname in mech_lst:
                if REF in picname:
                    #print('found mechanism',REF)
                    copy(mechanism_folder+'/'+picname,EGisub_path+'/mechanism'+'_'+REF+'.png')

#------------------------------------------------------------------------------
reactions_folder='1_reactions'
reaction_data=pd.read_excel('all_data_2020_Jan_1.xlsx',sheet_name='Identified CC')
mechanism_folder='mechanism_picture'
#get lst of mechanism_folder_files
mech_lst=[]
for path, subdirs, files in os.walk(mechanism_folder):
    for index,file in enumerate(files):
        #print(index)
        #print('files',file)
        #if '.png' not in file or '.jpg' not in file: continue
        mech_lst.append(file)
print('mech_lst',mech_lst) #<-DEBUGGING
#reaction_data=pd.read_excel('test.xlsx')

# Get the coordinate and reactions counting for the plot -----------------------



# This part to out unique lst file
'''
uni_lst_EG=get_uni_list(reaction_data,'LEAVE_A')
uni_lst_FB=get_uni_list(reaction_data,'RTFB')
uni_lst_LG=get_uni_list(reaction_data,'CAT_GROUP')
pd.DataFrame(uni_lst_EG).to_csv('uni_lst_EG.csv',index=None,header=['Item'])
pd.DataFrame(uni_lst_FB).to_csv('uni_lst_FB.csv',index=None,header=['Item'])
pd.DataFrame(uni_lst_LG).to_csv('uni_lst_LG.csv',index=None,header=['Item'])
#'''

uni_lst_EG=pd.read_csv('uni_lst_EG.csv')['Item'].values.tolist()
uni_lst_LG=pd.read_csv('uni_lst_LG.csv')['Item'].values.tolist()
uni_lst_FB=pd.read_csv('uni_lst_FB.csv')['Item'].values.tolist()

EG_coord=[]
LG_coord=[]
FB_coord=[]
EG_label=[]
LG_label=[]
FB_label=[]
count_coord=[]

for i,row in reaction_data.iterrows():
        FB_df=str(row['RTFB'])
        FB_uni_ind=uni_lst_FB.index(FB_df)


        LGs_df=str(row['CAT_GROUP'])
        if '.' in LGs_df: LGs=LGs_df.split('.')
        else: LGs=[LGs_df]
        for LG in LGs:
            if LG in uni_lst_LG:
                LG_uni_ind=uni_lst_LG.index(LG)
            else:
               LG=''

            EGs_df=str(row['LEAVE_A'])
            if '.' in EGs_df: EGs=EGs_df.split('.')
            else: EGs=[EGs_df]
            for EG in EGs:
                if EG in uni_lst_EG:
                    EG_uni_ind=uni_lst_EG.index(EG)
                else:
                    EG=''
                FB_coord.append(FB_uni_ind)
                FB_label.append(FB_df)
                LG_coord.append(LG_uni_ind)
                LG_label.append(LG)

clusters_coord=pd.DataFrame({'FB_C':FB_coord,
                     'LG_C':LG_coord,
                     'FB_L':FB_label,
                     'LG_L':LG_label,
                     })
clusters_coord_count=clusters_coord.groupby(['FB_C', 'LG_C','FB_L', 'LG_L']).size().reset_index(name='Count')


# Start plotting --------------------------------------------------------------

cluster_data=clusters_coord_count
unique_lst_all=[]
nor_lst_all=[]
norp_lst_all=[]

fig = plt.figure(figsize=(9, 9))
axs = fig.add_subplot(111)
axs.set_axisbelow(True)



#axs.set_xlim(-0.5, 15.5)
#axs.set_ylim(-1, 5.5)


axs.yaxis.grid(color='gainsboro', linestyle='dashed')
axs.xaxis.grid(color='gainsboro', linestyle='dashed')
axs.set_yticks(cluster_data['FB_C'])
axs.set_yticklabels(cluster_data['FB_L'],ha='right')
axs.set_xticks(cluster_data['LG_C'])
axs.set_xticklabels(cluster_data['LG_L'])

for tick in axs.yaxis.get_major_ticks():
    if 'Csp2-Csp2' in tick.label1.get_text():
        tick.label1.set_color('red') #set the color property
    elif 'Csp3-Csp2' in tick.label1.get_text() or 'Csp2-Csp3' in tick.label1.get_text():
        tick.label1.set_color('darkviolet')
    elif 'Csp3-Csp3' in tick.label1.get_text():
        tick.label1.set_color('navy')
    elif 'Csp3-Csp1' in tick.label1.get_text():
        tick.label1.set_color('blue')
    elif 'Csp2-Csp1' in tick.label1.get_text():
        tick.label1.set_color('darkgreen')
    elif 'Csp1-Csp1' in tick.label1.get_text():
        tick.label1.set_color('black')
# Make folder to store library of picture
try:
    oks.mkdir('2_library')
except:
    print('folder "library" already exists')

max_total_count=500
for index, row in cluster_data.iterrows():
    FB=row['FB_L']
    LG=row['LG_L']
    FBC=row['FB_C']
    LGC=row['LG_C']
    print(FB,LG)
# GO THROUGH EACH CLUSTER#-----------------------------------------------------------
    # Make folder for each
    FBns=FB.replace(' ','')
    LGns=LG.replace(' ','')
    cluster_path='2_library/'+FBns+'_'+LGns
    cluster_path=cluster_path.replace(' ','')
    #oks.mkdir(cluster_path)
    cluster_df=reaction_data[(reaction_data['RTFB']==FB) \
    &(reaction_data['CAT_GROUP']==LG)]
    print('number of reaction in cluster',len(cluster_df))
    total_count=len(cluster_df)
    if total_count==0: continue
    total_ref=len(set(cluster_df['FILE NAME'].values.tolist()))
    try:
        max_yield=int(cluster_df['Yield'].max())
    except:
        max_yield=0
    #min_yield=cluster_df['YIELD'].min()
    # These plot no need splitting (because scaffolds are the same)

# FIRST CASE: Spliting pie chart
    if FB != 'Alkyl-Alkyl Cross Coupling Csp3-Csp3' and \
    FB != 'Aryl-Aryl Cross Coupling Csp2-Csp2' and \
    FB != 'Alkyne-Alkyne Cross Coupling Csp1=Csp1' \
    and FB != 'Viny-Vinyl Cross Coupling Csp2-Csp2':

    # Get the unique pair and the number of reaction for each pair
        unique_pair_lst=get_uni_pair_list(cluster_df,'LEAVE_A',same_scaffold=False)
        nor_pair_lst=counting(cluster_df,unique_pair_lst,'LEAVE_A',pair_mode=True,same_scaffold=False)
        print('unique_pair',unique_pair_lst)
        print('nor_pair_lst',nor_pair_lst)

    # Rank the pair based on number of reaction
        pair_nor_df=pd.DataFrame({'pair':unique_pair_lst,'nor':nor_pair_lst})
        pair_nor_df=pair_nor_df.sort_values(by='nor', ascending=False) #<-sort by nor

    # Get only the top pair
        top_pair=pair_nor_df['pair'].values.tolist()[0]
        top_pair_nor=pair_nor_df['nor'].values.tolist()[0]
        print('top_pair',top_pair,'nor',top_pair_nor)

    # Split data to 2 part, the part contains the top pair and the rest:
        top_pair_df,the_rest_df=filtering_pair(cluster_df,'LEAVE_A',top_pair,same_scaffold=False)
        unique_lst1=get_uni_list(the_rest_df,'LEAVE_A1')
        unique_lst2=get_uni_list(the_rest_df,'LEAVE_A2')
        nor_lst1=counting(the_rest_df,unique_lst1,'LEAVE_A1')
        nor_lst2=counting(the_rest_df,unique_lst2,'LEAVE_A2')
    # Get number of elememnt
        item_nor_df1=pd.DataFrame({'item':unique_lst1,'nor':nor_lst1})
        item_nor_df1=item_nor_df1.sort_values(by='nor', ascending=True) #<-sort by nor
        unique_lst1=item_nor_df1['item'].values.tolist()
        nor_lst1=item_nor_df1['nor'].values.tolist()
    #
        item_nor_df2=pd.DataFrame({'item':unique_lst2,'nor':nor_lst2})
        item_nor_df2=item_nor_df2.sort_values(by='nor', ascending=True) #<-sort by nor
        unique_lst2=item_nor_df2['item'].values.tolist()
        nor_lst2=item_nor_df2['nor'].values.tolist()

        # This sum_all here is the sum of all leaving group, each reaction
        # will  have 2 leaving group, so the number of reaction for each
        # cluster must plus the top_pair_nor
        sum_all1=sum(nor_lst1)+top_pair_nor
        sum_all2=sum(nor_lst2)+top_pair_nor
        sum_all=sum_all1+sum_all2
        # top_pair_nor1 and 2 are euqal but their ratio is different.
        top_pair_norp1=top_pair_nor/sum_all1
        top_pair_norp2=top_pair_nor/sum_all2
        norp_lst1=[(x/sum_all1)/2 for x in nor_lst1]
        norp_lst2=[(x/sum_all2)/2 for x in nor_lst2]
        print('unique_lst1',unique_lst1)
        print('unique_lst2',unique_lst2)
        print('norp_lst1',norp_lst1)
        print('norp_lst2',norp_lst2)
        print('sum of reaction for each side',sum_all1,sum_all2)


        # Add number of reaction and yield
        axs.annotate(str(total_count)+'  ('+str(total_ref)+')',
                    (LGC, FBC),
                    textcoords="offset points", # how to position the text
                    xytext=(0,-20), # distance from text to points (x,y)
                    size=8,
                    #weight='bold',
                    ha='center')
        if max_yield>0:
            axs.annotate(str(max_yield)+' %',
                        (LGC, FBC),
                        textcoords="offset points", # how to position the text
                        xytext=(0,-26), # distance from text to points (x,y)
                        size=7,
                        style='italic',
                        ha='center')
    # # create 2 pie for top_pair
        print('top_pair_norp1','top_pair_norp2',top_pair_norp1,top_pair_norp2)
        TP_A1,TP_A2=top_pair.split('.')
        x1 = np.sin(2 * np.pi * np.linspace(0, top_pair_norp2))
        y1 = np.cos(2 * np.pi * np.linspace(0, top_pair_norp2))
        xy1 = np.row_stack([[0, 0], np.column_stack([x1, y1])])
        axs.scatter(LGC+0.05,FBC+0.1,s=(50+(sum_all/max_total_count)*(500)),\
        c=get_color(TP_A1),marker=xy1,linewidths=0.5)
        x1 = np.sin(2 * np.pi * np.linspace(1-top_pair_norp1, 1))
        y1 = np.cos(2 * np.pi * np.linspace(1-top_pair_norp1, 1))
        xy1 = np.row_stack([[0, 0], np.column_stack([x1, y1])])
        axs.scatter(LGC-0.05,FBC+0.1,s=(50+(sum_all/max_total_count)*(500)),\
        c=get_color(TP_A2),marker=xy1,linewidths=0.5)


# SECOND CASE: Single pie chart-------------------------------------------------
    else: #<= Csp3-Csp3,Csp2-Csp2,Csp1-Csp1,...

    # Get the unique pair and the number of reaction for each pair
        unique_pair_lst=get_uni_pair_list(cluster_df,'LEAVE_A')
        nor_pair_lst=counting(cluster_df,unique_pair_lst,'LEAVE_A',pair_mode=True)
        print('unique_pair',unique_pair_lst)
        print('nor_pair_lst',nor_pair_lst)

    # Rank the pair based on number of reaction
        pair_nor_df=pd.DataFrame({'pair':unique_pair_lst,'nor':nor_pair_lst})
        pair_nor_df=pair_nor_df.sort_values(by='nor', ascending=False) #<-sort by nor

    # Get only the top pair
        top_pair=pair_nor_df['pair'].values.tolist()[0]
        top_pair_nor=pair_nor_df['nor'].values.tolist()[0]
        print('top_pair',top_pair,'nor',top_pair_nor)

    # Split data to 2 part, the part contains the top pair and the rest:
        top_pair_df,the_rest_df=filtering_pair(cluster_df,'LEAVE_A',top_pair)

    # Now the same procedure for element instead of pair
        # Some different:
        # For the pair, number of pair is number of reactions
        # However, for the element, number of pair is 2 time number of reaction
        unique_lst=get_uni_list(the_rest_df,'LEAVE_A')
        nor_lst=counting(the_rest_df,unique_lst,'LEAVE_A')
        print('unique_lst',unique_lst)
        print('nor_lst',nor_lst)

        elmt_nor_df=pd.DataFrame({'elmt':unique_lst,'nor':nor_lst})
        elmt_nor_df=elmt_nor_df.sort_values(by='nor', ascending=True) #<-sort by nor
        unique_lst=elmt_nor_df['elmt'].values.tolist()
        nor_lst=elmt_nor_df['nor'].values.tolist()
        sum_all=sum(nor_lst)+top_pair_nor*2
        # This sum_all here is the sum of all leaving group, each reaction
        # will  have 2 leaving group, so the number of reaction has topair
        # must be time 2
        norp_lst=[(x/sum_all) for x in nor_lst]
        top_pair_norp=top_pair_nor/sum_all
        #norp in here is number of reaction (percentage)
    # Library generating part <- for making slide
        # EG_path=cluster_path+'/EG/'+EG
        # #oks.mkdir(EG_path)
        # library_generating(cluster_df,EG,EG_path,reactions_folder,mech_lst)
        #oks.mkdir(cluster_path+'/EG')

    # create 2 pie for top_pair
        TP_A1,TP_A2=top_pair.split('.')
        x1 = np.sin(2 * np.pi * np.linspace(0, top_pair_norp))
        y1 = np.cos(2 * np.pi * np.linspace(0, top_pair_norp))
        xy1 = np.row_stack([[0, 0], np.column_stack([x1, y1])])
        axs.scatter(LGC+0.05,FBC+0.05,s=(50+(sum_all/max_total_count)*(500)),\
        c=get_color(TP_A1),marker=xy1,linewidths=0.5)
        prev=top_pair_norp
        x1 = np.sin(2 * np.pi * np.linspace(1-top_pair_norp, 1))
        y1 = np.cos(2 * np.pi * np.linspace(1-top_pair_norp, 1))
        xy1 = np.row_stack([[0, 0], np.column_stack([x1, y1])])
        axs.scatter(LGC-0.05,FBC+0.05,s=(50+(sum_all/max_total_count)*(500)),\
        c=get_color(TP_A2),marker=xy1,linewidths=0.5)
        prev=top_pair_norp

fig.tight_layout()
fig.savefig('2_plot_top_bar_only.png', dpi=600)
