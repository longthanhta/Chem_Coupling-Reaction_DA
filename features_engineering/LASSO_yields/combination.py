import pandas as pd
data_df=pd.read_csv('ft.csv')

header_lst=list(data_df.columns.values)
combination_lst=[]
for header1 in header_lst:
    for header2 in header_lst:
        cbs=header1+'*'+header2#combination string
        cbs_i=header2+'*'+header1#combination string
        if cbs in combination_lst or cbs_i in combination_lst: continue
        else:
            combination_lst.append(cbs)
            combination_lst.append(cbs_i)
        comb_com=[]
        for index in range(len(data_df)):
            comb_com.append(data_df[header1][index]*data_df[header1][index])
        data_df[cbs]=comb_com
data_df.to_csv('ft.csv',index=None)
