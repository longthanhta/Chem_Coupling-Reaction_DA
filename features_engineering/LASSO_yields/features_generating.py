import pandas as pd


def get_en(elmt_check):
    if elmt_check=='H': return	2.2
    if elmt_check=='Li': return	0.98
    if elmt_check=='Be': return	1.57
    if elmt_check=='B': return	2.04
    if elmt_check=='C': return	0
    if elmt_check=='N': return	3.04
    if elmt_check=='O': return	3.44
    if elmt_check=='F': return	3.98
    if elmt_check=='Na': return	0.93
    if elmt_check=='Mg': return	1.31
    if elmt_check=='Al': return	1.61
    if elmt_check=='Si': return	1.9
    if elmt_check=='P': return	2.19
    if elmt_check=='S': return	2.58
    if elmt_check=='Cl': return	3.16
    if elmt_check=='K': return	0.82
    if elmt_check=='Ca': return	1
    if elmt_check=='Sc': return	1.36
    if elmt_check=='Ti': return	1.54
    if elmt_check=='V': return	1.63
    if elmt_check=='Cr': return	1.66
    if elmt_check=='Mn': return	1.55
    if elmt_check=='Fe': return	1.83
    if elmt_check=='Co': return	1.88
    if elmt_check=='Ni': return	1.91
    if elmt_check=='Cu': return	1.9
    if elmt_check=='Zn': return	1.65
    if elmt_check=='Ga': return	1.81
    if elmt_check=='Ge': return	2.01
    if elmt_check=='As': return	2.18
    if elmt_check=='Se': return	2.55
    if elmt_check=='Br': return	2.96
    if elmt_check=='Kr': return	3
    if elmt_check=='Rb': return	0.82
    if elmt_check=='Sr': return	0.95
    if elmt_check=='Y': return	1.22
    if elmt_check=='Zr': return	1.33
    if elmt_check=='Nb': return	1.6
    if elmt_check=='Mo': return	2.16
    if elmt_check=='Tc': return	1.9
    if elmt_check=='Ru': return	2.2
    if elmt_check=='Rh': return	2.28
    if elmt_check=='Pd': return	2.2
    if elmt_check=='Ag': return	1.93
    if elmt_check=='Cd': return	1.69
    if elmt_check=='In': return	1.78
    if elmt_check=='Sn': return	1.96
    if elmt_check=='Sb': return	2.05
    if elmt_check=='Te': return	2.1
    if elmt_check=='I': return	2.66
    if elmt_check=='Xe': return	2.6
    if elmt_check=='Cs': return	0.79
    if elmt_check=='Ba': return	0.89
    if elmt_check=='La': return	1.1
    if elmt_check=='Ce': return	1.12
    if elmt_check=='Pr': return	1.13
    if elmt_check=='Nd': return	1.14
    if elmt_check=='Sm': return	1.17
    if elmt_check=='Gd': return	1.2
    if elmt_check=='Dy': return	1.22
    if elmt_check=='Ho': return	1.23
    if elmt_check=='Er': return	1.24
    if elmt_check=='Tm': return	1.25
    if elmt_check=='Lu': return	1.27
    if elmt_check=='Hf': return	1.3
    if elmt_check=='Ta': return	1.5
    if elmt_check=='W': return	2.36
    if elmt_check=='Re': return	1.9
    if elmt_check=='Os': return	2.2
    if elmt_check=='Ir': return	2.2
    if elmt_check=='Pt': return	2.28
    if elmt_check=='Au': return	2.54
    if elmt_check=='Hg': return	2
    if elmt_check=='Tl': return	1.62
    if elmt_check=='Pb': return	2.33
    if elmt_check=='Bi': return	2.02
    if elmt_check=='Po': return	2
    if elmt_check=='At': return	2.2
    if elmt_check=='Ra': return	0.9
    if elmt_check=='Ac': return	1.1
    if elmt_check=='Th': return	1.3
    if elmt_check=='Pa': return	1.5
    if elmt_check=='U': return	1.38
    if elmt_check=='Np': return	1.36
    if elmt_check=='Pu': return	1.28
    if elmt_check=='Am': return	1.3
    if elmt_check=='Cm': return	1.3
    if elmt_check=='Bk': return	1.3
    if elmt_check=='Cf': return	1.3
    if elmt_check=='Es': return	1.3
    if elmt_check=='Fm': return	1.3
    if elmt_check=='Md': return	1.3
    if elmt_check=='No': return	1.3
def get_ham(sub_smi,pos,hc_df):
    match_row=hc_df[(hc_df['sub']==sub_smi)]
    if pos == 'p':
        result = match_row['p'].values
    if pos == 'o':
        result = match_row['p'].values
    if pos == 'm':
        result = match_row['m'].values
    if len(result)==0:
        return False
    else:
        return result[0]
def get_vol(sub_smi,vol_df):
    match_row=vol_df[(vol_df['sub']==sub_smi)]
    result1 = match_row['r1'].values
    result2 = match_row['r2'].values
    if len(result1) ==0 and len(result2)==0:
        return False, False
    else:
        return result1[0],result2[0]
if __name__ == '__main__':
    data_df=pd.read_excel('CC_CC_data_Ph_Ph_yield_sub.xlsx').fillna('')
    hc_df=pd.read_excel('hammet.xlsx')
    vol_df=pd.read_excel('sub_vol.xlsx')
    ft_dict= {'sub_1_hc_p':[], 'sub_1_r1_p':[],'sub_1_r2_p':[], \
            'sub_1_hc_m':[], 'sub_1_r1_m':[],'sub_1_r2_m':[], \
            'sub_1_hc_o':[], 'sub_1_r1_o':[],'sub_1_r2_o':[], \
            'sub_2_hc_p':[], 'sub_2_r1_p':[],'sub_2_r2_p':[], \
            'sub_2_hc_m':[], 'sub_2_r1_m':[],'sub_2_r2_m':[], \
            'sub_2_hc_o':[], 'sub_2_r1_o':[],'sub_2_r2_o':[], \
            'egc1_en':[],'egc2_en':[]}
    missing_info1_lst=[]
    missing_info2_lst=[]

    for row_i, row in data_df.iterrows():
        print(row_i)


        sub_1_hc_p_sum=0
        sub_2_hc_p_sum=0

        sub_1_hc_m_sum=0
        sub_2_hc_m_sum=0

        sub_1_hc_o_sum=0
        sub_2_hc_o_sum=0


        sub_1_r1_p_sum=0
        sub_1_r2_p_sum=0
        sub_2_r1_p_sum=0
        sub_2_r2_p_sum=0


        sub_1_r1_m_sum=0
        sub_1_r2_m_sum=0
        sub_2_r1_m_sum=0
        sub_2_r2_m_sum=0


        sub_1_r1_o_sum=0
        sub_1_r2_o_sum=0
        sub_2_r1_o_sum=0
        sub_2_r2_o_sum=0


        egc1,egc2=row['EG_class'].split('.')

        sub1=row['sub1']
        sub2=row['sub2']
        pos1=row['pos1']
        pos2=row['pos2']

#       split case
        if '.' in sub1: sub1_lst=sub1.split('.')
        else: sub1_lst=[sub1]
        if '.' in sub2: sub2_lst=sub2.split('.')
        else: sub2_lst=[sub2]
        if '.' in pos1: pos1_lst=pos1.split('.')
        else: pos1_lst=[pos1]
        if '.' in pos2: pos2_lst=pos2.split('.')
        else: pos2_lst=[pos2]

        missing_info=False
        for sub_i in range(len(sub1_lst)):
            sub_s=sub1_lst[sub_i]
            pos=pos1_lst[sub_i]

            if pos == 'm':
                sub_1_hc_m=get_ham(sub_s,pos,hc_df)
                sub_1_r1_m,sub_1_r2_m=get_vol(sub_s,vol_df)
                if sub_1_hc_m:
                    sub_1_hc_m_sum+=sub_1_hc_m
                else: missing_info=True
                if sub_1_r1_m:
                    sub_1_r1_m_sum+=sub_1_r1_m
                if sub_1_r2_m:
                    sub_1_r2_m_sum+=sub_1_r2_m
            if pos == 'p':
                sub_1_hc_p=get_ham(sub_s,pos,hc_df)
                sub_1_r1_p,sub_1_r2_p=get_vol(sub_s,vol_df)
                if sub_1_hc_p:
                    sub_1_hc_p_sum+=sub_1_hc_p
                else: missing_info=True
                if sub_1_r1_p:
                    sub_1_r1_p_sum+=sub_1_r1_p
                if sub_1_r2_p:
                    sub_1_r2_p_sum+=sub_1_r2_p
            if pos == 'o':
                sub_1_hc_o=get_ham(sub_s,pos,hc_df)
                sub_1_r1_o,sub_1_r2_o=get_vol(sub_s,vol_df)
                if sub_1_hc_o:
                    sub_1_hc_o_sum+=sub_1_hc_o*0.075
                else: missing_info=True
                if sub_1_r1_o:
                    sub_1_r1_o_sum+=sub_1_r1_o
                if sub_1_r2_o:
                    sub_1_r2_o_sum+=sub_1_r2_o
        missing_info1_lst.append(missing_info)

        missing_info=False
        for sub_i in range(len(sub2_lst)):
            sub_s=sub2_lst[sub_i]
            pos=pos2_lst[sub_i]

            if pos == 'm':
                sub_2_hc_m=get_ham(sub_s,pos,hc_df)
                sub_2_r1_m,sub_2_r2_m=get_vol(sub_s,vol_df)
                if sub_2_hc_m:
                    sub_2_hc_m_sum+=sub_2_hc_m
                else: missing_info=True
                if sub_2_r1_m:
                    sub_2_r1_m_sum+=sub_2_r1_m
                if sub_2_r2_m:
                    sub_2_r2_m_sum+=sub_2_r2_m

            if pos == 'p':
                sub_2_hc_p=get_ham(sub_s,pos,hc_df)
                sub_2_r1_p,sub_2_r2_p=get_vol(sub_s,vol_df)
                if sub_2_hc_p:
                    sub_2_hc_p_sum+=sub_2_hc_p
                else: missing_info=True
                if sub_2_r1_p:
                    sub_2_r1_p_sum+=sub_2_r1_p
                if sub_2_r2_p:
                    sub_2_r2_p_sum+=sub_2_r2_p

            if pos == 'o':
                sub_2_hc_o=get_ham(sub_s,pos,hc_df)
                sub_2_r1_o,sub_2_r2_o=get_vol(sub_s,vol_df)
                if sub_2_hc_o:
                    sub_2_hc_o_sum+=sub_2_hc_o*0.075
                else: missing_info=True
                if sub_2_r1_o:
                    sub_2_r1_o_sum+=sub_2_r1_o
                if sub_2_r2_o:
                    sub_2_r2_o_sum+=sub_2_r2_o


        missing_info2_lst.append(missing_info)

        egc1_en=get_en(egc1)
        egc2_en=get_en(egc2)

        ft_dict['sub_1_hc_p'].append(sub_1_hc_p_sum)
        ft_dict['sub_1_r1_p'].append(sub_1_r1_p_sum)
        ft_dict['sub_1_r2_p'].append(sub_1_r2_p_sum)

        ft_dict['sub_1_hc_m'].append(sub_1_hc_m_sum)
        ft_dict['sub_1_r1_m'].append(sub_1_r1_m_sum)
        ft_dict['sub_1_r2_m'].append(sub_1_r2_m_sum)

        ft_dict['sub_1_hc_o'].append(sub_1_hc_o_sum)
        ft_dict['sub_1_r1_o'].append(sub_1_r1_o_sum)
        ft_dict['sub_1_r2_o'].append(sub_1_r2_o_sum)

        ft_dict['sub_2_hc_p'].append(sub_2_hc_p_sum)
        ft_dict['sub_2_r1_p'].append(sub_2_r1_p_sum)
        ft_dict['sub_2_r2_p'].append(sub_2_r2_p_sum)

        ft_dict['sub_2_hc_m'].append(sub_2_hc_m_sum)
        ft_dict['sub_2_r1_m'].append(sub_2_r1_m_sum)
        ft_dict['sub_2_r2_m'].append(sub_2_r2_m_sum)

        ft_dict['sub_2_hc_o'].append(sub_2_hc_o_sum)
        ft_dict['sub_2_r1_o'].append(sub_2_r1_o_sum)
        ft_dict['sub_2_r2_o'].append(sub_2_r2_o_sum)

        ft_dict['egc1_en'].append(egc1_en)
        ft_dict['egc2_en'].append(egc2_en)
        #break
data_df['mi1']=missing_info1_lst
data_df['mi2']=missing_info2_lst
data_df.to_excel('CC_CC_data_Ph_Ph_yield_sub.xlsx',index=None)
ft_df=pd.DataFrame(ft_dict)
ft_df.to_csv('ft.csv',index=None)
    #print(data_df)
