import pandas as pd



#INPUT_FILES_______________________________________________________________________________________________________


input_excel_path = '__data_generated/result_16_3_2021_IsAr6_SubINFO_HetINFO.xlsx'

output_excel_path = '__data_generated/result_16_3_2021_Ar6nonHet-Ar6nonHet_SubINFO_HetINFO.xlsx'

#FUNCTIONS_______________________________________________________________________________________________________


#MAIN_______________________________________________________________________________________________________

if __name__ == '__main__':    
    input_df = pd.read_excel(input_excel_path)

    input_df_filtered = input_df.loc[   # both rct1 and rct2 are Ar6
                                        (input_df['IsAr6_rct1'] == 1) & \
                                        (input_df['IsAr6_rct2'] == 1) & \
                                        # no data in Hetero part
                                        pd.isna(input_df['Het_info_rct1']) & \
                                        pd.isna(input_df['Het_info_rct2'])
                                        ]

    input_df_filtered.to_excel(output_excel_path,index=None)

