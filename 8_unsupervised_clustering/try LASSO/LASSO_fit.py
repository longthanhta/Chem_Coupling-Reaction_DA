import pandas as pd
import plotly.express as px
from sklearn import linear_model
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import GridSearchCV
from sklearn.metrics import r2_score

if __name__ == '__main__':
    df_data = pd.read_excel('data_w_label_2D_l5000_p30_r200_i1000.xlsx',sheet_name='4 Vs 6')
    #df_features=df_data[['E_BDE','Nu_BDE','S_1_CNMR','S_2_CNMR','LIG_VB','PA']]
    df_features=df_data[['S_2_CNMR','Nu_BDE','S_1_CNMR','PA','S_BDE','LIG_VB']]

    #print('df_features',df_features)
    x=df_features

    target=df_data['E_BDE']
    #scaler = StandardScaler()
    #x=scaler.fit_transform(x)
    x_train, x_test, y_train, y_test= train_test_split(x, target, test_size=1)

    #x = np.asarray(x)
    #param_grid = {"alpha": np.logspace(-5, 5, 20)}
    #estimator = GridSearchCV(linear_model.Lasso(), cv=3,param_grid=param_grid)
    #estimator.fit(x,target)
    #best_alpha=estimator.best_estimator_.alpha
    reg = linear_model.Lasso(alpha=0.001)
    reg.fit(x_train,y_train)
    print(reg.coef_)
    ypred = reg.predict(x_test)
    ypred_train = reg.predict(x_train)

# Caculate Mean Absolute Error
    test_er=mean_squared_error(y_test, ypred)
    train_er=mean_squared_error(y_train, ypred_train)
    print('train error',np.sqrt(train_er))
    print('test error',np.sqrt(test_er))
    tets_r2=r2_score(y_test, ypred)
    train_r2=r2_score(y_train, ypred_train)
    print('train r^2',train_r2)
    print('test r^2',tets_r2)
    param_lst=reg.coef_
    out_df=pd.DataFrame(param_lst)
    out_df.to_csv('params.csv',index=None,header=None)
    out_df=pd.DataFrame(param_lst)
    out_df=out_df
    out_df['features']=df_features.columns.values
    out_df.to_csv('params_w_header.csv',index=None)
    pred_train_df=pd.DataFrame(ypred_train)
    pred_train_df.to_csv('predtrain.csv',index=None,header=None)

    fig = plt.figure()
    plt.scatter(ypred_train,y_train)
    plt.scatter(ypred,y_test)
    xl=np.arange(75,150)
    plt.plot(xl,xl,'k-')
    #plt.xlim(min(target)-10,100)
    #plt.ylim(min(target)-10,100)
    #plt.scatter(ypred,y_test)
#    plt.xlabel('predicted yield \n'+combine_label, fontsize=10)
    ax = fig.add_subplot(111)
    print(ypred)
    plt.show()

    # Ploting part
