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



if __name__ == '__main__':
    df_features = pd.read_csv('ft.csv')
    df_data = pd.read_excel('CC_CC_data_Ph_Ph_yield_sub.xlsx')
    print('df_features',df_features)
    x=df_features

    target=df_data['Yield']
    scaler = StandardScaler()
    scaler.fit_transform(x)
    #x_train, x_test, y_train, y_test= train_test_split(x, target, test_size=3)

    #x = np.asarray(x)
    #param_grid = {"alpha": np.logspace(-5, 5, 20)}
    #estimator = GridSearchCV(linear_model.Lasso(), cv=3,param_grid=param_grid)
    #estimator.fit(x,target)
    #best_alpha=estimator.best_estimator_.alpha
    reg = linear_model.Lasso(alpha=0.000001)
    reg.fit(x,target)
    print(reg.coef_)
    #ypred = reg.predict(x_test)
    ypred_train = reg.predict(x)

    #print('RMSE',np.sqrt(mean_squared_error(ypred,y_test)))

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
    plt.scatter(ypred_train,target)
    xl=np.arange(0,100)
    plt.plot(xl,xl,'k-')
    plt.xlim(min(target)-10,100)
    plt.ylim(min(target)-10,100)
    #plt.scatter(ypred,y_test)
    plt.show()


    # Ploting part
