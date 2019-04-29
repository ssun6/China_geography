import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import datasets,svm,metrics
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_curve,auc,accuracy_score
from sklearn.model_selection import StratifiedKFold, KFold,train_test_split
from sklearn.preprocessing import label_binarize
from sklearn.multiclass import OneVsRestClassifier
from itertools import cycle
from scipy import interp


df1 = pd.read_csv("/Users/shansun/Google Drive/China/data/China1/log_abundance_L6.csv",index_col=0, header = 0,delimiter=",")
map1 = pd.read_csv("/Users/shansun/Google Drive/China/Metadata_01082019/metadata_subset_01082019.csv",index_col=0, header = 0,delimiter=",")
map2=map1[["t1","t2"]]
df2=df1.transpose()
df=df2.join(map2)

#the prediction of province (t1)
X = df.iloc[:,0:1343]
X = np.c_[X]
y = pd.factorize(df['t1'])[0]
n_classes = 15

classifier = OneVsRestClassifier(RandomForestClassifier())
cv = StratifiedKFold(n_splits=6)

aucs = []
tps= []
for train, test in cv.split(X, y):
    for i in range(n_classes):
        y = label_binarize(y, classes=range(0,15))
        y_score= classifier.fit(X[train], y[train]).predict(X[test])
        fpr, tpr, _ = roc_curve(y[test][:, i], y_score[:, i])
        tp=accuracy_score(y[test][:, i], y_score[:, i])
        tps.append(tp)
        roc_auc = auc(fpr, tpr)
        aucs.append(roc_auc)
np.mean(tps)
np.std(tps)


from sklearn.ensemble import RandomForestRegressor
from sklearn.datasets import make_regression
from sklearn.metrics import mean_squared_error
pname=[11,21,23,31,32,33,37,41,42,43,45,52,53,55,61]


random_state = np.random.RandomState(0)
cv = KFold(n_splits=6)
regr = RandomForestRegressor(max_depth=2, random_state=0,n_estimators=100)

for j in [5]+list(range(9,67))+[69,70]+list(range(83,90)):
    print(map1.columns[j])
    map2=map1.iloc[:,[6,j]]
    df=df2.join(map2)
    df=df.dropna()
    
    k_rses=np.empty((15,6))
    pro_rses=np.empty((15,14))
    for m in range(15):
    
        df_n=df.loc[df["t1"]==pname[m],]

        X = df_n.iloc[:,0:1342]
        X1 = np.c_[X]
        y = df_n[map1.columns[j]]

        rses = []
        for train, test in cv.split(X1, y):
            ypred= regr.fit(X1[train], y[train]).predict(X1[test])
            rse=np.sqrt(mean_squared_error(ypred,y[test]))/np.mean(y[test])
            rses.append(rse)
        k_rses[m]=rses  
    
        rses1 = []
        for prov in np.delete(pname,m):
            df_m=df.loc[df["t1"]==prov,]
            X2 = df_m.iloc[:,0:1342]
            X3 = np.c_[X2]
            y1 = df_m[map1.columns[j]]
            ypred= regr.fit(X1, y).predict(X3)
            rse=np.sqrt(mean_squared_error(ypred,y1))/np.mean(y1)
            rses1.append(rse)
    
        pro_rses[m]=rses1
    all_rses=np.c_[k_rses, pro_rses]
    mname=map1.columns[j]
    np.savetxt("/Users/shansun/Documents/GitHub/RandomForest/RF_%s.csv"%mname, all_rses, delimiter=",")

random_state = np.random.RandomState(0)
cv = StratifiedKFold(n_splits=6)
Classifier= RandomForestClassifier()
from sklearn.metrics import accuracy_score
for j in [7,8,9,76,78,79,82]:
    print(map1.columns[j])
    map2=map1.iloc[:,[6,j]]
    df=df2.join(map2)
    df=df.dropna()
    
    k_tprs=np.empty((15,6))
    pro_tprs=np.empty((15,14))
    for m in range(15):
    
        df_n=df.loc[df["t1"]==pname[m],]
        X1 = df_n.iloc[:,0:1342]
        X = np.c_[X1]
        y = pd.factorize(df_n[map1.columns[j]])[0]
  
        tprs=[]
        for train, test in cv.split(X, y):       
            ypred = Classifier.fit(X[train], y[train]).predict(X[test])
            tpr= accuracy_score(y[test], ypred)
            tprs.append(tpr)
        k_tprs[m]=tprs
        
        
        tprs1 = []
        for prov in np.delete(pname,m):
            df_m=df.loc[df["t1"]==prov,]
            X2 = df_m.iloc[:,0:1342]
            X3 = np.c_[X2]
            y1 = pd.factorize(df_m[map1.columns[j]])[0]
            ypred = classifier.fit(X, y).predict(X3)
            tpr = accuracy_score(y1, ypred)
            tprs1.append(tpr)
    
        pro_tprs[m]=tprs1
    all_tprs=np.c_[k_tprs, pro_tprs]
    mname=map1.columns[j]
    np.savetxt("/Users/shansun/Documents/GitHub/RandomForest/RF_%s.csv"%mname, all_tprs, delimiter=",")

#shuffle the table and calculate the results for control
random_state = np.random.RandomState(0)
cv = StratifiedKFold(n_splits=6)
classifier = RandomForestClassifier()       
for j in [7,8,9,76,78,79,82]:
    print(map1.columns[j])
    map2=map1.iloc[:,[6,j]]
    df1=df2.join(map2)
    df1=df1.dropna()
    k_tprs=np.empty((150,6))
    pro_tprs=np.empty((150,14))
    n=0
    for i in range(10):   
        print(n)
        df=df1.copy()
        q=list(df[map1.columns[j]])#make a list of t2
        np.random.shuffle(q)#shuffle t2
        df[map1.columns[j]]=q #put the values back  
        for m in range(15):   
            df_n=df.loc[df["t1"]==pname[m],]
            X1 = df_n.iloc[:,0:1342]
            X = np.c_[X1]
            y = pd.factorize(df_n[map1.columns[j]])[0]
  
            tprs=[]
            for train, test in cv.split(X, y):       
                ypred = Classifier.fit(X[train], y[train]).predict(X[test])
                tpr= accuracy_score(y[test], ypred)
                tprs.append(tpr)
            k_tprs[n]=tprs
        
        
            tprs1 = []
            for prov in np.delete(pname,m):
                df_m=df.loc[df["t1"]==prov,]
                X2 = df_m.iloc[:,0:1342]
                X3 = np.c_[X2]
                y1 = pd.factorize(df_m[map1.columns[j]])[0]
                ypred = classifier.fit(X, y).predict(X3)
                tpr = accuracy_score(y1, ypred)
                tprs1.append(tpr)
            pro_tprs[n]=tprs1
            n=n+1
    all_tprs=np.c_[k_tprs, pro_tprs]
    mname=map1.columns[j]
    np.savetxt("/Users/shansun/Documents/GitHub/RandomForest/RF_shuffle_%s.csv"%mname, all_tprs, delimiter=",")

random_state = np.random.RandomState(0)
cv = KFold(n_splits=6)
egr = RandomForestRegressor(max_depth=2, random_state=0,n_estimators=100)

for j in [5]+list(range(9,67))+[69,70]+list(range(83,90)):
    print(map1.columns[j])
    map2=map1.iloc[:,[6,j]]
    df1=df2.join(map2)
    df1=df1.dropna()
    k_rses=np.empty((150,6))
    pro_rses=np.empty((150,14))
    n=0
    for i in range(10):   
        print(n)
        df=df1.copy()
        q=list(df[map1.columns[j]])#make a list of metadata
        np.random.shuffle(q)#shuffle metadata
        df[map1.columns[j]]=q #put the values back  
        for m in range(15):   
            df_n=df.loc[df["t1"]==pname[m],]
            X1 = df_n.iloc[:,0:1342]
            X = np.c_[X1]
            y = df_n[map1.columns[j]]
  
            rses=[]
            for train, test in cv.split(X, y):  
                ypred= regr.fit(X[train], y[train]).predict(X[test])
                rse=np.sqrt(mean_squared_error(ypred,y[test]))/np.mean(y[test])
                rses.append(rse)
            k_rses[n]=rses  
    
            rses1 = []
            for prov in np.delete(pname,m):
                df_m=df.loc[df["t1"]==prov,]
                X2 = df_m.iloc[:,0:1342]
                X3 = np.c_[X2]
                y1 = df_m[map1.columns[j]]
                ypred = regr.fit(X, y).predict(X3)
                rse=np.sqrt(mean_squared_error(ypred,y1))/np.mean(y1)
                rses1.append(rse)           
            pro_rses[n]=rses1
            n=n+1
    all_rses=np.c_[k_rses, pro_rses]
    mname=map1.columns[j]
    np.savetxt("/Users/shansun/Google Drive/China/figures/RF/RF_shuffle_%s.csv"%mname, all_rses, delimiter=",")
