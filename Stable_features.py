# import warnings filter
from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from numpy import arange
from sklearn.metrics import accuracy_score
import matplotlib.pyplot as plt
from sklearn.feature_selection import SelectFromModel
from sklearn.ensemble import ExtraTreesClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.ensemble import AdaBoostClassifier
from numpy import argmax
from sklearn import metrics
from imblearn.metrics import sensitivity_score
from imblearn.metrics import specificity_score
from sklearn.metrics import roc_auc_score
from sklearn.metrics import accuracy_score
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import MultinomialNB
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import LinearSVC, SVC
from sklearn.neighbors import KNeighborsClassifier
from sklearn.model_selection import cross_val_score
from sklearn.model_selection import cross_validate
from sklearn.model_selection import GridSearchCV
from sklearn.ensemble import RandomForestClassifier, VotingClassifier
from sklearn import preprocessing
from sklearn.preprocessing import MinMaxScaler
import pickle
from sklearn.metrics import confusion_matrix
from sklearn.utils import resample
from imblearn.over_sampling import SMOTE
from sklearn.metrics import classification_report
from sklearn import metrics
from sklearn.metrics import roc_curve, auc
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import cross_val_predict
from sklearn.model_selection import StratifiedKFold
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn import datasets
from sklearn.model_selection import train_test_split
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import accuracy_score
from sklearn.feature_selection import RFECV
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
import scipy.stats as ss
import itertools
from sklearn.metrics import fbeta_score, make_scorer
from sklearn.manifold import TSNE
from sklearn.multiclass import OneVsRestClassifier
from sklearn.model_selection import LeaveOneOut
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
import scipy
import numpy as np, scipy.stats as st
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.neural_network import MLPClassifier
from scipy import interp
from collections import Counter
from sklearn.model_selection import StratifiedKFold
from sklearn.ensemble import BaggingClassifier
from imblearn.ensemble import BalancedBaggingClassifier
from imblearn.ensemble import BalancedRandomForestClassifier
from imblearn.ensemble import EasyEnsembleClassifier
from imblearn.ensemble import RUSBoostClassifier
from imblearn.under_sampling import ClusterCentroids 
from sklearn.linear_model import RidgeCV, LassoCV, Ridge, Lasso
from imblearn.under_sampling import RepeatedEditedNearestNeighbours
from imblearn.under_sampling import EditedNearestNeighbours
from imblearn.under_sampling import RepeatedEditedNearestNeighbours
from imblearn.under_sampling import AllKNN
from imblearn.under_sampling import CondensedNearestNeighbour
from imblearn.under_sampling import TomekLinks
from imblearn.under_sampling import InstanceHardnessThreshold
from imblearn.under_sampling import NearMiss
from imblearn.under_sampling import NeighbourhoodCleaningRule
from imblearn.under_sampling import     OneSidedSelection 
from imblearn.under_sampling import RandomUnderSampler
from imblearn.metrics import classification_report_imbalanced
from imblearn.under_sampling import ClusterCentroids
from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import ADASYN
from imblearn.over_sampling import RandomOverSampler
from imblearn.over_sampling import SVMSMOTE
from imblearn.over_sampling import BorderlineSMOTE
from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import SVMSMOTE
from imblearn.over_sampling import SMOTENC
from imblearn.over_sampling import ADASYN
from imblearn.over_sampling import KMeansSMOTE
from imblearn.over_sampling import RandomOverSampler
from imblearn.over_sampling import BorderlineSMOTE
from imblearn.over_sampling import SMOTE
from imblearn.over_sampling import SVMSMOTE
from imblearn.over_sampling import SMOTENC
from imblearn.over_sampling import ADASYN
from imblearn.over_sampling import KMeansSMOTE
from imblearn.over_sampling import RandomOverSampler
from imblearn.over_sampling import BorderlineSMOTE
import pymrmr
import itertools
from imblearn.metrics import sensitivity_score
from imblearn.metrics import specificity_score
from sklearn.metrics import matthews_corrcoef


#make sure the following files are in the current working directory
gene_only=pd.read_table('DGE_RMA_MODS_day3day7.txt')
clinical_only=pd.read_table('DGE_RMA_MODS_day3day7_clinical.txt')
clinical_cols=clinical_only['Propensity']
labs= gene_only['Label']
labs=labs.astype('category').cat.codes
gene_only=gene_only.drop('Label',axis=1)
combined = pd.concat([gene_only.reset_index(drop=True), clinical_cols], axis=1)
deg=pd.read_table('DEG_day3day7.txt')
deg1=deg['X'].to_list()

def model_fit_and_validate_xgb(X,y,n,sampling,model,deg1):
    tprs = []
    base_fpr = np.linspace(0, 1, 101)
    idx = np.arange(0, len(y))

    plt.figure(figsize=(5, 5))
   
    rskf=RepeatedStratifiedKFold(n_splits=5,n_repeats=n)
    i=0
    auc=[];sensitivity=[];specificity=[];m=[];accuracy=[];cols=[]
    for train_index,test_index in rskf.split(X,y):
        i=i+1
        print(i,end="")
        with open('GSE66099_CC_RF_keep_track', 'w') as f:
            print(i, file=f)
        X_train,X_test=X.iloc[train_index],X.iloc[test_index]
        y_train,y_test=y.iloc[train_index],y.iloc[test_index]
        forest =  ExtraTreesClassifier(n_estimators=500,
                              random_state=42)
        scaler = MinMaxScaler()
        X_train = pd.DataFrame(scaler.fit_transform(X_train),columns = X_train.columns)
        X_test = pd.DataFrame(scaler.transform(X_test),columns = X_test.columns)
        forest.fit(X_train, y_train)
        importances = forest.feature_importances_
        std = np.std([tree.feature_importances_ for tree in forest.estimators_],
             axis=0)
        indices = np.argsort(importances)[::-1]
        rf_list=X_train.columns[indices][0:50] #Top 50 ranked RF features
        print("RF done")
        reg = LassoCV(random_state=42, n_jobs=110)
        reg.fit(X_train, y_train)
        coef = pd.Series(reg.coef_, index = X_train.columns)
        selected_f=X_train.columns[coef!=0] #Lasso features
        print("LAsso done")
        #Print the feature ranking
        df=X_train
        df.insert(loc=0, column='label', value=y_train)
        df['label']=df['label'].astype('category').cat.codes
        mrmr=pymrmr.mRMR(df, 'MIQ',20) #MRMR features
        merged_list = list(itertools.chain(*itertools.zip_longest(rf_list,selected_f,mrmr,deg1)))
        merged_list = [i for i in merged_list if i is not None]
        print("RF\n")
        print(rf_list)
        print("LASSO\n")
        print(selected_f)
        print("MRMR\n")
        print(mrmr)
        print("DEG\n")
        print(deg1)
       
        X_train_sel=X_train[set(merged_list)]
        #print("MRMR done") 
        
       # print(sampling)
        #sampling 
        X_samp,y_samp=sampling.fit_resample(X_train_sel,y_train.ravel())
        X_samp = pd.DataFrame(X_samp, columns = X_train_sel.columns)
        
        #model fitting
        if model=="BRF":
            sen, spe, auc_score,feat,acc,mcc,fpr,tpr,thresholds=BRF(X_samp,y_samp,X_test,y_test)
        elif model=="LOGIT":
            sen, spe, auc_score,feat,acc,mcc,fpr,tpr,thresholds=LOGIT(X_samp,y_samp,X_test,y_test)
        elif model=="SVM":
            sen, spe, auc_score,feat,acc,mcc=SVM(X_samp,y_samp,X_test,y_test)
        elif model=="ADABoost":
            sen, spe, auc_score,feat,acc,mcc=ADA(X_samp,y_samp,X_test,y_test)
        elif model=="ET":
            sen, spe, auc_score,feat,acc,mcc,fpr,tpr,thresholds=ET(X_samp,y_samp,X_test,y_test)
        elif model=="NaiveBayes":
            sen, spe, auc_score,feat,acc,mcc=NB(X_samp,y_samp,X_test,y_test)
        elif model=="LDA":
            sen, spe, auc_score,feat,acc,mcc=LDA(X_samp,y_samp,X_test,y_test)

        plt.plot(fpr, tpr, 'b', alpha=0.05)
        tpr = interp(base_fpr, fpr, tpr)
        tpr[0] = 0.0
        tprs.append(tpr)
        auc.append(auc_score)
        sensitivity.append(sen)
        specificity.append(spe)
        accuracy.append(acc)
        m.append(mcc)
        cols.append(feat)
        print("AUC", auc_score," ", "sen", sen, " ","spe",spe," ","acc",acc," ","mcc",mcc)
    tprs = np.array(tprs)
    mean_tprs = tprs.mean(axis=0)
    std = tprs.std(axis=0)
    tprs_upper = np.minimum(mean_tprs + std, 1)
    tprs_lower = mean_tprs - std  
    l_ci=mean_confidence_interval(auc)[1]
    r_ci=mean_confidence_interval(auc)[2]
    plt.plot(base_fpr, mean_tprs, 'b',color='darkorange',
         #lw=2,label=  '(AUC = %0.2f (CI: %0.2f - %0.2f)'  % (round(np.mean(auc),2),round(l_ci,2),round(r_ci,2)))
    plt.legend(loc='lower right')
    plt.fill_between(base_fpr, tprs_lower, tprs_upper, color='grey', alpha=0.3)
    plt.savefig('roc.pdf')
    plt.plot([0, 1], [0, 1],'r--')
    plt.xlim([-0.01, 1.01])
    plt.ylim([-0.01, 1.01])
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.axes().set_aspect('equal', 'datalim')
    plt.show()
    l=[np.mean(auc),mean_confidence_interval(auc)[1],mean_confidence_interval(auc)[2],np.mean(sensitivity),mean_confidence_interval(sensitivity)[1],mean_confidence_interval(sensitivity)[2],np.mean(spe),mean_confidence_interval(specificity)[1],mean_confidence_interval(specificity)[2],np.mean(m),mean_confidence_interval(m)[1],mean_confidence_interval(m)[2]]
    row=pd.Series(l,['auc','auc_lci','auc_rci','sen','sen_lci','sen_rci','spe','spe_lci','spe_rci','m','m_lci','m_rci'])
    df=pd.DataFrame()
    df=df.append([row],ignore_index=True)
    df.to_csv('class.csv',index=False)
    flat_list = [item for sublist in cols for item in sublist]
    x=Counter(flat_list)
    with open('Counter.txt', 'w') as f:
        print(x, file=f)
    return auc, sensitivity, specificity,accuracy,m,cols

def convert_to_labels(probs, thresh):
    return (probs >= thresh).astype('int')

def BRF(X_train,y_train,X_test,y_test):
    brf=RandomForestClassifier(class_weight="balanced_subsample")
    X_reduced=feature_reduction_using_RFECV(brf,X_train,y_train)
    param_grid = { 
    'n_estimators':[100,200,300,400,600,800,1000],
    'max_features': ['auto', 'sqrt', 'log2'],
    'max_depth' : [2,3,5,7,9,11],
}
    clf_brf = GridSearchCV(brf, param_grid, cv=3,scoring="roc_auc",n_jobs=-1)
    best_model_brf=clf_brf.fit(X_reduced,y_train)
    best_model_brf.fit(X_reduced, y_train)
    y_probs = best_model_brf.predict_proba(X_test[X_reduced.columns])[:,1]
    thresholds = arange(0, 1, 0.001)
    scores = [roc_auc_score(y_test, convert_to_labels(y_probs, t)) for t in thresholds]
    ix= argmax(scores)
    y_test_predictions = np.where(best_model_brf.predict_proba(X_test[X_reduced.columns])[:,1] > thresholds[ix], 1, 0)
    fp_rate, tp_rate, thresholds = roc_curve(y_test, y_probs)
    #auc_score=metrics.auc(fp_rate, tp_rate)
    auc_score=roc_auc_score(y_test,y_test_predictions)
    # y_test_predictions=best_model_dt.predict(X_test[X_reduced.columns])
    sen= sensitivity_score(y_test, y_test_predictions, pos_label=0)
    spe=specificity_score(y_test,y_test_predictions,pos_label=0)
    acc=accuracy_score(y_test,y_test_predictions)
    mcc=matthews_corrcoef(y_test, y_test_predictions)
    return sen, spe, auc_score, X_reduced.columns,acc,mcc,fp_rate,tp_rate,thresholds

def LOGIT(X_train,y_train,X_test,y_test):
    logit=LogisticRegression(class_weight="balanced")
    X_reduced=feature_reduction_using_RFECV(logit,X_train,y_train)
    penalty = ['l1', 'l2']
    # Create regularization hyperparameter space
    C = np.logspace(-3, 3, 10)
    # Create hyperparameter options
    hyperparameters = dict(C=C, penalty=penalty)
    clf_logit = GridSearchCV(logit, hyperparameters, cv=3, scoring="roc_auc",n_jobs=110)
    best_model_logit = clf_logit.fit(X_reduced,y_train)
    best_model_logit.fit(X_reduced, y_train)
    y_probs = best_model_logit.predict_proba(X_test[X_reduced.columns])[:,1]
    fp_rate, tp_rate, thresholds = roc_curve(y_test, y_probs)
    auc_score=metrics.auc(fp_rate, tp_rate)
    y_test_predictions=best_model_logit.predict(X_test[X_reduced.columns])
    sen= sensitivity_score(y_test, y_test_predictions, pos_label=0, average='binary')
    spe=specificity_score(y_test,y_test_predictions,pos_label=0,average='binary')
    acc=accuracy_score(y_test,y_test_predictions)
    mcc=matthews_corrcoef(y_test, y_test_predictions)
    return sen, spe, auc_score, X_reduced.columns,acc,mcc,fp_rate,tp_rate,thresholds

def ET(X_train,y_train,X_test,y_test):

    dt=ExtraTreesClassifier(class_weight="balanced_subsample")
    X_reduced=feature_reduction_using_RFECV(dt,X_train,y_train)
    parameters={ 'n_estimators':[100,200,300,400,600,800,1000,1200],'min_samples_split' : [2,3,4,5,6,7],'max_depth': [2,3,5,7], 'max_features':['auto','sqrt','log2']}
    clf_dt = GridSearchCV(dt, parameters, cv=3,scoring="roc_auc",n_jobs=110)
    best_model_dt=clf_dt.fit(X_reduced,y_train)
    best_model_dt.fit(X_reduced, y_train)
    y_probs = best_model_dt.predict_proba(X_test[X_reduced.columns])[:,1]
    thresholds = arange(0, 1, 0.001)
    scores = [roc_auc_score(y_test, convert_to_labels(y_probs, t)) for t in thresholds]
    ix= argmax(scores)
    y_test_predictions = np.where(best_model_dt.predict_proba(X_test[X_reduced.columns])[:,1] > thresholds[ix], 1, 0)
    fp_rate, tp_rate, thresholds = roc_curve(y_test, y_probs)
    #auc_score=metrics.auc(fp_rate, tp_rate)
    auc_score=roc_auc_score(y_test,y_test_predictions)
    # y_test_predictions=best_model_dt.predict(X_test[X_reduced.columns])
    sen= sensitivity_score(y_test, y_test_predictions, pos_label=0)
    spe=specificity_score(y_test,y_test_predictions,pos_label=0)
    acc=accuracy_score(y_test,y_test_predictions)
    mcc=matthews_corrcoef(y_test, y_test_predictions)
    return sen, spe, auc_score, X_reduced.columns,acc,mcc,fp_rate,tp_rate,thresholds

def feature_reduction_using_RFECV(model,X,y):
    feature_list=[]
    np.random.seed(315)
    
    rfecv = RFECV(estimator=model, step=1, cv=3, scoring='roc_auc',n_jobs=110)
    rfecv.fit(X, y)
    new_X=X[X.columns[rfecv.support_==True]]
    return new_X


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2., n-1)
    return m, m-h, m+h


##DONOT RUN Takes a very long time
cc=ADASYN(random_state=42)
auc, sensitivity, specificity,accuracy,m,cols=model_fit_and_validate(combined, labs, 7, cc, "BRF",deg1)
l=[np.mean(auc),mean_confidence_interval(auc)[1],mean_confidence_interval(auc)[2],np.mean(sensitivity),mean_confidence_interval(sensitivity)[1],mean_confidence_interval(sensitivity)[2],np.mean(specificity),mean_confidence_interval(specificity)[1],mean_confidence_interval(specificity)[2],np.mean(m),mean_confidence_interval(m)[1],mean_confidence_interval(m)[2]]
row=pd.Series(l,['auc','auc_lci','auc_rci','sen','sen_lci','sen_rci','spe','spe_lci','spe_rci','m','m_lci','m_rci'])
df=pd.DataFrame()
df=df.append([row],ignore_index=True)
#df.to_csv('/data/shayantan/Normalized_Data_for_training/GSE152075_results/GSE152075_CC_LOGIT_30.csv',index=False)
flat_list = [item for sublist in cols for item in sublist]
x=Counter(flat_list)
print(df)
