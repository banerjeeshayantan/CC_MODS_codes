# import warnings filter
from warnings import simplefilter
# ignore all future warnings
simplefilter(action='ignore', category=FutureWarning)
import pandas as pd
import numpy as np
from numpy import arange
from numpy import argmax
#from mlxtend.classifier import StackingClassifier
#from sklearn.ensemble import VotingClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import train_test_split
from numpy import arange
from sklearn.metrics import accuracy_score
import copy
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
#from sklearn.naive_bayes import MultinomialNB
#from sklearn.naive_bayes import GaussianNB
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
from sklearn.discriminant_analysis import QuadraticDiscriminantAnalysis
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
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler
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
#import xgboost as xgb
import itertools
from sklearn.naive_bayes import GaussianNB
from imblearn.metrics import sensitivity_score
from imblearn.metrics import specificity_score
from sklearn.metrics import balanced_accuracy_score
from sklearn.metrics import matthews_corrcoef

#Training data
gene_only_66099=pd.read_table('~/DGE_RMA_MODS_day3day7.txt')
clinical_only=pd.read_table('~/DGE_RMA_MODS_day3day7_clinical.txt')
clinical_cols=clinical_only['Propensity']
labs_train= gene_only_66099['Label']
labs_train=labs_train.astype('category').cat.codes
gene_only_66099=gene_only_66099.drop('Label',axis=1)
combined_train = pd.concat([gene_only_66099.reset_index(drop=True), clinical_cols], axis=1)
deg=pd.read_table('~/Top_feats_Day3Day7MODS_without_propensity.txt')
#deg=pd.read_table('Top_feats_Day3Day7MODS.txt')
genes_not_present = {'H1-2','PRG2','TCEAL9'} #genes not common across all four feature sets
deg_final = [ele for ele in deg if ele not in genes_not_present]
#Test data: 5882
gene_only=pd.read_table('~/External_Validation_Data_Day3Day7_5882.txt')
clinical_only=pd.read_table('~/External_Validation_Clinical_5882.txt')
clinical_cols=clinical_only['Propensity']
labs_test_5882= gene_only['Label']
labs_test_5882=labs_test_5882.astype('category').cat.codes
gene_only=gene_only.drop('Label',axis=1)
#combined_test_5882 = pd.concat([gene_only.reset_index(drop=True), clinical_cols], axis=1)
combined_test_5882=gene_only

#Test data: 144406
gene_only=pd.read_table('~/GSE144406_expression_data_0.96.txt')
labs_test_166640= gene_only['Label']
labs_test_166640=labs_test_166640.astype('category').cat.codes
gene_only=gene_only.drop('Label',axis=1)
combined_test_166640 = gene_only

#Test data: 10938_Prous
gene_only=pd.read_table('~/ext_vali_data_10938_Proulx.txt')
clinical_only=pd.read_table('~/ext_vali_data_clinical.txt')
clinical_cols=clinical_only['PRISM']
labs_test_10938= gene_only['Label']
labs_test_10938=labs_test_10938.astype('category').cat.codes
gene_only=gene_only.drop('Label',axis=1)
#combined_test = pd.concat([gene_only.reset_index(drop=True), clinical_cols], axis=1)
combined_test_10938=gene_only

#Test data: EMTAB-1548
gene_only=pd.read_table('~/E-MTAB-1548-all_genes.txt')
clinical_cols=gene_only['Propensity']
labs_test= gene_only['Mortality']
labs_test_1548=labs_test.astype('category').cat.codes
gene_only=gene_only.drop('Mortality',axis=1)
#gene_only=gene_only.drop('Propensity', axis=1)
#combined_test_1548 = pd.concat([gene_only.reset_index(drop=True), clinical_cols], axis=1)
combined_test_1548=gene_only

def model_tune_nb(X_samp, y_samp):
    nb = GaussianNB()
    grid_params = [{'var_smoothing': [1e-11, 1e-10, 1e-9]}]
    clf = GridSearchCV(nb, grid_params, verbose = 1, cv=3, n_jobs = 100)
    clf.fit(X_samp,y_samp)
    print("Tuning Done")
    best_model = clf.best_estimator_
    best_model.fit(X_samp, y_samp)
    return best_model

def model_tune_lda(X_samp, y_samp):
    lda = QuadraticDiscriminantAnalysis()
    #grid_params = [{'solver' : ['svd', 'lsqr', 'eigen']}]
    grid_params = [{'reg_param': [0.1, 0.2, 0.3, 0.4, 0.5]}]
    clf = GridSearchCV(lda, grid_params, verbose = 1, cv=3, n_jobs = 100)
    clf.fit(X_samp,y_samp)
    print("Tuning Done")
    best_model = clf.best_estimator_
    best_model.fit(X_samp, y_samp)
    return best_model

def model_tune_svm(X_samp, y_samp):
    svm = SVC(random_state=42, probability=True)
    grid_params= [
                  {'C': [1, 10, 100, 1000], 'kernel': ['linear']},
                  {'C': [1, 10, 100, 1000], 'gamma': [0.001, 0.0001], 'kernel': ['rbf']},
                 ]
    clf = GridSearchCV(svm, grid_params, verbose = 1, cv=3, n_jobs = 100)
    clf.fit(X_samp,y_samp)
    print("Tuning Done")
    best_model = clf.best_estimator_
    best_model.fit(X_samp, y_samp)
    return best_model

def model_tune_knn(X_samp, y_samp):
    knn=KNeighborsClassifier()
    grid_params = { 'n_neighbors' : [5,7,9,11,13,15],
               'weights' : ['uniform','distance'],
               'metric' : ['minkowski','euclidean','manhattan']}
    clf = GridSearchCV(knn, grid_params, verbose = 1, cv=3, n_jobs = 100)
    clf.fit(X_samp,y_samp)
    print("Tuning Done")
    best_model = clf.best_estimator_
    best_model.fit(X_samp, y_samp)
    return best_model

def model_tune_dt(X_samp, y_samp):
    dt=DecisionTreeClassifier(random_state=42)
    grid_params = {'max_depth':[2,4,6,8,10,12],
                   'min_samples_split' :np.linspace(0.1, 1.0, 10, endpoint=True),
                   'min_samples_leaf' : np.linspace(0.1, 0.5, 5, endpoint=True),
                   'max_features' :[None, 'sqrt', 'log2'],
                   'criterion': ["gini", "entropy"]}
    clf = GridSearchCV(dt, grid_params, verbose = 1, cv=3, n_jobs = 100)
    clf.fit(X_samp,y_samp)
    print("Tuning Done")
    best_model = clf.best_estimator_
    best_model.fit(X_samp, y_samp)
    return best_model

def model_tune_rf(X_samp, y_samp):
    rf=RandomForestClassifier(random_state=42)
    grid_params= { 
        'n_estimators':[10,20,21,30,50,100,200,500],
        'max_depth':[2,4,6,8,10],
                   'max_features': ['sqrt', 'log2'],
        'max_depth' : [2,3,5,7,9,11,13,15,17,19,21],
        'min_samples_leaf': [1,3,5],
        'min_samples_split': [2, 5,10,12],
    }
    clf = GridSearchCV(rf, grid_params, verbose = 1, cv=5, n_jobs = 100)
    clf.fit(X_samp,y_samp)
    print("Tuning Done")
    best_model = clf.best_estimator_
    best_model.fit(X_samp, y_samp)
    return best_model

def model_tune_et(X_samp, y_samp):
    rf=ExtraTreesClassifier(random_state=42)
    grid_params= { 
        'n_estimators':[5,10,20,21,30,50,100,200,500],
        'max_depth':[2,4,6,8,10],
                   'max_features': ['sqrt', 'log2'],
        'max_depth' : [2,3,5,7,9,11,13,15,17,19,21],
        'min_samples_leaf': [1,3,5],
        'min_samples_split': [2, 5,10,12],
    }
    clf = GridSearchCV(rf, grid_params, verbose = 1, cv=5, n_jobs = 100)
    clf.fit(X_samp,y_samp)
    print("Tuning Done")
    best_model = clf.best_estimator_
    best_model.fit(X_samp, y_samp)
    return best_model

def model_tune_ada(X_samp, y_samp):
    rf=AdaBoostClassifier(random_state=42)
    grid_params= {
        'n_estimators':[5,10,20,21,30,50,100,200,500],
        'learning_rate':[0.1,0.3,0.5,1,3,5,7,9],
        'algorithm':['SAMME', 'SAMME.R']
    }
    clf = GridSearchCV(rf, grid_params, verbose = 1, cv=5, n_jobs = 110)
    clf.fit(X_samp,y_samp)
    print("Tuning Done")
    best_model = clf.best_estimator_
    best_model.fit(X_samp, y_samp)
    return best_model

def print_metrics(name, X_test,y_test,X_samp, y_samp, threshold, best_model,feats):
    y_probs = best_model.predict_proba(X_test[X_samp.columns])[:,1]
    y_test_predictions = np.where(best_model.predict_proba(X_test[deg_final[0:feats]])[:,1] > threshold, 1, 0)
    fp_rate, tp_rate, thresh = roc_curve(y_test, y_probs)
    y_test_predictions = best_model.predict(X_test[X_samp.columns])
    auc = roc_auc_score(y_test, y_probs)
    sen= sensitivity_score(y_test, y_test_predictions, pos_label=1)
    spe=specificity_score(y_test,y_test_predictions,pos_label=1)
    acc=balanced_accuracy_score(y_test,y_test_predictions)
    mcc=matthews_corrcoef(y_test, y_test_predictions)
    #print("Name: ",name,"Sen: ",np.round(sen,3),", Spe: ",np.round(spe,3),", AUROC: ",np.round(auc,3),", Mcc: ",np.round(mcc,3))
    return(np.round(sen,3), np.round(spe,3), np.round(auc,3), np.round(mcc,3))

def convert_to_labels(pos_probs, threshold):
    return (pos_probs >= threshold).astype('int')

def draw_roc(y_test, y_probs, index, cl, cu, auc1):
    names = ['E-MTAB-10938 (Pediatric Validation Set)','GSE144406 (Pediatric Test Set)','E-MTAB-5882 (Adult Test Set)']

    fp_rate, tp_rate, thresh = roc_curve(y_test, y_probs)
    auc = roc_auc_score(y_test, y_probs)
    n_bootstraps = 1000
    rng_seed = 42  # control reproducibility
    bootstrapped_scores = []
    rng = np.random.RandomState(rng_seed)

    for k in range(n_bootstraps):
        # bootstrap by sampling with replacement on the prediction indices
        indices = rng.randint(0, len(y_probs), len(y_probs))
        if len(np.unique(y_test[indices])) < 2:
            continue
        score = roc_auc_score(y_test[indices], y_probs[indices])
        bootstrapped_scores.append(score)
                #print("Bootstrap #{} ROC area: {:0.3f}".format(k + 1, score))
    sorted_scores = np.array(bootstrapped_scores)
    sorted_scores.sort()
    confidence_lower = sorted_scores[int(0.05 * len(sorted_scores))]
    confidence_upper = sorted_scores[int(0.95 * len(sorted_scores))]
    plt.plot(fp_rate, tp_rate, label=names[index]+" [AUROC=" +str(round(auc1,3)) + 
                     " ("+ str(round(cl,3))+" - "+str(round(cu,3))+")]", linewidth = 4)


#Running validation analysis

# The best classifier-sampling technique combination was obtained using CV - ADA-RF
#This section is tuning only for number of features based on auroc as a metric only and using the validation set
auroc=[];  models=[]
sampling_names = ["RUS", "REDN", "CC", "IHT", "NearMiss", "EDN", "Tomek", "Allknn","CondensedNN", "OSS", "SMOTE", "ROS", "ADASYN", "KmeansSmote", "BorderlineSmote", "SVMSmote"]
#sampling_names = ["BorderlineSmote"]
#sampling_techniuqes=[BorderlineSMOTE(random_state=42)]
sampling_techniuqes = [RandomUnderSampler(random_state=42), RepeatedEditedNearestNeighbours(), ClusterCentroids(random_state=42), InstanceHardnessThreshold(random_state=42), NearMiss(), EditedNearestNeighbours(), TomekLinks(), AllKNN(), CondensedNearestNeighbour(random_state=42), OneSidedSelection(random_state=42), SMOTE(random_state=42), RandomOverSampler(random_state=42), ADASYN(random_state=42), KMeansSMOTE(random_state=42),BorderlineSMOTE(random_state=42), SVMSMOTE(random_state=42)]
classifier_names = ["ADA"]; model =[]
model_index = 0
df=pd.DataFrame()
feats = [5,10,15,20,25,30,35,45,50]; auroc=[]
#feats = [30]; auroc=[]

for n in range(0,len(feats)): 
    print("For ",feats[n]," genes: ") 
    #subsetiing n+1 genes and sampling technique -> REDN 
    X_train_sel = gene_only_66099[deg_final[0:feats[n]]] #GSE66099 as training set
    scaler = [MinMaxScaler(), StandardScaler(), RobustScaler()]
    scaler_names = ["minmax", "std", "robust"]
    for s in range(0, len(scaler)):
        print("Scaling scheme:", scaler_names[s])
        for samp in range(0,len(sampling_techniuqes)):
            print("Sampling: ", sampling_names[samp], "Scaling: ", scaler_names[s], "Classifier: ", classifier_names[0])
            X_train_sel = pd.DataFrame(scaler[s].fit_transform(X_train_sel),columns = X_train_sel.columns)
            X_samp,y_samp=sampling_techniuqes[samp].fit_resample(X_train_sel,labs_train)
            X_samp = pd.DataFrame(X_samp, columns = X_train_sel.columns)
            name_id = "Features_"+str(feats[n])+"_"+ scaler_names[s]+"_"+sampling_names[samp]+"_"+classifier_names[0]
            #tuning KNN classifier
            best_model = model_tune_ada(X_samp, y_samp)
            print(best_model)
            model.append(best_model)
            #print(best_model_et.get_params())
            print("Training done")
             #results on 10938 - validation set
            y_test=labs_test_10938
            X_test=combined_test_10938  #EMTAB10938 as validation 
            X_test = pd.DataFrame(scaler[s].fit_transform(X_test),columns = X_test.columns)
            y_probs = best_model.predict_proba(X_test[X_samp.columns])[:,1]
            thresholds = arange(0, 1, 0.001)
            scores = [roc_auc_score(y_test, convert_to_labels(y_probs, t)) for t in thresholds]
            ix= argmax(scores)
            print("Threshold: ", thresholds[ix])
            y_test_predictions = np.where(best_model.predict_proba(X_test[deg_final[0:feats[n]]])[:,1] > thresholds[ix], 1, 0)
            fp_rate, tp_rate, thresh = roc_curve(y_test, y_probs)
            y_test_predictions = best_model.predict(X_test[X_samp.columns])
            auc = roc_auc_score(y_test, y_probs)
            auroc.append(auc)
            print("AUROC: ", auc)
            sen= sensitivity_score(y_test, y_test_predictions, pos_label=1)
            spe=specificity_score(y_test,y_test_predictions,pos_label=1)
            acc=balanced_accuracy_score(y_test,y_test_predictions)
            mcc=matthews_corrcoef(y_test, y_test_predictions)
            l=[model_index, name_id, sen, spe, acc, mcc, auc, thresholds[ix]]
            row=pd.Series(l,['Model_index','name','sen','spe','acc','mcc','auc','threshold'])
            df=df.append([row],ignore_index=True)
            model_index = model_index + 1
            print("Sensitivity: ", sen, "Specificity: ", spe, "Balanced Accuracy Score: ", acc, "MCC: ", mcc,"AUC: ", auc)


            
scaler = [MinMaxScaler(), StandardScaler(), RobustScaler()]
Sen_5882=[]; Spe_5882=[];AUC_5882=[];MCC_5882=[];Sen_1548=[];Spe_1548=[];AUC_1548=[];MCC_1548=[];
Sen_144406=[];Spe_144406=[];AUC_144406=[];MCC_144406=[];

df = df.sort_values(by = ['auc'], ascending=False)
for i in range(0,df.shape[0]):
    mindex = df.Model_index.to_list()[i]
    feats = df.name[mindex].split("_")[1]
    sc = df.name[mindex].split("_")[2]
    sampler = df.name[mindex].split("_")[3]
    index_scaler = scaler_names.index(sc)
    index_sampler = sampling_names.index(sampler)
    c = sampling_techniuqes[index_sampler]
    scaling = scaler[index_scaler]
    threshold = df.threshold[mindex]
    print(c,scaling, threshold,feats)
    X_train_sel = gene_only_66099[deg_final[0:int(feats)]]
    X_train_model = X_train_sel
    X_train_model = pd.DataFrame(scaling.fit_transform(X_train_model),columns = X_train_model.columns)
    X_samp,y_samp=c.fit_resample(X_train_model,labs_train)
    X_samp = pd.DataFrame(X_samp, columns = X_train_model.columns)
    best_model = model[mindex]
    sen_5882, spe_5882, auc_5882, mcc_5882 = print_metrics("Adult cohort 5882|",pd.DataFrame(scaling.fit_transform(combined_test_5882),columns = combined_test_5882.columns),labs_test_5882, X_samp, y_samp, float(threshold), best_model,int(feats))
    Sen_5882.append(sen_5882);Spe_5882.append(spe_5882);AUC_5882.append(auc_5882);MCC_5882.append(mcc_5882);
    sen_1548, spe_1548, auc_1548, mcc_1548 = print_metrics("Adult cohort 1548|",pd.DataFrame(scaling.fit_transform(combined_test_1548),columns = combined_test_1548.columns),labs_test_1548, X_samp, y_samp, float(threshold), best_model,int(feats))
    Sen_1548.append(sen_1548);Spe_1548.append(spe_1548);AUC_1548.append(auc_1548);MCC_1548.append(mcc_1548);
    sen_144406, spe_144406, auc_144406, mcc_144406 = print_metrics("Adult cohort 144406|",pd.DataFrame(scaling.fit_transform(combined_test_166640),columns = combined_test_166640.columns),labs_test_166640, X_samp, y_samp, float(threshold), best_model,int(feats))
    Sen_144406.append(sen_144406);Spe_144406.append(spe_144406);AUC_144406.append(auc_144406);MCC_144406.append(mcc_144406);

df['Sen_5882'] = Sen_5882; df['Spe_5882'] = Spe_5882; df['AUC_5882']=AUC_5882; df['MCC_5882'] = MCC_5882
df['Sen_1548'] = Sen_1548; df['Spe_1548'] = Spe_1548; df['AUC_1548']=AUC_1548; df['MCC_1548'] = MCC_1548
df['Sen_144406'] = Sen_144406; df['Spe_144406'] = Spe_144406; df['AUC_144406']=AUC_144406; df['MCC_144406'] = MCC_144406

df = df.sort_values(by = ['auc'], ascending=False)
filename = '~/Ablation_study_with_10938_as_validation/'+classifier_names[0]+".csv"
df.to_csv(filename,index=False)

###################################################################### Reproducibility of external validation analysis ##############################################################################################################################
#based on the validation results, top 20 features, threshold 0f 0.488 and ET classifier gave the best results. Now testing the model on unseen cohorts: GSE144406 and E-MTAB-5882
feats = 20
X_train_sel = gene_only_66099[deg_final[0:feats]]
scaling = StandardScaler()
c = InstanceHardnessThreshold(random_state=42)
threshold = 0.488
X_train_model = X_train_sel
X_train_model = pd.DataFrame(scaling.fit_transform(X_train_model),columns = X_train_model.columns)
X_samp,y_samp=c.fit_resample(X_train_model,labs_train)
X_samp = pd.DataFrame(X_samp, columns = X_train_model.columns)
best_model = pickle.load(open("best_model_ET.pkl", "rb"))
sen_10938, spe_10938, auc_10938, mcc_10938 = print_metrics("Pediatric cohort 10938|", pd.DataFrame(scaling.fit_transform(combined_test_10938),columns = combined_test_10938.columns),labs_test_10938, X_samp, y_samp, threshold, best_model,feats)
sen_5882,spe_5882, auc_5882, mcc_5882 = print_metrics("Adult cohort 5882|",pd.DataFrame(scaling.fit_transform(combined_test_5882),columns = combined_test_5882.columns),labs_test_5882, X_samp, y_samp, threshold, best_model,feats)
sen_1548, spe_1548, auc_1548, mcc_1548 = print_metrics("Adult cohort 1548|", pd.DataFrame(scaling.fit_transform(combined_test_1548),columns = combined_test_1548.columns),labs_test_1548, X_samp, y_samp, threshold, best_model,feats)
sen_144406, spe_144406, auc_144406, mcc_144406 = print_metrics("Pediatric cohort 144406|", pd.DataFrame(scaling.fit_transform(combined_test_166640),columns = combined_test_166640.columns),labs_test_166640, X_samp, y_samp, threshold, best_model,feats)

plt.figure(figsize=(15,15))
classifier = ["ET"]
#best_model = pd.read_pickle('best_model_ADA.pkl')
feats = 20
scaling = StandardScaler()
c = InstanceHardnessThreshold(random_state=42)
threshold = 0.488
X_train_sel = gene_only_66099[deg_final[0:feats]]
X_train_model = pd.DataFrame(scaling.fit_transform(X_train_sel),columns = X_train_sel.columns)
X_samp,y_samp=c.fit_resample(X_train_model,labs_train)
X_samp = pd.DataFrame(X_samp, columns = X_train_model.columns)

X_test_10938 = pd.DataFrame(scaling.fit_transform(combined_test_10938),columns = combined_test_10938.columns)
X_test_144406 = pd.DataFrame(scaling.fit_transform(combined_test_166640),columns = combined_test_166640.columns)
X_test_5882 = pd.DataFrame(scaling.fit_transform(combined_test_5882),columns = combined_test_5882.columns)
#X_test_1548 = pd.DataFrame(scaling.fit_transform(combined_test_1548),columns = combined_test_1548.columns)
y_test_10938 = labs_test_10938
y_test_144406 = labs_test_166640
y_test_5882 = labs_test_5882
#y_test_1548 = labs_test_1548


y_probs_10938 = best_model.predict_proba(X_test_10938[X_samp.columns])[:,1]
y_probs_144406 = best_model.predict_proba(X_test_144406[X_samp.columns])[:,1]
y_probs_5882 = best_model.predict_proba(X_test_5882[X_samp.columns])[:,1]
#y_probs_1548 = best_model.predict_proba(X_test_1548[X_samp.columns])[:,1]
#use the AUROC values reported in lines 431 to 434
draw_roc(y_test_10938, y_probs_10938, index = 0 , cl = 0.739, cu = 0.753, auc1 = auc_10938) 
draw_roc(y_test_144406, y_probs_144406, index = 1, cl = 0.785, cu = 0.801, auc1 = auc_144406 )
draw_roc(y_test_5882, y_probs_5882, index = 2, cl = 0.772, cu = 0.789, auc1 = auc_5882)
#draw_roc(y_test_1548, y_probs_1548, index = 3 )

plt.plot([0,1], [0,1], color='orange', linestyle='--', linewidth = 3)

plt.xticks(np.arange(0.0, 1.1, step=0.2), fontsize = 28)
plt.xlabel("False Positive Rate", fontsize=28)

plt.yticks(np.arange(0.0, 1.1, step=0.2), fontsize = 28)
plt.ylabel("True Positive Rate", fontsize=28)

    #plt.title('ROC Curve', fontweight='bold', fontsize=15)
plt.legend(prop={'size':15, 'weight':'bold'}, loc='lower right')
plt.savefig('/home/shayantan/complex_course_new_actinium/Ablation_study_with_10938_as_validation/ROC_ext_modified'+classifier[0]+'.pdf')
plt.show()


