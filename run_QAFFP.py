import os,sys, os.path
from validation_metrics import *
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.rdBase
from rdkit import DataStructs
from rdkit.DataStructs import BitVectToText
from rdkit import DataStructs
from rdkit.Chem import Descriptors as Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from sklearn.model_selection import StratifiedKFold
from sklearn import preprocessing
import numpy as np
from sklearn import ensemble
from sklearn import metrics
from sklearn.utils import shuffle
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestRegressor
from scipy.stats import pearsonr

target_now = sys.argv[1]
type_tr=sys.argv[2]
repetition_now=int(sys.argv[3])
type_fps = sys.argv[4]
fp_length=int(sys.argv[5]) +1

suppl = Chem.SDMolSupplier("datasets/{}.sdf".format(target_now))
mols = [x for x in suppl if x is not None]
Kis = [x.GetProp("pIC50") for x in suppl if x is not None]
chembl_ids  = [x.GetProp("ChEMBL_ID") for x in suppl if x is not None]

if len(Kis) != len(mols):
    raise

# load descriptors (QAFFP)
if type_fps != "cp_binary_440":
    fil='./affps/{}_{}_models.csv'.format(target_now, type_fps)
    mol_ids= np.loadtxt(fil,delimiter=",",skiprows=1,dtype="|S440",usecols=0)
    mol_ids = [x.decode('UTF-8') for x in mol_ids]
    fps= np.loadtxt(fil,delimiter=",",skiprows=1,dtype="float",usecols=np.arange(1,fp_length ))
else:
    fil='./affps/{}_{}_models.csv'.format(target_now, type_fps)
    mol_ids= np.loadtxt(fil,delimiter=",",skiprows=1,dtype="|S440",usecols=0)
    mol_ids = [x.decode('UTF-8') for x in mol_ids]
    fil='./affps/{}_{}_models.csv_processed'.format(target_now, type_fps)
    fps= np.loadtxt(fil,delimiter=",",skiprows=1,dtype="float",usecols=np.arange(1,int(sys.argv[5])))

# guarantee same ordering
tokeep2=[]
for i in chembl_ids:
    tokeep2.append(fps[mol_ids.index(i)])

fps = np.asarray(tokeep2);

##--------------------------------------------------------
## calculate Morgan fingerprints
##--------------------------------------------------------
Morgan_fps = []; torm=[]
for i,sample in enumerate(mols):
    use = True
    try:
            fp = AllChem.GetMorganFingerprintAsBitVect(sample,2,nBits=1024)
    except:
            use = False
            torm.append(i)

    if use:
        Morgan_fps.append(fp)

Morgan_fps = np.array(Morgan_fps,dtype="float")
if Morgan_fps.shape[0] != len(Kis):
    raise "Dimensions do not match"

##--------------------------------------------------------
## physchem using mordred
##--------------------------------------------------------
if type_tr in ['physchem','QAFFP_physchem']:
    physchem = []
    descriptors = []
    
    # iterate over all available descriptots
    for D in Descriptors._descList:
        descriptors.append(D[0])   
    calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptors)
    
    for sample in mols:
        use = True
        try:
                pattern = calculator.CalcDescriptors(sample)
        except:
                use = False
               
        if use:
            physchem.append(pattern)
    physchem = np.array(physchem,dtype="float")    
    if physchem.shape[0] != len(Kis):
        raise "Dimensions do not match"

#--------------------------------------------------------
# center and scale
#--------------------------------------------------------
if type_tr == "physchem":
    fps = physchem
    del physchem
if type_tr == "Morgan":
    fps = Morgan_fps
    del Morgan_fps
if type_tr == "QAFFP":
    fps = fps
if type_tr == "QAFFP_Morgan":
    fps = np.column_stack((fps,Morgan_fps))
    del Morgan_fps
if type_tr == "QAFFP_physchem":
    fps = np.column_stack((fps,physchem))
    del physchem

fps = preprocessing.scale(fps)
Kis = np.asarray(Kis,dtype="float")

#--------------------------------------------------------
# Divide the original dataset into training and test set
#--------------------------------------------------------
for repetition in [repetition_now]:
    for test_size in [0.3]:
        if not os.path.isfile("results/test_"+type_tr+"_"+target_now+"_"+str(repetition)+"_"+str(test_size)):
             X_train, X_test, y_train, y_test, ids_train, ids_test = train_test_split(fps, Kis, chembl_ids, test_size=test_size, random_state=int(repetition))
             X_test = X_test.astype(float)
             X_train = X_train.astype(float)
             
             #--------------------------------------------------------
             # Train a Random Forest Model 
             # Do k-fold CV, each time calculating the alpha values
             #--------------------------------------------------------
             alphas = []; diffs=[]; ys=[]; zs=[]; stds=[]; pred_errors=[]; chembl_ids_CV=[]
             indexes = np.arange(len(y_train))
             np.random.seed(repetition)
             np.random.shuffle(indexes)
             stride= len(indexes)/10
             idx_used = []
             folds = []
             
             RF_model = RandomForestRegressor(n_estimators=100,random_state=23,n_jobs=1)
             RF_model.fit(X_train,y_train)
             ztest = RF_model.predict(X_test)
             diff = np.abs(y_test - ztest)
             
             f = open('results/test_'+type_tr+"_"+target_now+'_'+str(repetition)+"_"+str(test_size)+"_"+str(type_fps), 'w')
             f.write("obs\tpred\tdiff\tchembl_ids\n")
             for i in np.arange(0,len(y_test)):
                 f.write(str(y_test[i])+"\t"+str(ztest[i])+"\t"+str(diff[i])+"\t"+str(ids_test[i])+'\n')  
             f.close() 
             
             print(pearsonr(np.asarray(y_test), np.asarray(ztest)) )
             print (metrics.r2_score(np.asarray(y_test), np.asarray(ztest)) )
             metrics_now = ValidationMetrics(np.asarray(ztest),np.asarray(y_test))
             print ("RMSE: ", np.round(metrics_now.RMSE,decimals=2))
             print ("Q squared: ", np.round(metrics_now.Qsquared,decimals=2))
             print ("R squared: ", np.round(metrics_now.Rsquared,decimals=2))
             print ("R squared 0: ", np.round(metrics_now.Rsquared0,decimals=2))
             print ("Slope: ", np.round(metrics_now.slope,decimals=2))
