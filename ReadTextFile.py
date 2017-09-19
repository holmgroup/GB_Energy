#Imports
import numpy as np
import pandas as pd

#Read the text file
temp = pd.read_csv("fcc1-2.txt", delim_whitespace = True, skiprows = 6, header = None)

#Parse it into a readable column
main = pd.DataFrame(columns={0:'r1',1:'r2',2:'r3',4:'disorientation',5:'h1',6:'k1',7:'l1'})
i=0
while i < temp.shape[0]:
    #print temp.ix[i]
    if i==0:
        main =pd.DataFrame([temp.ix[i][0],temp.ix[i][1],temp.ix[i][2],temp.ix[i][4],temp.ix[i][5],temp.ix[i][6], temp.ix[i][7]])
    else:
        main = pd.concat([main,pd.DataFrame([temp.ix[i][0],temp.ix[i][1],temp.ix[i][2],temp.ix[i][4],temp.ix[i][5],temp.ix[i][6], temp.ix[i][7]])], axis=1)
    i=i+1
    try:
        while ~np.isnan(temp.ix[i][6]):
            i=i+1
        while np.isnan(temp.ix[i][6]):
            i=i+1
    except KeyError:
        continue

main = main.transpose().rename(columns={0:'r1',1:'r2',2:'r3',4:'disorientation',5:'h1',6:'k1',7:'l1'})
print main.head()
#Write to file
main.to_csv('GBEfeatures.csv')