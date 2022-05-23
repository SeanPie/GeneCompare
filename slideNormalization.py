import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import warnings
import copy
import time


#cleanup data aka remove outliers
warnings.filterwarnings("ignore")


def makeMatrix(dataset):
    # dsmatrix = pd.DataFrame()    
    # dsmatrix['y'] = dataset['y']
    # dsmatrix['x'] = dataset['x']
    # dsmatrix['tissue'] = dataset['tissue']
    return dataset.pivot('y', 'x', 'tissue')

#get 0 row counts
def getTop0Rows(dataset):
    dataset = dataset.sort_values(by=['y', 'x'])
    for index, row in dataset.iterrows():
        if row['tissue'] == 1:
            return row['y']

def getBottom0Rows(dataset):
    dataset = dataset.sort_values(by=['y'], ascending=[0])
    for index, row in dataset.iterrows():
        if row['tissue'] == 1:
            return row['y']

def getLeft0Rows(dataset):
    dataset = dataset.sort_values(by=['x', 'y'])
    for index, row in dataset.iterrows():
        if row['tissue'] == 1:
            return row['x']
            
def getRight0Rows(dataset):
    dataset = dataset.sort_values(by=['x'], ascending=[0])
    for index, row in dataset.iterrows():
        if row['tissue'] == 1:
            return row['x']

def moveRowstoTop(dataset, totalMovement):
    #remove rows
    newdf = dataset[dataset['y'] <= (dataset['y'].max() - totalMovement)]
    #adj cols
    newdf['y'] += totalMovement
    #addback old rows
    otherRows = dataset[dataset['y'] > (dataset['y'].max() - totalMovement)]
    otherRows['y'] -= (dataset['y'].max() - totalMovement) + 1

    allrows = pd.concat([otherRows, newdf])
   
    return allrows

def moveRowstoBottom(dataset, totalMovement):
    #remove rows
    newdf = dataset[dataset['y'] > totalMovement]
    #adj cols
    newdf['y'] -= totalMovement
    #addback old rows
    otherRows = dataset[dataset['y'] <= totalMovement]

    otherRows['y'] += (dataset['y'].max() - totalMovement) + 1

    allrows = pd.concat([newdf, otherRows])
   
    return allrows

def moveRowstoLeft(dataset, totalMovement):
    #remove rows
    newdf = dataset[dataset['x'] <= (dataset['x'].max() - totalMovement)]
    #adj cols
    newdf['x'] += totalMovement
    #addback old rows
    otherRows = dataset[dataset['x'] > (dataset['x'].max() - totalMovement)]
    otherRows['x'] -= (dataset['x'].max() - totalMovement) + 1

    allrows = pd.concat([otherRows, newdf])
   
    return allrows

def moveRowstoRight(dataset, totalMovement):
     #remove rows
    newdf = dataset[dataset['x'] > totalMovement]
    #adj cols
    newdf['x'] -= totalMovement
    #addback old rows
    otherRows = dataset[dataset['x'] <= totalMovement]

    otherRows['x'] += (dataset['x'].max() - totalMovement) + 1

    allrows = pd.concat([newdf, otherRows])
   
    return allrows

#Translation to center of matrix
def slidingWindow(dataset):
    bottomVal = getBottom0Rows(dataset)
    topVal = getTop0Rows(dataset)
    rightVal = getRight0Rows(dataset)
    leftVal = getLeft0Rows(dataset)

    if (dataset['y'].max() - bottomVal) > topVal:
        dataset = moveRowstoTop(dataset, round(((dataset['y'].max() - bottomVal) - topVal)/2) ) 
    else:
        dataset = moveRowstoBottom(dataset, round((topVal - (dataset['y'].max() - bottomVal))/2))


    if (dataset['x'].max() - rightVal) > leftVal:
        dataset = moveRowstoLeft(dataset, round(((dataset['x'].max() - rightVal) - leftVal)/2))
    else:
        dataset = moveRowstoRight(dataset, round((leftVal - (dataset['x'].max() - rightVal))/2))

    return dataset

def rotationofData2(dataset, radian):
    dataset2 = copy.deepcopy(dataset)
    dataset2['x'] =  round(dataset2['x'] * math.cos(radian) - dataset2['y'] * math.sin(radian), 0)
    dataset2['y'] = round(dataset2['y'] * math.cos(radian) + dataset2['x'] * math.sin(radian), 0)
    dataset2['x'] = dataset2['x'].astype(int)
    dataset2['y'] = dataset2['y'].astype(int)  
    dataset2 = dataset2.sort_values(by=['y', 'x'])
    
    dataset3 = dataset2.groupby(by=['x', 'y'])
    grouplist = []       
    for group_key, group_value in dataset3:
        dict1 = {'x': group_key[0], 'y': group_key[1], 'tissue': group_value['tissue'].max()}
        grouplist.append(dict1)

        
    newDF = pd.DataFrame(grouplist)
    return newDF

def dropLeftandRight(numRows, ds):
    halfNumRows = round(numRows/2)
    if numRows % 2 == 1:
        ds = ds[ds['x'] > (halfNumRows + 1)]
    else:
        ds = ds[ds['x'] > halfNumRows]
        
    ds = ds[ds['x'] <= (ds['x'].max() - halfNumRows)]
    return ds

def addLeftandRight(numRows, ds):
    for x in range(numRows):
        df = pd.DataFrame({'x' : [(ds['x'].max() + 1)], 'y': [0], 'tissue':[0]})
        ds = pd.concat([df, ds], ignore_index=True)
    return ds

def addTopandBot(numRows, ds):
    for x in range(numRows):
        df = pd.DataFrame({'x' : [0], 'y': [(ds['y'].max() + 1)], 'tissue':[0]})
        ds = pd.concat([df, ds], ignore_index=True)
    
    return ds

def dropTopandBot(numRows, ds):
    halfNumRows = round(numRows/2)
    if numRows % 2 == 1:
        ds = ds[ds['y'] > (halfNumRows + 1)]
    else:
        ds = ds[ds['y'] > halfNumRows]
    
    ds = ds[ds['y'] <= (ds['y'].max() - halfNumRows)]
    
    return ds

def zeroColumns(dataset):
    dataset = copy.deepcopy(dataset)
    
    #zero x values
    if dataset['x'].min() < 0:
        dataset['x'] += abs(dataset['x'].min()) 
    else:
        dataset['x'] -= abs(dataset['x'].min()) 
        
    #zero y values
    if dataset['y'].min() < 0:
        dataset['y'] += abs(dataset['y'].min()) 
    else:
        dataset['y'] -= abs(dataset['y'].min()) 
    return dataset

#just add to zero and then normalize for adding to columns

def reshapeDatatoSize(ds, x, y):
    ds1 = zeroColumns(ds)
    ds1 = slidingWindow(ds1)
    
    ydrop= y - ds1['y'].max()
    xdrop= x - ds1['x'].max()
    
    #drop columns to match base
    #     print("base: ", y, x)
    #     print(xdrop, x, ds1['x'].max())
    #drop columns
    if xdrop < 0:
        ds1 = dropLeftandRight(abs(xdrop), ds1)
    else:
        ds1 = addLeftandRight(xdrop, ds1)
    #         print("current: ", ds1['x'].max())   
        
    if ydrop < 0:
        ds1 = dropTopandBot(abs(ydrop), ds1)
    else:
        ds1 = addTopandBot(ydrop, ds1)   
    
    
    ds1 = zeroColumns(ds1)
    #     ds = dropLeftandRight(xdrop, ds)
    #     ds = dropTopandBot(ydrop, ds)

    return ds1

#gets the tissue value at the point
def getTissue(df, x, y):
    holding = df.loc[(df['y'] == y) & (df['x'] == x)]
    if holding.empty:
        return 0
    return int(holding['tissue'])

#collects surrounding tissues vals
def getSurroundingTissueAmount(dataset, row):
    count = 0
    for x in range(-2, 3):
        for y in range(-2, 3):
            count += getTissue(dataset, row['x'] + x, row['y'] + y) 
    
    return count

#removes valus that are not a part of the main tissue
def cleanSlide(dataset):
    for index, row in dataset.iterrows():
        if row['tissue'] == 1:
            #check if its apart of tissue
            if getSurroundingTissueAmount(dataset, row) < 5:
                dataset.at[index, 'tissue'] = 0
            
    return dataset

#changes x and y coordinates to the final spots
def getFinalDataset(dataset, degree):
    ds1 = copy.deepcopy(dataset)
    
    ds1 = slidingWindow(ds1)

    dsx2 = rotationofData2(ds1, math.radians(degree))
    dsx2 = reshapeDatatoSize(dsx2, dataset['x'].max(), dataset['y'].max())
    dsx2 = slidingWindow(dsx2)

    return dsx2

def comparison(ds, ds2, degree, bestDegreeAmountSame):
    count = 0
    for index, row in ds.iterrows():
        temp = ds2[ds2['x'] == row[0]][ds2['y'] == row[1]]
        try:
            if row[2] == temp['tissue']:
                count +=1
        except:
            continue
    if count > bestDegreeAmountSame:
        return count
    
    return 1


#main caller for the needed method
#returns degree for the best match of the datasets
def findMatch(dataset1, dataset2):
    ds1 = copy.deepcopy(dataset1)
    ds1 = cleanSlide(ds1)
    ds1 = slidingWindow(ds1)

    ds2 = copy.deepcopy(dataset2)
    ds2 = cleanSlide(ds2)
    ds2 = slidingWindow(ds2)

    bestDegreeAmountSame = 0
    bestdegree = 0
    matrix = makeMatrix(ds2)
    print("Running Rotations and Comparisons")
    for x in range(1,360):
        dsx2 = rotationofData2(ds1, math.radians(x))
        dsx2 = reshapeDatatoSize(dsx2, ds2['x'].max(), ds2['y'].max())
        try:
            newmatrix1 = makeMatrix(dsx2)
        except:
            print(dsx2)
            continue
        newmatrix1 = newmatrix1.fillna(0)

        nm = newmatrix1.lt(matrix)
        if sum(nm.sum()) > bestDegreeAmountSame:
            bestDegree = x
            bestDegreeAmountSame = sum(nm.sum())

    return getFinalDataset(ds1, bestdegree), ds2


