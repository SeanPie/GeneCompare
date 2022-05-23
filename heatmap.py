#!/usr/bin/env python3
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import argparse
import os
import csv
import cmath
import time

#function to rename ensembl gene name to more common shorthand gene name 
def genelookup():
    genelist = {}
    filepath1 = "{0}".format(args.Features)
    try:
        #file from outputs of SpaceRanger pipeline, first column is ensembl name and second column is gene same
        with open(filepath1, "r") as gene_file:
            for line in gene_file:
                line = line.rstrip()
                gene_split = line.split("\t")
                #makes sure its only copying gene expression gene names over
                if gene_split[2] == "Gene Expression":
                    #key is ensembl name, value is shorthand
                    genelist[gene_split[0]] = gene_split[1]
                else:
                    continue
    except IOError:
        print("Error reading file!")

    return genelist


#function to create lookup dictionary for y, x, genetype, and gene difference between the two datasets 
def dataget(genelist):
    matrix_dict = {}
    filepath2 = "{0}".format(args.CSV)
    try:
        #file that comes out of GeneCompare.py code with all above data
        with open(filepath2, "r") as dataset:
            #removes header line
            first_line = dataset.readline()
            for line in dataset:
                line = line.rstrip()
                line_split = line.split(",")
                #triple nested dictionary; first layer is key=Y coord, value is empty dictionary 
                matrix_dict.setdefault(int(line_split[0]), dict())
                #second layer is Key=X coord in last empty dictionary, value empty dictionary
                matrix_dict[int(line_split[0])].setdefault(int(line_split[1]), dict())
                #final dictionary Key= lookup value from gene list dictionary  and value is the floating point value of gene difference
                matrix_dict[int(line_split[0])][int(line_split[1])][genelist[line_split[2]]] = abs(float(line_split[3]))
    except IOError:
        print("Error reading file!")

    return matrix_dict

#function to build list of lists for heatmap from nested dictionary 
def heatmap_list(matrix_dict):
    list_data = []
    #if argument is gene name:
    if args.gene != "All":
        #create inner list for each y coord row
        for y in range(1, 37):
           list_data.append([])
           #for each inner list, add each x coordinates gene difference information one at a time from the nested dictionary
           for x in range(1, 126):
               list_data[-1].append(matrix_dict[y][x]["{0}".format(args.gene)]) # for a specific gene

    #if argument is "All"
    else:
        temp_list =[]
        #same as above
        for y in range(1, 37):
            list_data.append([])
            #same idea as above
            for x in range(1, 126):
                #for each x coordinate, take every gene difference value at that coordinate and take euclidean disance, and add that to inner list
                for g in matrix_dict[y][x]:
                     temp_list.append(matrix_dict[y][x][g] ** 2)
                     if len(temp_list) == len(matrix_dict[y][x]):
                         added = sum(temp_list)
                         num = cmath.sqrt(added)
                         list_data[-1].append(round(num.real, 2))
                         temp_list = []

    return list_data

#visualization of heatmap function
def getheatmap(list_data):
    #turn lists of lists into an array
    heat_map = np.array(list_data)
    fig, ax = plt.subplots() 
    #name plot accordingly 
    if args.gene == "All":
        plot_name = "Gene Expression Difference for All Genes"
    else:
        plot_name = "Gene Expression Difference for {0}".format(args.gene)
    #set axis names
    ax.set_title(plot_name)
    ax.set_xlabel('X Coordinates')
    ax.set_ylabel('Y Coordinates')
    #set colormap
    cmap = "YlGnBu"
    heatmap = ax.imshow(heat_map, cmap=cmap)
    fig.colorbar(heatmap, shrink=0.6, extend='both')
    fig = plt.show()



    return fig



#argparse arguments
parser = argparse.ArgumentParser()
parser.add_argument("gene",  help="You can specify a specific gene for heatmap generation, or All which will give you the average of all genes in the heatmap")
parser.add_argument("Features",  help="filepath to the features.tsv document outputted by Space Ranger Pipeline")
parser.add_argument("CSV", help="filepath to the output CSV file from GeneCompare.py") 
args = parser.parse_args()
start = time.time()
#run programs
#getheatmap(heatmap_list(dataget(genelookup())))
print(len(genelookup()))
end = time.time()
#prints time it took to ran
print((end-start)/60)

