#! /usr/bin/env python3
import csv, time, os
import slideNormalization as normalize
import pandas as pd

start = time.time()

#For the barcode dictionary
def makeBcDict(bc_file1, bc_file2):
    print ("Making Barcode Dictionary")
    bcDict = {}
    ds1 = pd.read_csv(bc_file1)
    ds2 = pd.read_csv(bc_file2)
    final_ds1, final_ds2 = normalize.findMatch(ds1, ds2)

    #adds barcodes column to normalized data frame
    bcList = ds1['barcodes'].tolist()
    final_ds1.insert(0, 'barcodes', bcList)

 
    for y_coord in range(38):
        for x_coord in range(128):
            oldBarcode = str(ds1[(ds1.x == x_coord) & (ds1.y == y_coord)]['barcodes'].values)
            newBarcode = str(final_ds1[(final_ds1.x == (x_coord+1)) & (final_ds1.y == y_coord)]['barcodes'].values)

    #Returns one dictionary assuming the barcodes are the same for both visium slides after theyve been matched.
    if oldBarcode != newBarcode:
        bcDict.setdefault(oldBarcode, newBarcode)
    bcDict[''] = ''
    return bcDict

#testDict = makeBcDict('12hr_bc.csv', '6weeks_barcodes.csv')

#Barcodes must be sorted and matching in both the matrix file and bc file
def bcSort(bc_file, matrix_file, bc_dict):
    print("Sorting Matrix File")
    bc_header = ['']
    new_bc_header = []

    with open(bc_file,'r') as barcode_file:
        bc_reader = csv.reader(barcode_file)
        next(bc_reader)
        #Creates a list of barcodes in the original order of barcode file
        for row in bc_reader:
            bc_header.append(row[0])
        #creates an ordered list of barcodes but replaces barcodes based on normalized slide
        for bc in bc_header:
            if bc in bc_dict:
                new_bc_header.append(bc_dict.get(bc))
            else:
                new_bc_header.append(bc)

    with open(matrix_file, 'r', newline='') as infile, open('temp.csv', 'w', newline='') as tempfile:
        # output dict needs a list for new column ordering
        writer_temp = csv.DictWriter(tempfile, fieldnames=bc_header)
        # reorder the header first
        writer_temp.writeheader()
        for row in csv.DictReader(infile):
            # writes the reordered rows to the new file
            writer_temp.writerow(row)

    with open(matrix_file[:-4]+'_Reordered.csv', 'w', newline='') as outfile, open('temp.csv', 'r+', newline='') as tempfile:
        writer_out = csv.writer(outfile)
        #outputs matrix file in the same order as barcode file but with barcodes from the normalized slide
        writer_out.writerow(new_bc_header)
        reader_out = csv.reader(tempfile)
        next(reader_out)
        for row in reader_out:
            writer_out.writerow(row)

    os.remove('temp.csv')

#bcSort('12hr_bc.csv','12hr_Raw_Matrix.csv', testDict)
#bcSort('6weeks_barcodes.csv', '6weeks_Raw_Matrix.csv', testDict)

#For the matrix dictionary
def makeGeneDict(matrix_file):
    print("Making Gene Dictionary")
    y_coord = 0
    x_coord = 0
    gene_tuple = []
    gene_list = []
    gene_dict = {}
    inner_dict = {}
    outer_dict = {}

    with open(matrix_file, "r", newline="") as m_file:
        nreader = csv.reader(m_file)
        next(nreader)
        #Loop makes a list of genes and expression values
        for row in (nreader):
            for i in row:
                if "E" in i:
                    continue
                else:
                    gene_tuple.append((row[0], i))
            gene_list.append(gene_tuple)
            gene_tuple = []

        #nested for loop to make the nested dictionary
        for expres_vals in range(len(gene_list[0])):
            for gene_row in range(len(gene_list)):
                gene_dict[gene_list[gene_row][expres_vals][0]] = gene_list[gene_row][expres_vals][1]
            inner_dict[x_coord] = dict(gene_dict)
            x_coord += 1
            gene_dict = {}
            outer_dict[y_coord] = dict(inner_dict)
            if x_coord >   127:
                inner_dict = {}
                y_coord += 1
                x_coord = 0

        return outer_dict

#12hrDict = makeGeneDict("12hr_Raw_Matrix_Reordered.csv")
#6wkDict = makeGeneDict("6weeks_Raw_Matrix_Reordered.csv")

#compares two matrix dictionaries. Outputs csv with difference value for each gene at every x,y coordinate
def GeneCompare(geneDict1, geneDict2):
    print("Running Gene Comparison")
    n = 3
    hm_rows = ["Y Coordinates", "X Coordinates", "Gene", "Expression Difference Value"]

    # # Produces a new matrix with added expression values based on sliding window of size n
    with open("hm_comparisonB.csv", "w+", newline="") as new_hm_file:
        matrix_writer = csv.writer(new_hm_file, delimiter=",")
        matrix_writer.writerow(hm_rows)
        #For every gene:
        for i in geneDict1[0][0]:
            for outY in range(n//2,(39-n//2)):
                for outX in range(n//2,(128-n//2)): 
                    hm_rows = []
                    hm_rows.append(outY)
                    hm_rows.append(outX)
                    hm_rows.append(i)
                    super_val = 0
                    super_val2 = 0
                    z = 0
                    #for values in window size n
                    for subY in range(n):
                        for subX in range(n):
                            z += 1
                            superX = (outX+subX-n//2)
                            superY = (outY+subY-n//2)
                            #if at edges of coordinates, continue
                            if superX > 127:
                                continue
                            if superY > 38:
                                continue
                            super_val += int(geneDict1[superY][superX][i])
                            super_val2 += int(geneDict2[superY][superX][i])
                            dif_val = (super_val2 - super_val)
                            #Writes row when z only when n is squared b/c this is when there is a full sub range
                            if z == n**2:
                                hm_rows.append(dif_val)
                                matrix_writer.writerow(hm_rows)

def GridTransform(hm_csv):
    first_line = ["Y Coordinates","X Coordinates","Gene","Expression Difference Value"]
    y_val_odd = 1
    y_val_even = 2
    y_val = 1
    x_val = 1

    with open(hm_csv,'r') as in_csv, open(hm_csv[:-4]+'_Visium.csv','w') as out_csv:
        write = csv.writer(out_csv)
        write.writerow(first_line)

        read = csv.reader(in_csv)
        next(read)
        #Converts csv file with coordinates into the original shape of the visium alternating grid
        for row in in_csv:
            row = row.rstrip('\n')
            row = row.split(',')
            if x_val > 126:
                x_val = 1
                y_val_even += 2
                y_val_odd += 2
            if y_val_even > 74:
                y_val_even = 2
                y_val_odd = 1

            if x_val % 2 == 0:
                y_val = y_val_even
            else:
                y_val = y_val_odd
            ls = []
            ls.append(y_val)
            ls.append(x_val)
            ls.append(row[2])
            ls.append(row[3])
            write.writerow(ls)
            x_val += 1



GeneCompare(12hrDict, 6wkDict)
#GridTransform("hm_comparison.csv")
end = time.time()
print ('time:', (end - start)/60)

