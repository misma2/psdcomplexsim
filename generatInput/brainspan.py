"""Ezt itt a Brainspan adatok kinyerésére vonatkozó script"""
import csv
from numpy import genfromtxt

def find_regions(_metadata_,_subject_):
    print(_metadata_)
    tmp=[]
    for data in _metadata_:
        if (data['donor_name']==_subject_):
            tmp.append(data['structure_acronym'])
    return tmp



def read_metadata(_filename_):
    tmp=[]
    with open(_filename_, mode='r') as csv_file:
        csv_reader = csv.DictReader(csv_file)
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                line_count += 1
            tmp.append(row)
            line_count += 1
        print(f'Processed {line_count} lines.')
    return tmp

def find_row_num(_csv_reader_,_gene_symbol_):
    for row in _csv_reader_:
        if row["gene_symbol"]==_gene_symbol_:
            tmp= row["row_num"]
    return tmp
def find_col_num(_csv_reader_,_donor_name_,_structure_acronym_):
    tmp="can't find"
    for row in _csv_reader_:
        if row["structure_acronym"]==_structure_acronym_ and row["donor_name"]==_donor_name_:
            tmp= row["column_num"]
    return tmp


def main():
    dir='D:\phd\genes_matrix_csv'
    my_data = genfromtxt(dir+'/expression_matrix.csv', delimiter=',')
    rows_metadata=read_metadata(dir+"/rows_metadata.csv")
    cols_metadata=read_metadata(dir+"/columns_metadata.csv")
    subjects=[]
    for col in cols_metadata:
        if col["donor_name"] not in subjects:
            subjects.append(col["donor_name"])
    for subject_ in subjects:
        subject=subject_
        brain_regions= find_regions(cols_metadata,subject) #["Ocx","M1C-S1C","AMY","MGE","STC","URL","CGE","DTH","MFC","DFC","OFC","LGE","ITC","HIP","VFC","PCx"]
        prots=["GRIN2B","GRIA1","DLG4","SYNGAP1","DLGAP1","SHANK3","HOMER1"]
        domains = {"GRIN2B": [["CNMDAR", 1]],
                   "GRIA1": [["CAMPA", 1]],
                   "DLG4": [["PDZtandem", 1], ["PDZ3", 1], ["SH3psd95", 1], ["GK", 1]],
                   "SYNGAP1": [["CSYNGAP", 1]],
                   "GKAP1": [["GBRS", 2], ["GH", 1], ["CGKAP", 1]],
                   "SHANK3": [["SH3", 1], ["PDZsh", 1], ["SAM", 1], ["virtualEVH1", 1]],
                   "HOMER1": [["EVH1", 1]]
                   }
        PPI = [["GRIN2B", "DLG4"], ["GRIA1", "DLG4"], ["DLG4", "SYNGAP1"], ["DLG4", "GKAP1"], ["GKAP1", "SHANK3"],
               ["SHANK3", "HOMER1"]]
        DDI = [["CNMDAR", "PDZtandem"], ["CAMPA", "PDZtandem"], ["CSYNGAP", "PDZ3"], ["GK", "GBRS"], ["CGKAP", "PDZsh"],
               ["virtualEVH1", "EVH1"]]


        import os
        directory = subject

        if not os.path.exists(directory):
            os.makedirs(directory)
        if False:
            fgo = open(subject+"/fgo.csv", "w")
            for prot in prots:
                fgo.write(prot+"\t7268\n")
            fgo.close()

            fcom = open(subject+"/fcom.csv", "w")
            for prot in prots:
                fcom.write(prot+"\tPSD\n")
            fcom.close()

            fdoms = open(subject+"/fdoms.csv", "w")
            for prot in prots:
                for domain in domains[prot]:
                    fdoms.write(prot+"\t"+domain+"\n")
            fdoms.close()

            fppi = open(subject+"/fppi.csv", "w")
            for interaction in PPI:
                fppi.write(interaction[0]+"\t"+interaction[1]+"\n")
            fppi.close()

            fddi = open(subject+"/fddi.csv", "w")
            for interaction in DDI:
                fddi.write(interaction[0]+"\t"+interaction[1]+"\n")
            fddi.close()
        import math
        for brain_region in brain_regions:
            region_num=find_col_num(cols_metadata,subject,brain_region)
            f = open(subject+"/abd_data_"+subject+"_"+brain_region+".csv", "w")
            for prot in prots:
                prot_row=find_row_num(rows_metadata,prot)
                print(brain_region)
                print(prot_row)
                print(region_num)
                print(my_data[int(prot_row)-1][int(region_num)])
                f.write(prot+"\t"+str(int(math.floor(my_data[int(prot_row)-1][int(region_num)]*300/66.82)))+"\n")
            f.close()
        """
        for brain_region in brain_regions:
            region_num=find_col_num(cols_metadata,subject,brain_region)
            #f = open(subject+"/abd_data_"+subject+"_"+brain_region+".csv", "w")
            szamok=[]
            for prot in ["DLGAP1","GKAP1"]:
                prot_row=find_row_num(rows_metadata,prot)
                print(brain_region)
                print(prot_row)
                print(region_num)
                szamok.append(my_data[int(prot_row)-1][int(region_num)])
                #f.write(prot+"\t"+str(int(my_data[int(prot_row)-1][int(region_num)]*100))+"\n")
            #f.close()
            print("kulonbseg")
            if(maxi<abs(szamok[0]-szamok[1])):
                maxi=szamok[0]-szamok[1]
        print(maxi)
        """


if __name__== "__main__":
    main()

