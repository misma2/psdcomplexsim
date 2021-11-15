import json
from tkinter import filedialog

import numpy as np

############Memory Error Handling#######################
#
# https://airbrake.io/blog/python/memoryerror
#


################My Codes###############


def makeDirectory(_dirName_):
    import os
    dirName = _dirName_
    try:
        # Create target Directory
        os.mkdir(dirName)
        print("Directory ", dirName, " Created ")
    except FileExistsError:
        print("Directory ", dirName, " already exists")


def complexMatrix(complex):
    size = len(complex["structure"])
    matrix = np.zeros((size, size))
    for binding in complex["bindings"]:
        matrix[binding["first_p"], binding["second_p"]] = 1
        matrix[binding["second_p"], binding["first_p"]] = 1
    return matrix


def columnSimilar(_column, column_):
    if (np.sum(_column) == np.sum(column_)):
        return True
    else:
        return False


def getLens(vectors):
    lens = []
    for vector in vectors:
        lens.append(len(vector))
    return lens.sort()


def find_indices(lst, condition):
    return [i for i, elem in enumerate(lst) if condition(elem)]


def getChains(complex, protein):
    protIndices = find_indices(complex["structure"], lambda e: e == protein)
    pairs = []
    for binding in complex["bindings"]:
        if (binding["first_p"] in protIndices and binding["second_p"] in protIndices):
            pairs.append([binding["first_p"], binding["second_p"]])
    pairs = np.array(pairs)
    prots = {}
    for pair in pairs:
        for prot in pair:
            if (prot in prots):
                prots[prot] += 1
            else:
                prots[prot] = 1
    terminals = {}
    for key in prots:
        if (prots[key] == 1):
            terminals[key] = True
    chains = []
    for key in terminals:
        pid = key
        chain = [pid]
        if (terminals[key]):
            neighbours = 2
            while neighbours != 1:
                for pair in pairs:
                    if (pid == pair[0] and pair[1] not in chain):
                        pid = pair[1]
                        chain.append(pid)
                        neighbours = prots[pid]
                        if (neighbours == 1):
                            terminals[pid] = False
                    if (pid == pair[1] and pair[0] not in chain):
                        pid = pair[0]
                        chain.append(pid)
                        neighbours = prots[pid]
                        if (neighbours == 1):
                            terminals[pid] = False
            chains.append(np.array(chain))
    chains = np.array(chains)
    """
    neighbourslists=[]
    for chain in chains:
        neighbourlist=[]
        for prot in chain:
            neighbourlist.append(prots[prot])
        neighbourslists.append(neighbourlist)
    """

    return chains


def getProtPosition(complex, chaintype, prottype):
    protIndices = find_indices(complex["structure"], lambda e: e == prottype)
    chains = getChains(complex, chaintype)
    protPositions = []
    for chain in chains:
        protPos = []
        for pos in chain:
            found = False
            for binding in complex["bindings"]:
                if (pos == binding["first_p"]):
                    if (binding["second_p"] in protIndices):
                        protPos.append(True)
                        found = True
                if (pos == binding["second_p"]):
                    if (binding["first_p"] in protIndices):
                        protPos.append(True)
                        found = True
            if (not found):
                protPos.append(False)
        protPositions.append(protPos)
    return protPositions


from itertools import permutations


def complexEqual(_complex, complex_):
    isthesame = False
    if (_complex["structure"] != complex_["structure"]):
        # print("Not the same")
        pass
    else:
        if ((complexMatrix(_complex) == complexMatrix(complex_)).all()):
            # print("The two complexes are the same. :) ")
            isthesame = True
        else:
            if (len(complex_["structure"]) >= 10):
                if (getLens(getChains(_complex, "SHANK3")) == getLens(getChains(complex_, "SHANK3"))):
                    # print("The two complexes may be the same")
                    return True
                else:
                    isthesame = False
                    # print("The complex is obviously not the same")
            else:
                # print("Length: "+str(len(complex_["structure"])))
                # print("May be the same with isomorph matrices")
                similarCols = []
                for i in range(len(complex_["structure"])):
                    for j in range(len(complex_["structure"])):
                        if (complex_["structure"][i] == complex_["structure"][j] and i != j):
                            # if(columnSimilar(complexMatrix(complex_)[:, i], complexMatrix(complex_)[:, j])):
                            similarCols.append(np.array([i, j]))
                similarCols = np.array(similarCols)[:int(len(similarCols) / 2), :]
                # print(complex_["structure"])
                if (len(similarCols >= 10)):
                    return True
                else:
                    # print("SimilarCols:"+str(len(similarCols)))
                    perms = list(permutations(range(0, len(similarCols))))
                    for perm in perms:
                        temp = complex_
                        for similarCol in perm:
                            col0 = similarCols[similarCol][0]
                            col1 = similarCols[similarCol][1]
                            for binding in temp["bindings"]:
                                if (binding["first_p"] == col0):
                                    binding["first_p"] = col1
                                else:
                                    if (binding["first_p"] == col1):
                                        binding["first_p"] = col0
                                if (binding["second_p"] == col0):
                                    binding["second_p"] = col1
                                else:
                                    if (binding["second_p"] == col1):
                                        binding["second_p"] = col0

                            if ((complexMatrix(_complex) == complexMatrix(temp)).all()):
                                isthesame = True
                                break
                            else:
                                temp = temp
                        else:
                            if (isthesame):
                                break
                            else:
                                continue
    return isthesame


def complexEqualNeighbors(_complex, complex_):
    if (_complex["structure"] == complex_["structure"]):
        _nodes = []
        nodes_ = []
        protID = 0
        for protein in _complex["structure"]:
            _node = {}
            _node["ID"] = protID
            _node["color"] = protein
            _node["neighbors"] = []
            for binding in _complex["bindings"]:
                if (binding["first_p"] == protID):
                    _node["neighbors"].append(_complex["structure"][binding["second_p"]])
                if (binding["second_p"] == protID):
                    _node["neighbors"].append(_complex["structure"][binding["first_p"]])
            _nodes.append(_node)
            node_ = {}
            node_["ID"] = protID
            node_["color"] = protein
            node_["neighbors"] = []
            for binding in complex_["bindings"]:
                if (binding["first_p"] == protID):
                    node_["neighbors"].append(complex_["structure"][binding["second_p"]])
                if (binding["second_p"] == protID):
                    node_["neighbors"].append(complex_["structure"][binding["first_p"]])
            nodes_.append(node_)
            protID += 1
        assigned = []
        for _node in _nodes:
            for node_ in nodes_:
                isAssigned = False
                if (_node["neighbors"].sort() == node_["neighbors"].sort() and node_["ID"] not in assigned):
                    assigned.append(node_["ID"])
                    isAssigned = True
                    break
            if (not isAssigned):
                return False
        if (len(assigned) == len(_nodes)):
            return True


    else:
        return False


def generateMatrices(complex):
    matrices = []
    similarCols = []
    for i in range(len(complex["structure"])):
        for j in range(len(complex["structure"])):
            if (complex["structure"][i] == complex["structure"][j] and i != j):
                # if(columnSimilar(complexMatrix(complex_)[:, i], complexMatrix(complex_)[:, j])):
                similarCols.append(np.array([i, j]))
    similarCols = np.array(similarCols)[:int(len(similarCols) / 2), :]

    from itertools import permutations
    perms = list(permutations(range(0, len(similarCols))))
    # print(perms)
    for perm in perms:
        temp = complex
        for similarCol in perm:
            col0 = similarCols[similarCol][0]
            col1 = similarCols[similarCol][1]
            for binding in temp["bindings"]:
                if (binding["first_p"] == col0):
                    binding["first_p"] = col1
                else:
                    if (binding["first_p"] == col1):
                        binding["first_p"] = col0
                if (binding["second_p"] == col0):
                    binding["second_p"] = col1
                else:
                    if (binding["second_p"] == col1):
                        binding["second_p"] = col0

            matrices.append(complexMatrix(temp))
            temp = temp
    # complex["matrices"]=matrices


def newIDsofaSimulation(file):
    with open(file) as f:
        data = json.load(f)
    f.close()

    tempComplexes = data["complexes"]

    new_complexes = []
    newID = 0
    counter = 0

    for i in range(len(tempComplexes)):
        if (tempComplexes[i]["ID"] == -1):
            pass
        else:
            for j in range(len(tempComplexes)):
                if (tempComplexes[j]["ID"] != -1):
                    if (tempComplexes[i]["ID"] != tempComplexes[j]["ID"]):
                        if (complexEqual2(tempComplexes[i], tempComplexes[j])):
                            counter += 1
                            tempComplexes[i]["occurence"][0]["abundance"] += tempComplexes[j]["occurence"][0][
                                "abundance"]
                            print("Warning:complex with ID: " + str(
                                tempComplexes[i]["ID"]) + " is the same as the complex with ID:" + str(
                                tempComplexes[j]["ID"]))
                            tempComplexes[j]["ID"] = -1

            new_complexes.append(tempComplexes[i])
            tempComplexes[i]["ID"] = -1
    del tempComplexes
    newIDcounter = 0
    for new_complex in new_complexes:
        new_complex["ID"] = int(newIDcounter)
        new_complex["occurence"][0]["abundance"] = int(new_complex["occurence"][0]["abundance"])
        for binding in new_complex["bindings"]:
            binding["first_p"] = int(binding["first_p"])
            binding["second_p"] = int(binding["second_p"])
        newIDcounter += 1

    print("Finished checking same complexes in the same simulation")
    return new_complexes


def plotChains(chains, GKAPs, chainprotname, GKAPname, chaincolor, GKAPcolor):
    # chains=getChains(data["complexes"][3],"SHANK3")
    # GKAPs=getProtPosition(data["complexes"][3],"SHANK3","GKAP1")
    # print(chains)
    x = []
    y = []
    for i in range(len(chains)):
        for j in range(chains[i].size):
            x.append(j)
            y.append(i)

    xx = []
    yy = []
    for i in range(len(GKAPs)):
        for j in range(len(GKAPs[i])):
            if (GKAPs[i][j]):
                xx.append(j)
                yy.append(i)

    import matplotlib.pyplot as plt
    shank = plt.scatter(x, y, color=chaincolor)
    gkap = plt.scatter(xx, yy, color=GKAPcolor)
    plt.legend((shank, gkap), (chainprotname + " protein alone", chainprotname + " protein with " + GKAPname),
               loc='upper right',
               ncol=3,
               fontsize=8)
    plt.title("Distribution of " + GKAPname + " on " + chainprotname + " chains in the complex:XX")
    # plt.grid(color = 'green', linestyle = '--', linewidth = 0.5)

    # Show the major grid lines with dark grey lines
    plt.grid(b=True, which='major', color='#666666', linestyle='-')

    # Show the minor grid lines with very faint and almost transparent grey lines
    plt.minorticks_on()
    plt.grid(b=True, which='minor', color='#999999', linestyle='--', alpha=0.8)

    plt.ylabel("Different chains")
    plt.xlabel("Position in a chain")
    plt.show()


import glob


def saveNewcomplexesFiles(save, folder):
    print("saving New Complexes Files started")
    print(folder)
    for i in range(5):
        i += 1
        if (save):
            files = sorted(glob.glob(folder + str(i) + "/*"))
        else:
            files = sorted(glob.glob(folder + "Averaged/*.json"))  # [100:]
        for file in files:
            newComplexes = newIDsofaSimulation(file)
            newData = {"complexes": newComplexes}
            if (save):
                with open(folder + "newIDaSim/" + str(i) + "/" + file.split('\\')[1], 'w') as outfile:
                    outfile.write(json.dumps(newData, indent=4, sort_keys=True))
                outfile.close()
    print("saving New Complexes Files finished")


import math


def calculateAverages(folder_selected):
    print("calculateAverages started")
    files = sorted(glob.glob(folder_selected + "newIDaSim/" + str(1) + "/*"))  # [100:]

    for file in files:
        with open(file) as f:
            _data = json.load(f)
        f.close()
        for c in _data["complexes"]:
            c["eredetibolvan"] = int(c["ID"])
        print("length of data: " + str(len(_data["complexes"])) + " in file: " + file)
        for i in range(2, 5):
            # print(file+"   "+str(i))
            with open(folder_selected + "newIDaSim/" + str(i) + "/" + file.split('\\')[1]) as f:
                data_ = json.load(f)
            f.close()
            for complex_ in data_["complexes"]:
                found = False
                for _complex in _data["complexes"]:
                    if (complexEqual2(_complex, complex_)):
                        _complex["occurence"][0]["abundance"] += complex_["occurence"][0]["abundance"]
                        found = True
                        # print(str(complex_["ID"])+" Were found")
                        break
                if (not found):
                    # print(str(complex_["ID"])+" Were not found")
                    complex_["appendelt"] = 1
                    _data["complexes"].append(complex_)
        del data_
        IDcounter = 0
        for complex in _data["complexes"]:
            complex["newID"] = int(IDcounter)
            IDcounter += 1
            complex["occurence"][0]["abundance"] = int(math.floor(complex["occurence"][0]["abundance"] / 5))
        with open(folder_selected + "/Averaged/" + file.split('\\')[1], 'w') as outfile:
            outfile.write(json.dumps(_data, indent=4, sort_keys=True))
        outfile.close()
        del _data
        print("calculateAverages finished")


def CreateNum(_complex_):
    prots = ["GRIN2B", "GRIA1", "DLG4", "SYNGAP1", "DLGAP1", "SHANK3", "HOMER1"]
    import numpy as np
    tmp = np.zeros(len(prots))
    count = 0
    # print(_complex_["structure"])
    for prot in prots:
        for prot_ in _complex_["structure"]:
            if (prot == prot_):
                tmp[count] += 1
        count += 1
    _complex_["num"] = tmp.tolist()


def createAllComplexesVector(folder):
    files = sorted(glob.glob(folder))
    allComplex = []
    for file in files:
        with open(file) as f:
            _data = json.load(f)
        f.close()
        for _complex in _data["complexes"]:
            found = False
            for complex_ in allComplex:
                if (complexEqual2(_complex, complex_)):
                    found = True
                    break
            if (not found):
                allComplex.append(_complex)
    primaryComplexes = []
    secondaryComplexes = []
    for complex in allComplex:
        CreateNum(complex)
        num = np.array(complex["num"])
        if (num.max() < 3):
            primaryComplexes.append(complex)
        else:
            secondaryComplexes.append(complex)

    primaryComplexesSorted = sorted(primaryComplexes, key=lambda k: k['num'])
    secondaryComplexesSorted = sorted(secondaryComplexes, key=lambda k: k['num'])
    allComplexSorted = primaryComplexesSorted + secondaryComplexesSorted
    newID = 0
    for complex in allComplexSorted:
        complex["ID"] = newID
        complex["occurence"][0]["abundance"] = int(0)
        newID += 1
    return allComplexSorted


def IDasjustion(folder_selected):
    print("ID adjustion started")
    folder = folder_selected + "Averaged/*.json"
    allComplexes = createAllComplexesVector(folder)
    files = sorted(glob.glob(folder))
    for file in files:
        _complexes = allComplexes
        with open(file) as f:
            data_ = json.load(f)
        f.close()
        for _complex in _complexes:
            for complex_ in data_["complexes"]:
                if (complexEqual2(_complex, complex_)):
                    _complex["occurence"][0]["abundance"] = complex_["occurence"][0]["abundance"]
        print(str(len(data_["complexes"])))
        data_["complexes"] = _complexes
        print("utanna: " + str(len(data_["complexes"])))
        with open(folder_selected + "Averaged/AveragedSorted/" + file.split('\\')[1], 'w') as outfile:
            outfile.write(json.dumps(data_, indent=4, sort_keys=True))
        outfile.close()
    print("ID adjustion finished")


def getLeaves(complex):
    psd95s = find_indices(complex["structure"], lambda e: e == "DLG4")
    psd95sNeighbours = []
    for psd95 in psd95s:
        leaves = []
        psd95N = {}
        for binding in complex["bindings"]:
            if (binding["first_p"] == psd95):
                leaves.append(complex["structure"][binding["second_p"]])
            if (binding["second_p"] == psd95):
                leaves.append(complex["structure"][binding["first_p"]])
        psd95N["leaves"] = leaves
        psd95N["ID"] = psd95
        psd95sNeighbours.append(psd95N)
    return psd95sNeighbours


def getPSD95s(complex):
    protIndices = find_indices(complex["structure"], lambda e: e == "DLGAP1")
    psd95Indices = find_indices(complex["structure"], lambda e: e == "DLG4")
    psd95Ns = getLeaves(complex)
    GKAPs = []
    for prot in protIndices:
        psd95s = []
        GKAP = {}
        for binding in complex["bindings"]:
            if (binding["first_p"] == prot and binding["second_p"] in psd95Indices):
                psd95s.append(complex["structure"][binding["second_p"]])
            if (binding["second_p"] == prot and binding["first_p"] in psd95Indices):
                psd95s.append(complex["structure"][binding["first_p"]])
        GKAP["psd95s"] = psd95s
        GKAP["ID"] = prot
        GKAPs.append(GKAP)
    for GKAP in GKAPs:
        leaves = []
        for psd95 in GKAP["psd95s"]:
            for psd95_ in psd95Ns:
                if (psd95 == psd95_["ID"]):
                    leaves.append(psd95_["leaves"])
        GKAP["leaves"] = leaves
    return GKAPs


def assignLeaves(_complex, complex_):  # Ez valójában a neighboursal megvolt PSD-95 szinten, de GKAP szinten nem
    _leaves = getPSD95s(_complex)
    leaves_ = getPSD95s(complex_)
    assigned = []
    for _leaf in _leaves:
        isAssigned = False
        for leaf_ in leaves_:
            if (_leaf["leaves"].sort() == leaf_["leaves"].sort() and leaf_["ID"] not in assigned):
                assigned.append(leaf_["ID"])
                isAssigned = True
                break
        if (not isAssigned):
            return False
    return True


def assignPositions(_complex, complex_, chainttype, prottype):
    _pos = getProtPosition(_complex, chainttype, prottype)
    pos_ = getProtPosition(complex_, chainttype, prottype)
    for _chain in _pos:
        found = False
        for chain_ in pos_:
            if ((_chain == chain_ or _chain == chain_[::-1]) and chain_[0] != "assigned"):
                found = True
                chain_[0] = "assigned"
        if (not found):
            return False
    return True


def complexEqual2(_complex, complex_):
    if (complexEqualNeighbors(_complex, complex_)):
        if (assignLeaves(_complex, complex_)):
            if (getLens(getChains(_complex, "SHANK3")) == getLens(getChains(complex_, "SHANK3"))):
                if (assignPositions(_complex, complex_, "SHANK3", "DLGAP1")):
                    if (assignPositions(_complex, complex_, "SHANK3", "HOMER1")):
                        return True
                    else:
                        return False
                else:
                    return False
            else:
                return False
        else:
            return False
    else:
        return False


def createColors(region_names):
    ##Generate color for each Region Type##
    region_types = []
    for region in region_names:
        region_types.append(region.split("\\")[1].split("_")[3][:-4])
    print(region_types)

    region_types_set = sorted(list(set(region_types)))
    print(region_types_set)

    from random import randint
    colors = []
    n = len(region_types_set)
    for i in range(n):
        colors.append('#%06X' % randint(0, 0xFFFFFF))
    return [colors, region_types_set, region_names, region_types]


def calcPCA(input):
    from numpy import genfromtxt
    with open(input, 'r', encoding='utf-8-sig') as f:
        inputData = genfromtxt(f, dtype=float, delimiter=';')

    ##Do the PCA##
    from sklearn.decomposition import PCA

    #X = np.transpose(np.array(np.delete(inputData, 0, 1)))
    X=np.transpose(inputData)
    print(X)

    pca = PCA(n_components=2)
    pca.fit(X)

    X_pca = pca.transform(X)
    print("original shape:   ", X.shape)
    print("transformed shape:", X_pca.shape)
    print("PCA components:")
    print(np.argmax(pca.components_,axis=1))
    print("variance: ", pca.explained_variance_ratio_)
    vmi=pca.components_
    vmi[0][np.argmax(pca.components_,axis=1)[0]]=0
    vmi[1][np.argmax(pca.components_, axis=1)[1]] = 0
    print(np.argmax(vmi, axis=1))

    np.savetxt(input+"PCA.csv", X_pca, delimiter=",")

    # X_new = pca.inverse_transform(X_pca)
    # plt.scatter(X[:, 0], X[:, 1], alpha=0.2)
    return X_pca


def plotPCA(folder_selected, input, output, colors, region_types_set, region_names, region_types):
    import matplotlib.pyplot as plt
    fig = plt.figure()
    # fig, (ax1, ax2,ax3) = plt.subplots(1, 2)
    afont = {'fontname': 'Arial'}
    # fig.suptitle('Principal Component Analyses',**afont)
    ax1 = plt.subplot(221)
    ax2 = plt.subplot(222)
    ax3 = plt.subplot(212)
    # ax1.plot(x, y)
    # ax2.plot(x, -y)
    print("inputunk:")
    print(input)
    for i in range(len(region_names)):
        ax1.scatter(input[i, 0], region_types_set.index(region_types[i]), color=colors[region_types_set.index(region_types[i])])#input[i, 1]
        ax2.scatter(output[i, 0], region_types_set.index(region_types[i]), color=colors[region_types_set.index(region_types[i])])#output[i, 1]

    import matplotlib.patches as mpatches
    pops = []
    for i in range(len(colors)):
        # legend
        pops.append(mpatches.Patch(color=colors[i], label=region_types_set[i]))
    # plt.axis('equal');
    # Shrink current axis's height by 10% on the bottom
    """
    box = ax2.get_position()
    ax2.set_position([box.x0, box.y0 + box.height * 0.1,
                     box.width, box.height * 0.9])
    box = ax1.get_position()
    ax1.set_position([box.x0, box.y0 + box.height * 0.1,
                      box.width, box.height * 0.9])           
    """

    ax1.set_xlabel("Principal Component 1", **afont)
    ax1.set_ylabel("Principal Component 2", **afont)
    ax2.set_xlabel("Principal Component 1", **afont)
    ax2.set_ylabel("Principal Component 2", **afont)
    ax1.set_title("A", **afont)
    ax2.set_title("B", **afont)

    # fig.subplots_adjust(bottom=0, wspace=0.33)
    ax3.axis("off")
    import matplotlib.font_manager as font_manager

    font = font_manager.FontProperties(family='Arial',
                                       weight='normal',
                                       style='normal', size=12)
    plt.legend(handles=pops,  # bbox_to_anchor=(0, 0),
               loc='center', ncol=5, numpoints=1, borderaxespad=0.1, prop=font)  # facecolor="plum"
    fig.tight_layout()

    from PIL import Image
    from io import BytesIO
    # (1) save the image in memory in PNG format
    png1 = BytesIO()
    fig.savefig(png1, format='png')

    # (2) load this image into PIL
    png2 = Image.open(png1)
    (width, height) = png2.size
    pngresized = png2.resize((int(math.floor(width * 1.5)), int(math.floor(height * 1.5))), Image.ANTIALIAS)
    (width, height) = pngresized.size

    # (3) save as TIFF
    pngresized.save(folder_selected + 'Figures/Fig5_linearproba.tiff', dpi=(600, 600))
    png1.close()


def plotComplexIDs(folder_selected):
    prots = ["GRIN2B", "GRIA1", "DLG4", "SYNGAP1", "DLGAP1", "SHANK3", "HOMER1"]
    with open(folder_selected + "Averaged/AveragedSorted/abd_data_H376.IIA.50_AMY.json") as f:
        data = json.load(f)
    f.close()
    all_nums = []
    IDs = []
    for complex in data["complexes"]:
        nums = np.array(complex["num"])
        if (True):  # nums.max()<3):
            all_nums.append(nums)
            IDs.append(complex["ID"])

    all_nums = np.array(all_nums)

    np.savetxt(folder_selected + "CSVs/complexIDs.csv", all_nums, delimiter=";")

    import matplotlib.pyplot as plt
    from matplotlib import colors
    ################  Complex IDs  ##############################
    fig, axs = plt.subplots(1)
    plt.tight_layout()
    print(all_nums.max())
    print(IDs)
    bd = -0.5
    bounds = []
    for i in range(int(all_nums.max()) + 2):
        bounds.append(bd)
        bd += 1
    cmap = colors.ListedColormap(['white', 'navajowhite', 'yellow', 'orange', 'red', 'purple', 'pink', 'blue'])
    # bounds = [-0.5, 0.5, 1.5, 2.5, 3.5, 4.5,5.5,6.5,7.5]
    norm = colors.BoundaryNorm(bounds, cmap.N)

    cax1 = axs.matshow(np.transpose(all_nums), cmap=plt.cm.get_cmap('rainbow', all_nums.max() + 1),
                       interpolation='none', norm=norm)
    # axs.tick_params(axis='both', which='minor', labelsize=5)YlOrRd

    # axs.set_xticks(np.arange(len(all_nums)))
    axs.set_xticks(IDs)
    axs.tick_params(axis='x', which='both', labelsize=8, rotation=90)
    axs.tick_params(axis='x', which='major', grid_color="black", labelsize=8, rotation=90, grid_linewidth=1)

    temp = axs.xaxis.get_ticklabels()
    temp = list(set(temp) - set(temp[::2]))
    for label in temp:
        label.set_visible(False)
    temp = axs.get_xgridlines()
    temp = list(set(temp) - set(temp[::2]))
    for grid in temp:
        grid.set_color('white')

    axs.set_yticks(np.arange(len(prots)))
    prots2 = ["NMDAR", "AMPAR", "PSD-95", "SynGAP1", "GKAP", "Shank3", "Homer1"]
    axs.set_yticklabels(prots2)
    axs.set_aspect(4.8)
    axs.grid()
    fig.colorbar(cax1, orientation="vertical", ticks=range(math.floor(all_nums.max()) + 1), label='Protein Count',
                 shrink=.5)
    plt.tight_layout()

    # plt.title("Complexes", pad=8)
    # plt.show()

    from PIL import Image
    from io import BytesIO
    # (1) save the image in memory in PNG format
    png1 = BytesIO()
    fig.savefig(png1, format='png')

    # (2) load this image into PIL
    png2 = Image.open(png1)
    (width, height) = png2.size
    png3 = png2.resize((int(math.floor(width * 1.5)), int(math.floor(height * 1.5))), Image.ANTIALIAS)
    # (3) save as TIFF
    png3.save(folder_selected + 'Figures/Fig2.tiff', dpi=(600, 600))
    png1.close()
    print("plotComplexIDs were saved")


def saveComplexAbundanceMatrix(folder_selected):
    files = sorted(glob.glob(folder_selected + "Averaged/AveragedSorted/*.json"))
    complexAbundanceMatrix = []
    for file in files:
        with open(file) as f:
            data = json.load(f)
        f.close()
        abundancevector = []
        for complex in data["complexes"]:
            abundancevector.append(complex["occurence"][0]["abundance"])
        abundancevector = np.array(abundancevector)
        complexAbundanceMatrix.append(abundancevector)
    complexAbundanceMatrix = np.array(complexAbundanceMatrix)
    np.savetxt(folder_selected + "CSVs/output.csv", np.transpose(complexAbundanceMatrix), delimiter=';')


def calcDistances(input, output):
    from sklearn.metrics.pairwise import euclidean_distances
    from numpy import genfromtxt
    with open(input, 'r', encoding='utf-8-sig') as f:
        inputData = np.transpose(genfromtxt(f, dtype=int, delimiter=';'))
    with open(output, 'r', encoding='utf-8-sig') as f:
        outputData = np.transpose(genfromtxt(f, dtype=int, delimiter=';'))

    abd_distances = euclidean_distances(inputData, inputData)
    com_distances = euclidean_distances(outputData, outputData)

    print(abd_distances)
    print("254hez közeli")
    print(np.argsort(abd_distances[254]))
    print("339hez közeli")
    print(np.argsort(abd_distances[339]))
    abd_distances = (abd_distances - np.min(abd_distances)) / np.ptp(abd_distances)
    print(com_distances)
    com_distances = (com_distances - np.min(com_distances)) / np.ptp(com_distances)
    distances = abs(abd_distances - com_distances)

    result = np.where(distances == np.amax(distances))
    print(result)

    abd_distances[result[0][0]][result[0][0]] = 100
    result2 = np.where(abd_distances == np.amin(abd_distances[result[0][0]]))
    print(result2)
    # print(result[0][0])

    return distances


def plotDistances(distances, folder_selected):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt

    with open(folder_selected + "CSVs/M.csv", 'r', encoding='utf-8-sig') as f:
        clusterDistances = genfromtxt(f, dtype=float, delimiter=',')

    fig = plt.figure()
    plt.tight_layout()
    # fig, (ax1, ax2,ax3) = plt.subplots(1, 2)
    afont = {'fontname': 'Arial'}
    # fig.suptitle('Principal Component Analyses', **afont)
    ax1 = plt.subplot(121)

    cax1 = ax1.matshow(distances, cmap="rainbow")
    fig.colorbar(cax1, orientation="horizontal", label='Normalized distances', shrink=.5)
    ax2 = plt.subplot(122)

    ax1.set_xlabel("Brain Regions", **afont)
    ax1.set_ylabel("Brain regions", **afont)
    ax1.set_title("A", **afont)  # Distances of the distances between brain regions' inputs and outputs", pad=8)

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    sns.set_theme(style="white")

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(clusterDistances, cmap="rainbow", vmax=distances.max(),
                square=True, linewidths=.5, cbar_kws={"shrink": .5, "orientation": "horizontal"}, annot=True)
    # ax2.matshow(clusterDistances, cmap="rainbow")
    ax2.set_xlabel("Input Clusters", **afont)
    ax2.set_ylabel("Output Clusters", **afont)
    ax2.set_title("B", **afont)

    # plt.savefig("outputs/plots/" + simulation_name + "/heatmapdistances.png")
    # plt.show()
    """
    from string import ascii_letters

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    sns.set_theme(style="white")

    # Compute the correlation matrix
    #corr = d.corr()

    # Generate a mask for the upper triangle
    mask = np.triu(np.ones_like(distances, dtype=bool))

    # Set up the matplotlib figure
    f, ax = plt.subplots(figsize=(100, 100))

    # Generate a custom diverging colormap
    #cmap = sns.diverging_palette(230, 20, as_cmap=True)

    # Draw the heatmap with the mask and correct aspect ratio
    sns.heatmap(distances, mask=mask, cmap="rainbow", vmax=distances.max(),
                square=True, linewidths=.5, cbar_kws={"shrink": .1}) #mask=mask, center=0
    import matplotlib.pylab as pylab
    params = {'legend.fontsize': 100,
              'figure.figsize': (15, 5),
              'axes.labelsize': 'x-large',
              'axes.titlesize': 500,
              'xtick.labelsize': 'x-large',
              'ytick.labelsize': 'x-large'}
    pylab.rcParams.update(params)
    plt.title("Changes in Cross-Distances")
    plt.xlabel("Brain Regions")
    plt.ylabel("Brain Regions")
    #plt.show()
    #f.savefig("test_hidpi_tiff.tifff", dpi=500)
    # save figure
    """
    from PIL import Image
    from io import BytesIO
    # (1) save the image in memory in PNG format
    png1 = BytesIO()
    fig.savefig(png1, format='png')

    # (2) load this image into PIL
    png2 = Image.open(png1)
    (width, height) = png2.size
    png3 = png2.resize((int(math.floor(width * 1.5)), int(math.floor(height * 1.5))), Image.ANTIALIAS)
    # (3) save as TIFF
    png3.save(folder_selected + 'Figures/Fig4.tiff', dpi=(600, 600))
    png1.close()


def plotHeatmapAbundances(input, output, folder_selected):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    with open(input, 'r', encoding='utf-8-sig') as f:
        inputData = genfromtxt(f, dtype=float, delimiter=';')
    with open(output, 'r', encoding='utf-8-sig') as f:
        outputData = genfromtxt(f, dtype=int, delimiter=';')
    fig = plt.figure()
    plt.tight_layout()
    afont = {'fontname': 'Arial'}
    # fig.suptitle('Principal Component Analyses', **afont)
    ax1 = plt.subplot(211)
    cax1 = ax1.imshow(inputData, cmap="rainbow", aspect='auto')
    fig.colorbar(cax1, orientation="vertical", label='Protein Abundance', shrink=.5)

    prots = ["NMDAR", "AMPAR", "PSD-95", "SynGAP1", "GKAP", "Shank3", "Homer1"]
    ax1.set_yticks(np.arange(len(prots)))
    ax1.set_yticklabels(prots)

    ax1.set_xlabel("Brain regions", **afont)
    ax1.set_title("A", **afont)
    ax2 = plt.subplot(212)

    cax2 = ax2.imshow(outputData, cmap="rainbow", aspect='auto')
    fig.colorbar(cax2, orientation="vertical", label='Complex Abundance', shrink=.5)

    ax2.set_xlabel("Brain regions", **afont)
    ax2.set_ylabel("Complex IDs", **afont)
    ax2.set_title("B", **afont)

    from PIL import Image
    from io import BytesIO
    # (1) save the image in memory in PNG format
    png1 = BytesIO()
    fig.savefig(png1, format='png')

    # (2) load this image into PIL
    png2 = Image.open(png1)
    (width, height) = png2.size
    png3 = png2.resize((int(math.floor(width * 1.5)), int(math.floor(height * 1.5))), Image.ANTIALIAS)
    # (3) save as TIFF
    png3.save(folder_selected + 'Figures/Fig3.tiff', dpi=(600, 600))
    png1.close()


def plotPies(input, output, regions, elseboundary, colors2, folder_selected):
    import matplotlib.pyplot as plt
    from numpy import genfromtxt
    with open(input, 'r', encoding='utf-8-sig') as f:
        inputData = np.transpose(genfromtxt(f, dtype=float, delimiter=';'))
    with open(output, 'r', encoding='utf-8-sig') as f:
        outputData = np.transpose(genfromtxt(f, dtype=int, delimiter=';'))
    fig = plt.figure()
    afont = {'fontname': 'Arial'}

    labels = np.arange(7)
    # sizes = [15, 30, 45, 10]
    colors = {
        0: "blue",
        1: "orange",
        2: "green",
        3: "red",
        4: "purple",
        5: "brown",
        6: "pink"
    }

    ax1 = plt.subplot(231)
    ax1.pie(inputData[regions[0]],
            labels=labels,
            colors=[colors[key] for key in labels])
    ax2 = plt.subplot(232)
    ax2.pie(inputData[regions[1]],
            labels=labels,
            colors=[colors[key] for key in labels])
    ax3 = plt.subplot(233)
    ax3.pie(inputData[regions[2]],
            labels=labels,
            colors=[colors[key] for key in labels])
    ax1.set_title("A", **afont)
    ax2.set_title("B", **afont)
    ax3.set_title("C", **afont)


    ax4 = plt.subplot(234)

    abper = []
    abnev = []
    egyeb = 0
    count = 0

    for item in outputData[regions[0]]:
        if (item < elseboundary):
            egyeb += item
        else:
            abper.append(item)
            abnev.append(count)
        count += 1

    abper.append(egyeb)
    abnev.append("other")

    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = abnev
    sizes = abper
    print(sizes)
    explode = (0, 0, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
    ax4.pie(sizes, labels=labels, colors=[colors2[key] for key in labels],  # autopct='%1.1f%%',
            shadow=False, startangle=90)
    ax4.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax4.set_title("D")
    ax5 = plt.subplot(235)

    abper = []
    abnev = []
    egyeb = 0
    count = 0

    for item in outputData[regions[1]]:
        if (item < elseboundary):
            egyeb += item
        else:
            abper.append(item)
            abnev.append(count)
        count += 1

    abper.append(egyeb)
    abnev.append("other")

    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = abnev
    sizes = abper
    print(sizes)
    explode = (0, 0, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
    ax5.pie(sizes, labels=labels, colors=[colors2[key] for key in labels],  # autopct='%1.1f%%',
            shadow=False, startangle=90)
    ax5.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax5.set_title("E")

    ax6 = plt.subplot(236)

    abper = []
    abnev = []
    egyeb = 0
    count = 0

    for item in outputData[regions[2]]:
        if (item < elseboundary):
            egyeb += item
        else:
            abper.append(item)
            abnev.append(count)
        count += 1

    abper.append(egyeb)
    abnev.append("other")

    # Pie chart, where the slices will be ordered and plotted counter-clockwise:
    labels = abnev
    sizes = abper
    print(sizes)
    explode = (0, 0, 0, 0)  # only "explode" the 2nd slice (i.e. 'Hogs')
    ax6.pie(sizes, labels=labels, colors=[colors2[key] for key in labels],  # autopct='%1.1f%%',
            shadow=False, startangle=90)
    ax6.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
    ax6.set_title("F")

    plt.show()


def changeFig(figname, savefig):
    from PIL import Image
    from io import BytesIO

    # (2) load this image into PIL
    png2 = Image.open(figname)
    (width, height) = png2.size
    q=1560/width
    png3 = png2.resize((int(math.floor(width * q)), int(math.floor(height * q))), Image.ANTIALIAS)
    # (3) save as TIFF
    png3.save(savefig, dpi=(300, 300))
    png2.close()


def tSNE(input,output):
    import numpy as np
    from numpy import genfromtxt
    with open(input, 'r', encoding='utf-8-sig') as f:
        inputData = np.transpose(genfromtxt(f, dtype=float, delimiter=';'))
    with open(output, 'r', encoding='utf-8-sig') as f:
        outputData = np.transpose(genfromtxt(f, dtype=int, delimiter=';'))
    from sklearn.manifold import TSNE
    tSNEs=[]
    for matrix in [inputData,outputData]:
        X_embedded = TSNE(n_components=2).fit_transform(matrix)
        tSNEs.append(X_embedded)
    return tSNEs

def calcSizes(abundance, complexes,T,folder_selected):
    from numpy import genfromtxt
    with open(abundance, 'r', encoding='utf-8-sig') as f:
        abundanceData = np.transpose(genfromtxt(f, dtype=int, delimiter=';'))
    with open(complexes, 'r', encoding='utf-8-sig') as f:
        complexesData = np.transpose(genfromtxt(f, dtype=int, delimiter=';'))

    complexSizes = np.sum(complexesData, axis=0)
    complexCount = np.sum(abundanceData, axis=1)

    averageSizes = abundanceData.dot(complexSizes) * np.reciprocal(complexCount.astype(float))
    np.savetxt(folder_selected+"CSVs/averageSizes.csv",averageSizes,delimiter=';')
    #print(np.max(complexesData, axis=0))
    #Fact=str(np.max(complexesData, axis=0))
    #T.insert('1.0',Fact)

from tkinter import *


def runProcessing(folder_selected,root):

    makeDirectory(folder_selected + "Figures")
    makeDirectory(folder_selected + "CSVs")
    makeDirectory(folder_selected + "newIDaSim")
    for i in range(1, 6):
        makeDirectory(folder_selected + "newIDaSim/" + str(i))
    makeDirectory(folder_selected + "Averaged")
    makeDirectory(folder_selected + "Averaged/AveragedSorted")
    saveNewcomplexesFiles(True, folder_selected)
    calculateAverages(folder_selected)
    saveNewcomplexesFiles(False, folder_selected)
    IDasjustion(folder_selected)
    files = sorted(glob.glob(folder_selected + "Averaged/*.json"))

    plotComplexIDs(folder_selected)
    saveComplexAbundanceMatrix(folder_selected)
    plotDistances(calcDistances(folder_selected + "CSVs/input.csv", folder_selected + "CSVs/output.csv"),
                  folder_selected)


    files = sorted(glob.glob(folder_selected + "Averaged/*.json"))
    colorsandtypes = createColors(region_names=files)

    #[inputtSNE,outputtSNE]=tSNE(folder_selected+"CSVs/input.csv",folder_selected+'CSVs/output.csv')
    """
    plotPCA(folder_selected, inputtSNE,
            outputtSNE, colors=colorsandtypes[0], region_types_set=colorsandtypes[1],
            region_names=colorsandtypes[2], region_types=colorsandtypes[3])
    """
    inputRound="../inputFiles/generatInput/inputRound.csv" #calcPCA(folder_selected + "CSVs/output.csv")
    allRNAs="D:\phd\genes_matrix_csv\expression_matrix.csv"
    plotPCA(folder_selected, calcPCA(folder_selected+"CSVs/input.csv"),calcPCA(folder_selected+'CSVs/output.csv') , colors = colorsandtypes[0], region_types_set = colorsandtypes[1], region_names = colorsandtypes[2], region_types = colorsandtypes[3])
    #plotHeatmapAbundances("input.csv", folder_selected + "CSVs/output.csv", folder_selected)
    #calcSizes(folder_selected + "CSVs/output.csv", folder_selected + "CSVs/complexIDs.csv",root)
    """
    distances=calcDistances(folder_selected + "CSVs/input.csv", folder_selected + "CSVs/output.csv")
    from numpy import linalg as LA
    norms=LA.norm(distances,axis=0)
    print(np.argsort(norms))
    print(np.argmax
          (norms))
    """
    calcPCA(folder_selected + "CSVs/input.csv")
    calcPCA(folder_selected + "CSVs/output.csv")
    #kMeans(calcPCA(folder_selected+"CSVs/output.csv"),folder_selected)
    print("Finished")

def runPlotChains(selected_folder,cID):
    files=glob.glob(selected_folder+"Averaged/AveragedSorted/*.json")
    with open(files[0]) as f:
        data_ = json.load(f)
    f.close()

    plotChains(getChains(data_["complexes"][cID], "SHANK3"), getProtPosition(data_["complexes"][cID], "SHANK3", "DLGAP1"),
               "Shank3", "GKAP", "brown", "purple")


def kMeans(X, folder_selected):
    # Instantiate the KMeans models
    #
    from sklearn import datasets
    from sklearn.cluster import KMeans
    from sklearn.metrics import silhouette_score
    km = KMeans(n_clusters=3, random_state=42)
    #
    # Fit the KMeans model
    #
    km.fit_predict(X)
    #
    # Calculate Silhoutte Score
    #
    score = silhouette_score(X, km.labels_, metric='euclidean')
    #
    # Print the score
    #
    print('Silhouetter Score: %.3f' % score)

    from yellowbrick.cluster import SilhouetteVisualizer
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(2, 2, figsize=(15, 8))
    for i in [2, 3, 4, 5]:
        '''
        Create KMeans instance for different number of clusters
        https://dzone.com/articles/kmeans-silhouette-score-explained-with-python-exam
        '''
        km = KMeans(n_clusters=i, init='k-means++', n_init=10, max_iter=100, random_state=42)
        q, mod = divmod(i, 2)
        '''
        Create SilhouetteVisualizer instance with KMeans instance
        Fit the visualizer
        '''
        visualizer = SilhouetteVisualizer(km, colors='yellowbrick', ax=ax[q - 1][mod])
        visualizer.ax.set_xlabel("silhouette score")
        visualizer.ax.set_ylabel("simulations")
        visualizer.fit(X)
    from PIL import Image
    from io import BytesIO
    # (1) save the image in memory in PNG format
    png1 = BytesIO()
    fig.savefig(png1, format='png')

    # (2) load this image into PIL
    png2 = Image.open(png1)
    (width, height) = png2.size
    png3 = png2.resize((int(math.floor(width * 1.5)), int(math.floor(height * 1.5))), Image.ANTIALIAS)
    # (3) save as TIFF
    png3.save(folder_selected + 'Figures/S4.tiff', dpi=(600, 600))
    png1.close()
    #fig.show()

def runChangeFig():
    file = filedialog.askopenfilename()
    savefile = filedialog.asksaveasfilename()
    changeFig(file, savefile)


def main():
    # Memory Handling#
    try:
        root = Tk()
        root.geometry("500x250")
        frame = Frame(root)
        frame.pack()

        leftframe = Frame(root)
        leftframe.pack(side=LEFT)

        rightframe = Frame(root)
        rightframe.pack(side=RIGHT)

        label = Label(frame, text="CytoCast PostProcessing")
        label.pack()
        Output = Text(root, height=5,
                      width=25,
                      bg="light cyan")
        Output.pack(padx=1,pady=6)
        Output.insert('1.0',"Output")
        from tkinter import filedialog
        folder_selected = filedialog.askdirectory() + "/"
        button2 = Button(root, text="Run", command=lambda : runProcessing(folder_selected,Output), fg="red", font="Verdana 14",
                         bd=2, bg="light blue", relief="groove")
        button2.pack(padx=1, pady=1)

        L1 = Label(root, text="Plot Pie Charts of the regions with ID below")
        L1.pack(padx=1, pady=2)
        E1 = Entry(root, bd=4)
        E1.pack(padx=1, pady=3)
        E2 = Entry(root, bd=4)
        E2.pack(padx=1, pady=3)
        E3 = Entry(root, bd=4)
        E3.pack(padx=1, pady=3)

        from random import randint
        colors = {}
        for i in range(0,94):
            # colors.append()
            colors[i] = '#%06X' % randint(0, 0xFFFFFF)
        colors["other"] = "pink"

        button3 = Button(leftframe, text="Plot Pies",
                         command=lambda : plotPies("input.csv", folder_selected + "CSVs/output.csv",
                                                  [int(E1.get()), int(E2.get()), int(E3.get())], 2, colors,folder_selected))
        button3.pack(padx=3, pady=3)

        button1 = Button(leftframe, text="Change Figure", command=runChangeFig,
                         fg="red", font="Verdana 14",
                         bd=2, bg="light blue", relief="groove")
        button1.pack(padx=4, pady=4)

        L2 = Label(root, text="Complex ID to plot Chains")
        L2.pack(padx=5, pady=6)
        E4 = Entry(root, bd=4)
        E4.pack(padx=5, pady=6)
        buttonChain = Button(leftframe, text="Plot Chains", command=lambda: runPlotChains(folder_selected,int(E4.get())),
                         fg="red", font="Verdana 14",
                         bd=2, bg="light blue", relief="groove")
        buttonChain.pack(padx=5, pady=5)
        root.bind('<Escape>', lambda e:root.quit())
        root.title("CytoCast PostProcessing")
        root.mainloop()
        return 0
    except MemoryError as error:
        # Output expected MemoryErrors.
        log_exception(error)
    except Exception as exception:
        # Output unexpected Exceptions.
        log_exception(exception, False)


if __name__ == "__main__":
    main()
