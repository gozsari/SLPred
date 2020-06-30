import numpy as np
def read_iFeatures():
    file_feat = open("iFeatures_list", "r")
    set_iFeatures = set()
    for line in file_feat:
        set_iFeatures.add(line.strip())
    return set_iFeatures

def read_points_file(filename):
    pts = []
    prot_id_list = list()
    with open(filename, "r") as file_r:
        for line in file_r:
            if line[0] == "#":
                continue
            list_line = line.strip("\n").split()
            prot_id = list_line[0].split("|")[1].strip()
            prot_id_list.append(prot_id)
            pt = list_line[1:]
            #print(pt)
            ls = [float(value) for value in pt]
            pts.append(ls)
    return prot_id_list, pts

def read_data(directory,pos_class, coverage, FVTYPE):
    pos_prot_id_list, pts_0 = read_points_file(directory + pos_class +"_iFeatures_"+ coverage+"/"
                                + pos_class +"_"+coverage +"_positive_" + FVTYPE + ".fv")

    neg_prot_id_list, pts_1 = read_points_file(directory + pos_class +"_iFeatures_"+ coverage+"/"
                                + pos_class +"_"+coverage +"_negative_" + FVTYPE + ".fv")

    x = pts_0 + pts_1
    labels = [0] * len(pts_0) + [1] * len(pts_1)
    x = np.array(x)
    print(x.shape)
    return pos_prot_id_list, neg_prot_id_list, x, labels

def read_points_spmap(filename):
    pts = []
    prot_id_list = list()
    with open(filename, "r") as file_r:
        for line in file_r:
            list_line = line.strip("\n").split("\t")
            prot_id = list_line[0][1:].strip()
            prot_id_list.append(prot_id)
            pt = list_line[1:]
            #print(pt)
            ls = [float(value) for value in pt]
            pts.append(ls)
    return prot_id_list, pts
def read_spmap_features(directory,pos_class,coverage, FVTYPE):
    pos_prot_id_list, pts_0 = read_points_spmap(directory +"SPMAP/"
                                               + pos_class + "_" + coverage + "_positive_" + FVTYPE + ".fv")

    neg_prot_id_list, pts_1 = read_points_spmap(directory + "SPMAP/"
                                                    + pos_class + "_" + coverage + "_negative_" + FVTYPE + ".fv")

    x = pts_0 + pts_1
    labels = [0] * len(pts_0) + [1] * len(pts_1)
    x = np.array(x)
    print(x.shape)
    return pos_prot_id_list, neg_prot_id_list, x, labels

def read_points_pssm(filename):
    pts = []

    with open(filename, "r") as file_r:
        count = 0
        for line in file_r:
            if count == 0:
                count += 1
                continue
            list_line = line.strip("\n").split(",")

            pt = list_line
            ls = list()
            for value in pt:
                if value == "-inf":
                    value = "-9999999"
                elif value == "inf":
                    value = "9999999"

                ls.append(float(value))
            pts.append(ls)
    return pts
def read_pssm_features(directory,pos_class,coverage, FVTYPE):
    pts_0 = read_points_pssm(directory +pos_class+"_PSSM_features_"+coverage+"/"
                             + pos_class + "_" + coverage + "_positive_" + FVTYPE + ".csv")
    pts_1 = read_points_pssm(directory + pos_class+"_PSSM_features_"+coverage+"/"
                                 + pos_class + "_" + coverage + "_negative_" + FVTYPE + ".csv")
    x = pts_0 + pts_1
    labels = [0] * len(pts_0) + [1] * len(pts_1)
    x = np.array(x)
    print(x.shape)
    return x, labels