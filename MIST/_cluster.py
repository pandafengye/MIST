import os
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings('ignore')
def cluster(refdir,output,cutoff):
 #   print(cutoff)
    if refdir==None:
        print("refdir is null, please enter the file path of the sequence!")
    # print("######################################   Cluster  #######################################")
#    if not os.path.exists(fastapath):
    file_list=os.listdir(refdir)
    fileObject = open(output+'/_MIST_cluster_sample_List', 'w')
    for i in file_list:
        fileObject.write(i)
        fileObject.write('\n')
    fileObject.close()
    pwd=os.getcwd()+'/'
#    print(pwd)
    os.chdir(refdir)
    print("fastANI --rl "+pwd+output+"/_MIST_cluster_sample_List --ql "+pwd+output+"/_MIST_cluster_sample_List -o " +pwd+output+"/_MIST_cluster_sample_Distance")
    os.system("fastANI --rl "+pwd+output+"/_MIST_cluster_sample_List --ql " +pwd+output+"/_MIST_cluster_sample_List -o " +pwd+output+"/_MIST_cluster_sample_Distance")
    CRO={};point=[]
    cutoff=cutoff.split(",")
    for co in cutoff:
#        print(co)
        if float(co)<=0 or float(co)>=1:
            print("Please enter the correct cutoff value,0<cutoff<1")
        else:
            Co=100*float(co)
            #print("rm -rf " +pwd+output+"/_Alpha_Cluster_sample_Distance_Filter")
            os.system("rm -rf " +pwd+output+"/_MIST_cluster_sample_Distance_Filter")
            os.system("cat "+pwd+output+"/_MIST_cluster_sample_Distance | awk -F'\t' '$3>"+str(Co)+"{print}' >> "+pwd+output+"/_MIST_cluster_sample_Distance_Filter")
            data = pd.read_table(pwd+output+"/_MIST_cluster_sample_Distance_Filter", header=None, usecols=[0,1])
            b = []
            N = data.as_matrix()  # 转换成矩阵
            for i in range(len(N)):
                b.append(N[i][0])
                b.append(N[i][1])
            point = list(set(b))  # 去重
            G = nx.Graph()#创建空图
            G.add_nodes_from(point)#加入结点
            edglist = []
            for i in range(len(N)):
                for j in range(1,len(N[0])):
                    edglist.append((N[i][0], N[i][j]))
            G = nx.Graph(edglist)
        # nx.spring_layout(G)
            position = nx.spring_layout(G)
            c=nx.connected_components(G)
            list1=[];list2=[]
            for i,z in enumerate(c):
                for a in z:
                    list1.append(i);list2.append(a)
            List=zip(list2,list1)
            CRO[co]=dict(List)
    df = pd.DataFrame.from_dict(CRO, orient='index')
    print(df.T)
    os.system("rm -rf " +pwd+output+"/_MIST_cluster_sample_Distance_Filter")
    os.system("rm -rf " +pwd+output+"/_MIST_cluster_sample_Distance")
    os.system("rm -rf " +pwd+output+"/_MIST_cluster_sample_List")
#    os.chdir(output)
    df.T.to_csv(pwd+output+"/_MIST_ref_cluster.csv")
