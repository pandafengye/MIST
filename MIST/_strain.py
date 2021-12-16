import os,math,re,sys,datetime,heapq
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from numpy import nan
import matplotlib.pyplot as plt
from sklearn import datasets,ensemble,naive_bayes
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib
import warnings,datetime
sys.setrecursionlimit(100000) #
warnings.filterwarnings('ignore')
#from collections import Counter

def as_num(x):
    y = '{:.10f}'.format(x)  # .10f 
    return y

def getlistnum(li, num):  # This function is to count each element of the list
    li = list(li)
    set1 = set(range(num))
    dict1 = {}
    for item in set1:
        dict1.update({item: li.count(item)})
    return dict1


def list_add(a, b):
    c = []
    for i in range(len(a)):
        c.append(str(a[i]) + str(b[i]))
    return c

def strain(indexpath,read_length,cluster_output,output,threads=8,single_end=None,pair_1=None,pair_2=None,genome_size=None):
    if not os.path.exists(output + "/_MIST_strain/"):
        os.makedirs(output+ "/_MIST_strain/")
    if not os.path.exists(output + "/_MIST_strain/_MIST_map_alignment/"):
        os.makedirs(output+ "/_MIST_strain/_MIST_map_alignment/")
    """Map the fastq file to the reference genome using bowtie2"""

    index_list=os.listdir(indexpath)
    if pair_1==None and pair_2==None and single_end!=None:
        for i in index_list:
            os.system("bowtie2 --local  -p {} -x {} -U {} -S {}".format(threads,indexpath + i + "/" + i, single_end,
                                                                             output + "/_MIST_strain/_MIST_map_alignment/" + i + ".fq.sam"))
    if pair_1!=None and pair_2!=None and single_end==None:
        os.system("cat {} {} > {} ".format(pair_1, pair_2,output+"/_MIST_strain/_MIST_map_merge.fq"))
        for i in index_list:
            os.system("bowtie2 --local  -p {} -x {} -U {} -S {}".format(threads,indexpath + i + "/" + i, output+"/_MIST_strain/_MIST_map_merge.fq",
                                                                            output + "/_MIST_strain/_MIST_map_alignment/" + i + ".fq.sam"))
    os.system("rm -rf {}/_MIST_strain/_MIST_map_merge.fq".format(output))
    """Statistics Mismatch number"""
    DICT = {}
    alignment_filePath=output + "/_MIST_strain/_MIST_map_alignment/"
    files = os.listdir(alignment_filePath)  # Get all file names under the folder
    for file in files:  
        #print(file)
        if not os.path.isdir(file):  
            f = open(alignment_filePath + "/" + file, 'r')
            result = []
            key = []
            value = []
            reg = r'NM:i:(.{1})'
            wordreg = re.compile(reg)
            for line in f:
                # f.decode("utf8","ignore")
                li = line.strip()
                if not li.startswith("@"):  # Skip @ previous three lines
                    a = li.split()  
                    b = a[0:2];
                    c1 = a[5:6];  # Extract the first two and fifth columns
                    c2 = re.findall(wordreg, li)  # Find the value behind NM: i:
                    key.append(b);
                    value.append(c1 + c2)

            Mismatch_matrix = []  # Calculate Mismatch
            zimu_regex = re.compile(r'[a-zA-z]')
            for i, j in zip(value, key):

                mismatch_matrix = {}
                zimu_list = zimu_regex.findall(i[0])
                mismatch = 0
                if i[0] == "*":
                    mismatch = read_length
                else:
                    mismatch = int(i[1])
                mismatch_matrix[j[0]] = mismatch
                Mismatch_matrix.append(mismatch_matrix)
            #print(Mismatch_matrix)
            Dict0 = {}
            for i in Mismatch_matrix:
                if list(i.keys())[0] not in Dict0:
                    #     print(DICT[list(b.keys())[0]])
                    Dict0[list(i.keys())[0]] = list(i.values())[0]
                else:
                  #  print(list(i.values())[0])
                   # print(Dict0[list(i.keys())[0]])
                    Dict0[list(i.keys())[0]] = list(i.values())[0] + Dict0[list(i.keys())[0]]
            Dict0 = dict(sorted(Dict0.items(), key=lambda d: d[0]))
        DICT['.'.join(file.split(".")[0:2])] = Dict0
    import pandas as pd
    mismatch_result = pd.DataFrame.from_dict(DICT, orient='index')
    mismatch_result=mismatch_result.T
    Count=0
    for i in range(len(mismatch_result)):# Line by line
        if any(mismatch_result.iloc[i]<10):# Look for reads with mismatch less than 10, which is a complete comparison
            Count=Count+1# The number of the complete comparison
    print("There are a total of {} sequences, a total of {} mapped!".format(len(mismatch_result),Count))
    mismatch_result.to_csv(output+"/_MIST_strain/_MIST_map_Mismatch_matrix.csv")
    """Separate each cluster  """
    cluster=pd.read_csv(cluster_output,index_col=0)# read data
    Similarity_all=cluster.columns.values.tolist()
    #print(Similarity_all)
    if genome_size!=None:
        for i in range(cluster.shape[1]):
            print("Level {}, start building prob matrix\n".format(Similarity_all[i]))
    #         print(Similarity_all[i])
    #        start = datetime.datetime.now()
    ##############################################################################################
            dict_cluster=cluster.iloc[:,[i]].T.to_dict(orient='records')[0]
            for key in dict_cluster:
                dict_cluster[key] = str(dict_cluster[key])
            p_matrix=mismatch_result.copy()
            count = len(p_matrix.index)#;print(count)
            p_matrix.rename(columns=dict_cluster, inplace = True)
            List1=sorted(set(p_matrix.columns))#print(List1)
    #         List1=list(map(int,List1))
            P_matrix1=pd.DataFrame()
            for j in List1:
                if isinstance(p_matrix[j],pd.Series):
                    df=p_matrix[j].to_frame()
                if isinstance(p_matrix[j],pd.DataFrame):
                    df=p_matrix[j]
                P_matrix1[j]=df.min(axis=1)
            data1 = (0.05**P_matrix1.values)*((1-0.05)**(0-P_matrix1.values))
            P_matrix2 = pd.DataFrame(data1)
    #         print(P_matrix1.columns)
            P_matrix2.columns = P_matrix1.columns
            P_matrix2.insert(0, 'reads', P_matrix1.index)
        #        print(" {} Done\n ".format(Similarity_all[i]))
        #    print(P_matrix1.head())
            """Count the minimum number of Mismatch, the number of unique minimum values, 
               and the total Mismatch and total similarity of each cluster"""
            DT = {}
            data = P_matrix1
            Total_Min = {};
            Unique_Min = {};
            Total_Mismatch = {};
            Total_Similarity = {};
            Total_Depth={}
            values = data.min(axis=1)
            Keys = data.columns
            for k in Keys:
                Min = 0;
                Mismatch = 0;
                Similarity = 0
                temp = data[k]
                for j, z in zip(temp, values):
                    # print(j)
                    if j == z:  #
                        if z < read_length/20 :
                            Mismatch = Mismatch+z
                            Min = Min + 1
                            if pair_1!=None and pair_2!=None and single_end==None:
                                Similarity = 1 - Mismatch / (2*read_length * Min)
                            else:
                                Similarity = 1 - Mismatch / (read_length * Min)
                Total_Min[k] = Min;
                Total_Mismatch[k] = Mismatch
                Total_Similarity[k] = Similarity;
           
     #       print(Total_Similarity)
            #         print(Total_Min)
            for j in range(len(values)):
                unique_min = 0
                d = defaultdict(list)
                #         print(data1.iloc[0])
                for k, vy in [(v, y) for y, v in zip(data.columns, data.iloc[j])]:
                    d[k].append(vy)
                min_count = sorted(d.items(), key=lambda d: d[0], reverse=True)[0]
                #             print(min_count)
                if len(min_count[1]) == 1:
                    unique_min = unique_min + 1
                    if min_count[1][0] in Unique_Min.keys():
                        unique_min = Unique_Min[min_count[1][0]] + 1
                #             print(unique_min)
                Unique_Min[min_count[1][0]] = unique_min
            for k in Keys:
                if k not in Unique_Min.keys():
                    Unique_Min[k] = 0
            for k in Keys:
                Depth = read_length*Total_Min[k]/genome_size
                Total_Depth[k]=Depth
            DT["Total_Min"] = Total_Min
            DT["Unique_Min"] = Unique_Min
            DT["Total_Mismatch"] = Total_Mismatch
            DT["Similarity"] = Total_Similarity
            DT['Depth']=Total_Depth
            Last_result = pd.DataFrame.from_dict(DT, orient='index')
            #print(Last_result.T)
            """Calculate standard deviation and Pvalue"""
            data=Last_result.T
    #         print(data.sort_values(by="Total_Min" , ascending=False) )
     #####################################################################################################################
            if i==0:
                b=list(data.index)
                result={}
                result["Shared_best_reads"]=data['Total_Min']-data['Unique_Min']
                result["Unique_best_reads"]=data['Unique_Min']
                result["Similarity"]=data["Similarity"]
                result['Depth']=data['Depth']
               # print(result)
                Result=pd.DataFrame.from_dict(result,orient='index')
                D2=Result.T.round(4)
    #             print(D2)
         #       print("Count Done ！！！")
                print("Level {}, start calculating abundance\n".format(Similarity_all[i]))
                """计算 abundance"""
                matrix_flag = False
                max = 0.0
                # Get Q-matrix matrix
                data = P_matrix2

                d=list(data.columns)
                d.remove('reads')
    #             print(d)
                w=list(set(d)-set(b))
    #             print(w)
                for x in w:
                    del data[x]
                data=data.fillna(0)
                data1 = np.array(data)
                output_matrix = np.delete(data1, 0, axis=1)
                #print(output_matrix[:10])
                col_num = output_matrix.shape[1]
                m = output_matrix.shape[0]
            #    print(col_num)
                out_data = output_matrix[np.nonzero(output_matrix)]
        #         print(out_data)
                min_data = np.min(out_data)
                output_matrix[output_matrix == 0] = min_data
                # Find the index of the row maximum and the maximum value respectively
                arr=pd.DataFrame(output_matrix)
                arr['max_value']=arr.max(axis=1)
                arr['max_index']=np.argmax(output_matrix,axis=1)
                max_count = getlistnum(arr['max_index'],col_num)    # Get the number of maximum values in each column
                #print(max_count)
                x_ori = []
                for k  in max_count:        # Calculate the proportion of the maximum value of Q-value in each column as the initial value of the objective function
                    x = max_count[k]/m
                    x_ori.append(x)
                #x_ori = [0.1,0.1,0.2,0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
                new_list = [w for w in x_ori if w > 0]
                for w in range(len(x_ori)):
                    if x_ori[w] == 0:
                        x_ori[w] = min(new_list)
                # Set the objective function
                n=col_num
                # s = n
                a1=["x[" for x in range(0,n)]
                b1=[str(w) for w in range(0,n)]
                d1=["]" for w in range(0,n)]
                c=list_add(a1,b1)
                final=list_add(c,d1)
                ar=np.array(final)
                matrix_list=output_matrix.tolist()
                signal=["*" for w in range(0,n)]
                sign=[0 for w in range(0,n)]
                error_modify=len(matrix_list)
                name_list=[x for x in range(0, error_modify)]
                for j in range(0, error_modify):
                    temp = list_add(matrix_list[j], signal)
                    sign = list_add(temp, ar)
                    temp1 = "+".join(sign)
                    name_list[j] = "-math.log(" + temp1 + ")"
                func = "".join(name_list)
                f1 = "lambda x : " + func
                # Get function constraints
                func2 = "+".join(final)+"-1"
                var = [w for w in range(0, 2*n)]
                for k in range(0, n):
                    var[k] = "{'type': 'ineq', 'fun': %s x: x[%d]-0}," % ("lambda", k)
                    var[k+n] = "{'type': 'ineq', 'fun': %s x: 1- x[%d]}," % ("lambda", k)
                test = "".join(var)
                test = test[:-1]
                cons1 = "{'type': 'eq', 'fun': %s x: %s}," % ("lambda", func2)
                cons1 = "".join(cons1)
                cons2 = cons1 + test
                # Set the boundary value
                lines = []
                for w in range(0, n):
                    lines.append((0, 1))
                bound = tuple(lines)
                # Minimize the objective function
                cons = eval(cons2)
                x0 = np.array(x_ori)  # Set initial value
                res = minimize(eval(f1), x0, method='SLSQP', jac=False, bounds=bound, constraints=cons) 
               # print(res)
                result_data = []
                for j in range(0,col_num):
                    result_data.append(as_num(res.x[j]))
                boot_file=D2
                boot_file.insert(0, 'Abundance', result_data)
                boot_file["Abundance"]=boot_file["Abundance"].astype('float64')
                boot_file=boot_file.round(4)
                boot_file = boot_file[["Abundance",'Unique_best_reads','Shared_best_reads',"Similarity",'Depth']]
                boot_file=boot_file[boot_file["Abundance"]>0.05]
                boot_file.insert(0,"Cluster",boot_file.index)
                boot_file.reset_index(drop = True,inplace=True)
                boot_file['Unique_best_reads']=boot_file['Unique_best_reads'].astype(int)
                boot_file['Shared_best_reads']=boot_file['Shared_best_reads'].astype(int)
                #print(boot_file)
                boot_file.to_csv(output+'/_MIST_strain/_MIST_{}_measure.csv'.format(Similarity_all[i]))
                print("Level {}, done\n".format(Similarity_all[i]))
                
            else:
                dict_={k: list(set(g[Similarity_all[i]].tolist())) for k,g in cluster.groupby(Similarity_all[i-1])}
               # print(dict_)
                all_=set(cluster[Similarity_all[i-1]].values)
                BF=boot_file.sort_values(by="Abundance" , ascending=False)
                top_=list(map(int,BF["Cluster"].values.tolist()[:5]))
    #            print(top_)
                w=list(set(all_)-set(top_))
    #             data=data.drop(w,axis=1)
                for a in range(len(w)):
                    del dict_[w[a]] # Delete the cluster with a smaller percentage
    #             print(dict_)
                DF=pd.DataFrame();
        #         Filter_top=list()
                for key in dict_:
        #             Filter_top.append(dict_[key])
    #                 print(list(map(str,dict_[key])))
                    df=data.ix[list(map(str,dict_[key])), :]
                    if df.shape[0]>=5: 
                        df.sort_values(by="Total_Min" , ascending=False,inplace=True)
                        df1=df.iloc[:10]#;print(df1)
                        DF=DF.append(df1)
                    else:
    #                     print(df)
                        DF=DF.append(df)
    #             print(DF.sort_index())
                data=DF.sort_index()
                b=list(data.index)
                result={}
                result["Shared_best_reads"]=data['Total_Min']-data['Unique_Min']
                result["Unique_best_reads"]=data['Unique_Min']
                result["Similarity"]=data["Similarity"]
                result['Depth']=data['Depth']
               # print(result)
                Result=pd.DataFrame.from_dict(result,orient='index')
                D2=Result.T.round(4)
    #             print(D2)
         #       print("Count Done ！！！")
                print("Level {}, start calculating abundance\n".format(Similarity_all[i]))
                """计算 abundance"""
                matrix_flag = False
                max = 0.0
                #  get Q-matrix
                data = P_matrix2

        ################################################################
                d=list(data.columns)
                d.remove('reads')
    #             print(d);print(b)
                w=list(set(d)-set(b))
                data=data.drop(w,axis=1)
            ##############################################################
    #             print(data.columns)
                data=data.fillna(0)
                data1 = np.array(data)
                output_matrix = np.delete(data1, 0, axis=1)
                #print(output_matrix[:10])
                col_num = output_matrix.shape[1]
                m = output_matrix.shape[0]
            #    print(col_num)
                out_data = output_matrix[np.nonzero(output_matrix)]
        #         print(out_data)
                min_data = np.min(out_data)
                output_matrix[output_matrix == 0] = min_data
                # Find the index of the row maximum and the maximum value respectively
                arr=pd.DataFrame(output_matrix)
                arr['max_value']=arr.max(axis=1)
                arr['max_index']=np.argmax(output_matrix,axis=1)
                max_count = getlistnum(arr['max_index'],col_num)    # Get the number of maximum values in each column
    #             print(arr)
                x_ori = []
                for k  in max_count:        # Calculate the proportion of the maximum value of Q-value in each column as the initial value of the objective function
                    x = max_count[k]/m
                    x_ori.append(x)
                #x_ori = [0.1,0.1,0.2,0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
                new_list = [w for w in x_ori if w > 0]
                for w in range(len(x_ori)):
                    if x_ori[w] == 0:
                        x_ori[w] = min(new_list)
                # Set the objective function
                n=col_num
                # s = n
                a1=["x[" for x in range(0,n)]
                b1=[str(w) for w in range(0,n)]
                d1=["]" for w in range(0,n)]
                c=list_add(a1,b1)
                final=list_add(c,d1)
                ar=np.array(final)
                matrix_list=output_matrix.tolist()
                signal=["*" for w in range(0,n)]
                sign=[0 for w in range(0,n)]
                error_modify=len(matrix_list)
                name_list=[x for x in range(0, error_modify)]
                for j in range(0, error_modify):
                    temp = list_add(matrix_list[j], signal)
                    sign = list_add(temp, ar)
                    temp1 = "+".join(sign)
                    name_list[j] = "-math.log(" + temp1 + ")"
                func = "".join(name_list)
                f1 = "lambda x : " + func
                # Get function constraints
                func2 = "+".join(final)+"-1"
                var = [w for w in range(0, 2*n)]
                for k in range(0, n):
                    var[k] = "{'type': 'ineq', 'fun': %s x: x[%d]-0}," % ("lambda", k)
                    var[k+n] = "{'type': 'ineq', 'fun': %s x: 1- x[%d]}," % ("lambda", k)
                test = "".join(var)
                test = test[:-1]
                cons1 = "{'type': 'eq', 'fun': %s x: %s}," % ("lambda", func2)
                cons1 = "".join(cons1)
                cons2 = cons1 + test
                # Set the boundary value
                lines = []
                for w in range(0, n):
                    lines.append((0, 1))
                bound = tuple(lines)
                # Minimize the objective function
                cons = eval(cons2)
                x0 = np.array(x_ori)  # Set initial value
                res = minimize(eval(f1), x0, method='SLSQP', jac=False, bounds=bound, constraints=cons)  # Minimize the objective function
                result_data = []
                for j in range(0,col_num):
                    result_data.append(as_num(res.x[j]))
                boot_file=D2
                boot_file.insert(0, 'Abundance', result_data)
                boot_file["Abundance"]=boot_file["Abundance"].astype('float64')
                boot_file=boot_file.round(4)
                boot_file = boot_file[["Abundance",'Unique_best_reads','Shared_best_reads',"Similarity",'Depth']]
                boot_file=boot_file[boot_file["Abundance"]>0.05]
                boot_file.insert(0,"Cluster",boot_file.index)
                boot_file.reset_index(drop = True,inplace=True)
                boot_file['Unique_best_reads']=boot_file['Unique_best_reads'].astype(int)
                boot_file['Shared_best_reads']=boot_file['Shared_best_reads'].astype(int)
               # print(boot_file)
                boot_file.to_csv(output+'/_MIST_strain/_MIST_{}_measure.csv'.format(Similarity_all[i]))
                print("Level {}, done\n".format(Similarity_all[i]))
    else:
        for i in range(cluster.shape[1]):
     #       print(i)
            print("Level {}, start building prob matrix\n".format(Similarity_all[i]))
    ##############################################################################################
            dict_cluster=cluster.iloc[:,[i]].T.to_dict(orient='records')[0]
            for key in dict_cluster:
                dict_cluster[key] = str(dict_cluster[key])
            p_matrix=mismatch_result.copy()
            count = len(p_matrix.index)#;print(count)
            p_matrix.rename(columns=dict_cluster, inplace = True)
            List1=sorted(set(p_matrix.columns))#print(List1)
    #         List1=list(map(int,List1))
            P_matrix1=pd.DataFrame()
            for j in List1:
                if isinstance(p_matrix[j],pd.Series):
                    df=p_matrix[j].to_frame()
                if isinstance(p_matrix[j],pd.DataFrame):
                    df=p_matrix[j]
                P_matrix1[j]=df.min(axis=1)
            data1 = (0.05**P_matrix1.values)*((1-0.05)**(0-P_matrix1.values))
            P_matrix2 = pd.DataFrame(data1)
    #         print(P_matrix1.columns)
            P_matrix2.columns = P_matrix1.columns
            P_matrix2.insert(0, 'reads', P_matrix1.index)
        #        print(" {} Done\n ".format(Similarity_all[i]))
        #    print(P_matrix1.head())
            """Count the minimum number of Mismatch, the number of unique minimum values, 
               and the total Mismatch and total similarity of each cluster"""
            DT = {}
            data = P_matrix1
            Total_Min = {};
            Unique_Min = {};
            Total_Mismatch = {};
            Total_Similarity = {}
            values = data.min(axis=1)
            Keys = data.columns
            for k in Keys:
                Min = 0;
                Mismatch = 0;
                Similarity = 0
                temp = data[k]
                for j, z in zip(temp, values):
                    # print(j)
                    if j == z:  #
                        if z<read_length/20:
                            Mismatch = Mismatch+z
                            Min = Min + 1
                            if pair_1!=None and pair_2!=None and single_end==None:
                                Similarity = 1 - Mismatch / (2*read_length * Min)
                            else:
                                Similarity = 1 - Mismatch / (read_length * Min)
            
                Total_Min[k] = Min;
                Total_Mismatch[k] = Mismatch;
                Total_Similarity[k] = Similarity;
           
     #       print(Total_Similarity)
            #         print(Total_Min)
            for j in range(len(values)):
                unique_min = 0
                d = defaultdict(list)
                #         print(data1.iloc[0])
                for k, vy in [(v, y) for y, v in zip(data.columns, data.iloc[j])]:
                    d[k].append(vy)
                min_count = sorted(d.items(), key=lambda d: d[0], reverse=True)[0]
                #             print(min_count)
                if len(min_count[1]) == 1:
                    unique_min = unique_min + 1
                    if min_count[1][0] in Unique_Min.keys():
                        unique_min = Unique_Min[min_count[1][0]] + 1
                #             print(unique_min)
                Unique_Min[min_count[1][0]] = unique_min
            for k in Keys:
                if k not in Unique_Min.keys():
                    Unique_Min[k] = 0
            DT["Total_Min"] = Total_Min
            DT["Unique_Min"] = Unique_Min
            DT["Total_Mismatch"] = Total_Mismatch
            DT["Similarity"] = Total_Similarity
            Last_result = pd.DataFrame.from_dict(DT, orient='index')
           # print(Last_result.T)
            """Calculate standard deviation and Pvalue"""
            data=Last_result.T
    #         print(data.sort_values(by="Total_Min" , ascending=False) )
     #####################################################################################################################
            if i==0:
                b=list(data.index)
                result={}
                result["Shared_best_reads"]=data['Total_Min']-data['Unique_Min']
                result["Unique_best_reads"]=data['Unique_Min']
                result["Similarity"]=data["Similarity"]
               # print(result)
                Result=pd.DataFrame.from_dict(result,orient='index')
                D2=Result.T.round(4)
                print("Level {}, start calculating abundance\n".format(Similarity_all[i]))
                """计算 abundance"""
                matrix_flag = False
                max = 0.0
                # Get Q-matrix matrix
                data = P_matrix2

                d=list(data.columns)
                d.remove('reads')
    #             print(d)
                w=list(set(d)-set(b))
    #             print(w)
                for x in w:
                    del data[x]
                data=data.fillna(0)
                data1 = np.array(data)
                output_matrix = np.delete(data1, 0, axis=1)
                #print(output_matrix[:10])
                col_num = output_matrix.shape[1]
                m = output_matrix.shape[0]
            #    print(col_num)
                out_data = output_matrix[np.nonzero(output_matrix)]
        #         print(out_data)
                min_data = np.min(out_data)
                output_matrix[output_matrix == 0] = min_data
                # Find the index of the row maximum and the maximum value respectively
                arr=pd.DataFrame(output_matrix)
                arr['max_value']=arr.max(axis=1)
                arr['max_index']=np.argmax(output_matrix,axis=1)
                max_count = getlistnum(arr['max_index'],col_num)    # Get the number of maximum values in each column
                #print(max_count)
                x_ori = []
                for k  in max_count:        # Calculate the proportion of the maximum value of Q-value in each column as the initial value of the objective function
                    x = max_count[k]/m
                    x_ori.append(x)
                #x_ori = [0.1,0.1,0.2,0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
                new_list = [w for w in x_ori if w > 0]
                for w in range(len(x_ori)):
                    if x_ori[w] == 0:
                        x_ori[w] = min(new_list)
                # Set the objective function
                n=col_num
                # s = n
                a1=["x[" for x in range(0,n)]
                b1=[str(w) for w in range(0,n)]
                d1=["]" for w in range(0,n)]
                c=list_add(a1,b1)
                final=list_add(c,d1)
                ar=np.array(final)
                matrix_list=output_matrix.tolist()
                signal=["*" for w in range(0,n)]
                sign=[0 for w in range(0,n)]
                error_modify=len(matrix_list)
                name_list=[x for x in range(0, error_modify)]
                for j in range(0, error_modify):
                    temp = list_add(matrix_list[j], signal)
                    sign = list_add(temp, ar)
                    temp1 = "+".join(sign)
                    name_list[j] = "-math.log(" + temp1 + ")"
                func = "".join(name_list)
                f1 = "lambda x : " + func
                # Get function constraints
                func2 = "+".join(final)+"-1"
                var = [w for w in range(0, 2*n)]
                for k in range(0, n):
                    var[k] = "{'type': 'ineq', 'fun': %s x: x[%d]-0}," % ("lambda", k)
                    var[k+n] = "{'type': 'ineq', 'fun': %s x: 1- x[%d]}," % ("lambda", k)
                test = "".join(var)
                test = test[:-1]
                cons1 = "{'type': 'eq', 'fun': %s x: %s}," % ("lambda", func2)
                cons1 = "".join(cons1)
                cons2 = cons1 + test
                # Set the boundary value
                lines = []
                for w in range(0, n):
                    lines.append((0, 1))
                bound = tuple(lines)
                # Minimize the objective function
                cons = eval(cons2)
                x0 = np.array(x_ori)  # Set initial value
                res = minimize(eval(f1), x0, method='SLSQP', jac=False, bounds=bound, constraints=cons) 
               # print(res)
                result_data = []
                for j in range(0,col_num):
                    result_data.append(as_num(res.x[j]))
                boot_file=D2
                boot_file.insert(0, 'Abundance', result_data)
                boot_file["Abundance"]=boot_file["Abundance"].astype('float64')
                boot_file=boot_file.round(4)
                boot_file = boot_file[["Abundance",'Unique_best_reads','Shared_best_reads',"Similarity"]]
                boot_file=boot_file[boot_file["Abundance"]>0.05]
                boot_file.insert(0,"Cluster",boot_file.index)
                boot_file.reset_index(drop = True,inplace=True)
                boot_file['Unique_best_reads']=boot_file['Unique_best_reads'].astype(int)
                boot_file['Shared_best_reads']=boot_file['Shared_best_reads'].astype(int)
            #    print(boot_file)
                boot_file.to_csv(output+'/_MIST_strain/_MIST_{}_measure.csv'.format(Similarity_all[i]))
                print("Level {}, done\n".format(Similarity_all[i]))
                
            else:
                dict_={k: list(set(g[Similarity_all[i]].tolist())) for k,g in cluster.groupby(Similarity_all[i-1])}
                all_=set(cluster[Similarity_all[i-1]].values)
                #print(all)
                BF=boot_file.sort_values(by="Abundance" , ascending=False)
                top_=list(map(int,BF["Cluster"].values.tolist()[:5]))
    #            print(top_)
                w=list(set(all_)-set(top_))
    #             data=data.drop(w,axis=1)
                for a in range(len(w)):
                    del dict_[w[a]] # Delete the cluster with a smaller percentage
    #             print(dict_)
                DF=pd.DataFrame();
        #         Filter_top=list()
                for key in dict_:
        #             Filter_top.append(dict_[key])
    #                 print(list(map(str,dict_[key])))
                    df=data.ix[list(map(str,dict_[key])), :]
                    if df.shape[0]>=5: 
                        df.sort_values(by="Total_Min" , ascending=False,inplace=True)
                        df1=df.iloc[:10]#;print(df1)
                        DF=DF.append(df1)
                    else:
    #                     print(df)
                        DF=DF.append(df)
   #             print(DF.sort_index())
                data=DF.sort_index()
                b=list(data.index)
                result={}
                result["Shared_best_reads"]=data['Total_Min']-data['Unique_Min']
                result["Unique_best_reads"]=data['Unique_Min']
                result["Similarity"]=data["Similarity"]
               # print(result)
                Result=pd.DataFrame.from_dict(result,orient='index')
                D2=Result.T.round(4)
                #print(D2)
         #       print("Count Done ！！！")
                print("Level {}, start calculating abundance\n".format(Similarity_all[i]))
                """计算 abundance"""
                matrix_flag = False
                max = 0.0
                #  get Q-matrix
                data = P_matrix2

        ################################################################
                d=list(data.columns)
                d.remove('reads')
    #             print(d);print(b)
                w=list(set(d)-set(b))
                data=data.drop(w,axis=1)
            ##############################################################
    #             print(data.columns)
                data=data.fillna(0)
                data1 = np.array(data)
                output_matrix = np.delete(data1, 0, axis=1)
                #print(output_matrix[:10])
                col_num = output_matrix.shape[1]
                m = output_matrix.shape[0]
            #    print(col_num)
                out_data = output_matrix[np.nonzero(output_matrix)]
        #         print(out_data)
                min_data = np.min(out_data)
                output_matrix[output_matrix == 0] = min_data
                # Find the index of the row maximum and the maximum value respectively
                arr=pd.DataFrame(output_matrix)
                arr['max_value']=arr.max(axis=1)
                arr['max_index']=np.argmax(output_matrix,axis=1)
                max_count = getlistnum(arr['max_index'],col_num)    # Get the number of maximum values in each column
    #             print(arr)
                x_ori = []
                for k  in max_count:        # Calculate the proportion of the maximum value of Q-value in each column as the initial value of the objective function
                    x = max_count[k]/m
                    x_ori.append(x)
                #x_ori = [0.1,0.1,0.2,0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
                new_list = [w for w in x_ori if w > 0]
                for w in range(len(x_ori)):
                    if x_ori[w] == 0:
                        x_ori[w] = min(new_list)
                # Set the objective function
                n=col_num
                # s = n
                a1=["x[" for x in range(0,n)]
                b1=[str(w) for w in range(0,n)]
                d1=["]" for w in range(0,n)]
                c=list_add(a1,b1)
                final=list_add(c,d1)
                ar=np.array(final)
                matrix_list=output_matrix.tolist()
                signal=["*" for w in range(0,n)]
                sign=[0 for w in range(0,n)]
                error_modify=len(matrix_list)
                name_list=[x for x in range(0, error_modify)]
                for j in range(0, error_modify):
                    temp = list_add(matrix_list[j], signal)
                    sign = list_add(temp, ar)
                    temp1 = "+".join(sign)
                    name_list[j] = "-math.log(" + temp1 + ")"
                func = "".join(name_list)
                f1 = "lambda x : " + func
                # Get function constraints
                func2 = "+".join(final)+"-1"
                var = [w for w in range(0, 2*n)]
                for k in range(0, n):
                    var[k] = "{'type': 'ineq', 'fun': %s x: x[%d]-0}," % ("lambda", k)
                    var[k+n] = "{'type': 'ineq', 'fun': %s x: 1- x[%d]}," % ("lambda", k)
                test = "".join(var)
                test = test[:-1]
                cons1 = "{'type': 'eq', 'fun': %s x: %s}," % ("lambda", func2)
                cons1 = "".join(cons1)
                cons2 = cons1 + test
                # Set the boundary value
                lines = []
                for w in range(0, n):
                    lines.append((0, 1))
                bound = tuple(lines)
                # Minimize the objective function
                cons = eval(cons2)
                x0 = np.array(x_ori)  # Set initial value
                res = minimize(eval(f1), x0, method='SLSQP', jac=False, bounds=bound, constraints=cons)  # Minimize the objective function
                result_data = []
                for j in range(0,col_num):
                    result_data.append(as_num(res.x[j]))
                #D2=D2.dropna(axis=0,how='all')  
                boot_file=D2
                #print(boot_file)
                #print(result_data)
                boot_file.insert(0, 'Abundance', result_data)
                boot_file["Abundance"]=boot_file["Abundance"].astype('float64')
                boot_file=boot_file.round(4)
                boot_file = boot_file[["Abundance",'Unique_best_reads','Shared_best_reads',"Similarity"]]
                boot_file=boot_file[boot_file["Abundance"]>0.05]
                boot_file.insert(0,"Cluster",boot_file.index)
                boot_file.reset_index(drop = True,inplace=True)
                boot_file['Unique_best_reads']=boot_file['Unique_best_reads'].astype(int)
                boot_file['Shared_best_reads']=boot_file['Shared_best_reads'].astype(int)
               # print(boot_file)
                boot_file.to_csv(output+'/_MIST_strain/_MIST_{}_measure.csv'.format(Similarity_all[i]))
                print("Level {}, done\n".format(Similarity_all[i]))
    print("All done \n")

