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

def subspecies(cluster_output,mismatch_matrix_output,read_length,output):
    if not os.path.exists(output):
        os.makedirs(output)
    """Separate each cluster  """
    cluster=pd.read_csv(cluster_output,index_col=0)# read data
    Similarity_all=cluster.columns.values.tolist()
#    f_time=open(output+"/_MIST_measure_time.txt","a")
    for i in range(cluster.shape[1]):
        print("Level {}, start building prob matrix\n".format(Similarity_all[i]))
#         print(Similarity_all[i])
#        start = datetime.datetime.now()
##############################################################################################
        dict_cluster=cluster.iloc[:,[i]].T.to_dict(orient='records')[0]
        for key in dict_cluster:
            dict_cluster[key] = str(dict_cluster[key])
        p_matrix=pd.read_csv(mismatch_matrix_output,index_col=0)
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
                    Mismatch = Mismatch+z
                    Min = Min + 1
                   # print(Mismatch);print(Min)
                    Similarity = 1 - Mismatch / (read_length * Min)
            Total_Min[k] = Min;
            Total_Mismatch[k] = Mismatch;
            Total_Similarity[k] = Similarity
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
        #print(Last_result.T)
        """Calculate standard deviation and Pvalue"""
        data=Last_result.T
#         print(data.sort_values(by="Total_Min" , ascending=False) )
 #####################################################################################################################
        if i==0:
            b=list(data.index)
            result={}
            result1={};result2={}
            _ROOT = os.path.abspath(os.path.dirname(__file__))
            for a in range(len(data.index)):
                # Load model files and generate model objects
                new_model1 = joblib.load(_ROOT+"/STD_model.pkl")
                new_pred_data1 = [[Similarity_all[i], data.iloc[a,3], data.iloc[a,1]]]
                # Use the generated model to predict new samples
                test_result1 = new_model1.predict(new_pred_data1)
                if test_result1[0]>1:
                    test_result1[0]=1
                    #print(test_result1)
                    result1[data.index[a]]=test_result1[0]
                elif test_result1[0]<0:
                    #print(test_result1)
                    test_result1[0]=0
                    result1[data.index[a]]=test_result1[0]
                else:
                    result1[data.index[a]]=test_result1[0]
                # Load model files and generate model objects
                new_model2 = joblib.load(_ROOT+"/P_model.pkl")

                new_pred_data2 = [[Similarity_all[i], data.iloc[a,3], data.iloc[a,1]]]
                # Use the generated model to predict new samples
                test_result2 = new_model2.predict(new_pred_data2)
                if test_result2[0]<0:
                    test_result2[0]=0
                    result2[data.index[a]]=test_result2[0]
                else:
                    result2[data.index[a]]=test_result2[0]
#             print("P value :")
#             print(result2)
            Pct_STDEV_Pct_Average=result1
            result["P-value"]=result2
            result["Similarity"]=data["Similarity"]
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
            result_data = []
            for j in range(0,col_num):
                result_data.append(as_num(res.x[j]))
#            re1 = map(result_data.index, heapq.nlargest(5, result_data))  # Find the five largest indexes
                                                                          # nsmallest is the opposite of nlargest, seeking the smallest
#             top_=list(re1)
#            re2 = heapq.nlargest(5, result_data)  # 求最大的五个元素
#             print(list(re1));print(re2)
            #print(result)
            #Result=pd.DataFrame.from_dict(result,orient='index')
           # Result=Result.T
            LC={};UC={}
            for x in range(len(Pct_STDEV_Pct_Average)):
                X=float(result_data[x])
                Y=list(Pct_STDEV_Pct_Average.values())[x]
     #           print(X,Y)
                Z=X-1.96*((Y*X)/10);W=X+1.96*((Y*X)/10)
                key=list(Pct_STDEV_Pct_Average.keys())[x]
                LC[key]=Z;UC[key]=W
        #    print(LC,UC)
#            print(result_data)
            boot_file=D2
            boot_file.insert(0,'LowerConfidence', list(LC.values()))
            boot_file.insert(0,'UpperConfidence', list(UC.values()))
            boot_file.insert(0, 'Abundance', result_data)
            boot_file["Abundance"]=boot_file["Abundance"].astype('float64')
            boot_file=boot_file.round(4)
            boot_file = boot_file[["Abundance","LowerConfidence","UpperConfidence","P-value","Similarity"]]
            boot_file=boot_file[boot_file["Abundance"]>0.05]
            boot_file.insert(0,"Cluster",boot_file.index)
            boot_file.reset_index(drop = True,inplace=True)
#            print(boot_file)
            boot_file.to_csv(output+'/_MIST_{}_measure.csv'.format(Similarity_all[i]))
            print("Level {}, done\n".format(Similarity_all[i]))
#            end = datetime.datetime.now()
#            spend_time=end-start
#            print("Level {} spend time: {}".format(Similarity_all[i],str(spend_time)))
#            f_time.write("Level {} spend time: {}".format(Similarity_all[i],str(spend_time)))
#            f_time.write("\n")
            
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
            result1={};result2={}
            _ROOT = os.path.abspath(os.path.dirname(__file__))
            for a in range(len(data.index)):
                # Load model files and generate model objects
                new_model1 = joblib.load(_ROOT+"/STD_model.pkl")
                new_pred_data1 = [[Similarity_all[i], data.iloc[a,3], data.iloc[a,1]]]
                # Use the generated model to predict new samples
                test_result1 = new_model1.predict(new_pred_data1)
                if test_result1[0]>1:
                    test_result1[0]=1
                    #print(test_result1)
                    result1[data.index[a]]=test_result1[0]
                elif test_result1[0]<0:
                    #print(test_result1)
                    test_result1[0]=0
                    result1[data.index[a]]=test_result1[0]
                else:
                    result1[data.index[a]]=test_result1[0]
                # Load model files and generate model objects
                new_model2 = joblib.load(_ROOT+"/P_model.pkl")
                new_pred_data2 = [[Similarity_all[i], data.iloc[a,3], data.iloc[a,1]]]
                # Use the generated model to predict new samples
                test_result2 = new_model2.predict(new_pred_data2)
                if test_result2[0]<0:
                    test_result2[0]=0
                    result2[data.index[a]]=test_result2[0]
                else:
                    result2[data.index[a]]=test_result2[0]
    #        print(result1)
            Pct_STDEV_Pct_Average=result1
            result["P-value"]=result2
            result["Similarity"]=data["Similarity"]
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
            
#             w=list(map(int,w));print(w)
#             print(data.columns)
            data=data.drop(w,axis=1)
#             for x in w:
#                 del data[x]
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
#             re1 = map(result_data.index, heapq.nlargest(10, result_data))  
#            re2 = heapq.nlargest(5, result_data)  # 
#            top_=list(re1);print(top_)
            #print(result)
            #Result=pd.DataFrame.from_dict(result,orient='index')
           # Result=Result.T
            LC={};UC={}
            for x in range(len(Pct_STDEV_Pct_Average)):
                X=float(result_data[x])
                Y=list(Pct_STDEV_Pct_Average.values())[x]
     #           print(X,Y)
                Z=X-1.96*((Y*X)/10);W=X+1.96*((Y*X)/10)
                key=list(Pct_STDEV_Pct_Average.keys())[x]
                LC[key]=Z;UC[key]=W
        #    print(LC,UC)
#             print(result_data)
            boot_file=D2
            boot_file.insert(0,'LowerConfidence', list(LC.values()))
            boot_file.insert(0,'UpperConfidence', list(UC.values()))
            boot_file.insert(0, 'Abundance', result_data)
            boot_file["Abundance"]=boot_file["Abundance"].astype('float64')
            boot_file=boot_file.round(4)
            boot_file = boot_file[["Abundance","LowerConfidence","UpperConfidence","P-value","Similarity"]]
            boot_file=boot_file[boot_file["Abundance"]>0.05]
            boot_file.insert(0,"Cluster",boot_file.index)
            boot_file.reset_index(drop = True,inplace=True)
#            print(boot_file)
            boot_file.to_csv(output+'/_MIST_{}_measure.csv'.format(Similarity_all[i]))
            print("Level {}, done\n".format(Similarity_all[i]))
#            end = datetime.datetime.now()
#            spend_time=end-start
#            print("Level {} spend time: {}".format(Similarity_all[i],str(spend_time)))
#            f_time.write("Level {} spend time: {}".format(Similarity_all[i],str(spend_time)))
#            f_time.write("\n")
            
    print("All done \n")
#    f_time.close()

