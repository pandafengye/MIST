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
import warnings
sys.setrecursionlimit(30000) #
warnings.filterwarnings('ignore')
def as_num(x):
    y = '{:.10f}'.format(x)  # .10f 保留10位小数
    return y
def getlistnum(li, num):  # 这个函数就是要对列表的每个元素进行计数
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
def measure(cluster_output,mismatch_matrix_output,read_length,output):
    """将每个聚类分开"""
    cluster=pd.read_csv(cluster_output,index_col=0)#读取数据
    Similarity_all=cluster.columns.values.tolist()
    
    for i in range(cluster.shape[1]):
        print("Level {}, start building prob matrix\n".format(Similarity_all[i]))
#         print(Similarity_all[i])
        dict_cluster=cluster.iloc[:,[i]].T.to_dict(orient='records')[0]
        for key in dict_cluster:
            dict_cluster[key] = str(dict_cluster[key])
        p_matrix=pd.read_csv(mismatch_matrix_output,index_col=0)
        p_matrix.rename(columns=dict_cluster, inplace = True)
        List1=sorted(set(p_matrix.columns))
        P_matrix1=pd.DataFrame()
        for j in List1:
            if isinstance(p_matrix[j],pd.Series):
                df=p_matrix[j].to_frame()
            if isinstance(p_matrix[j],pd.DataFrame):
                df=p_matrix[j]
            P_matrix1[j]=df.min(axis=1)
        data1 = (0.01**P_matrix1.values)*(0.99**(0-P_matrix1.values))
        P_matrix2 = pd.DataFrame(data1)
        P_matrix2.insert(0, 'reads', P_matrix1.index)
#        print(" {} Done\n ".format(Similarity_all[i]))
#         print(P_matrix2.head())
        """统计每个聚类的Mismatch最小值个数、唯一最小值个数以及总Mismatch和总的相似性"""
 #       print("Count the number of Mismatch minimums, the number of unique minimums,the total Mismatch, and the total similarity for {} cluster.\n".format(Similarity_all[i]))
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
                    Min = Min + 1;
                    Similarity = 1 - Mismatch / (read_length * Min)
            Total_Min[k] = Min;
            Total_Mismatch[k] = Mismatch;
            Total_Similarity[k] = Similarity

        #         print(Total_Min)
        for j in range(len(values)):
            unique_min = 0
            d = defaultdict(list)
            #         print(data1.iloc[0])
            for k, va in [(v, a) for a, v in zip(data.columns, data.iloc[j])]:
                d[k].append(va)
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
#         print(Last_result.T)
        """计算标准差以及Pvalue """
        data=Last_result.T
        result={}
        result1={};result2={}
        _ROOT = os.path.abspath(os.path.dirname(__file__))
        for a in range(len(data.index)):
            # 加载模型文件，生成模型对象
            new_model1 = joblib.load(_ROOT+"/STD_model.pkl")
            new_pred_data1 = [[Similarity_all[i], data.iloc[a,1], data.iloc[a,0]]]
            # 使用加载生成的模型预测新样本
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
            # 加载模型文件，生成模型对象
            new_model2 = joblib.load(_ROOT+"/P_model.pkl")
            
            new_pred_data2 = [[Similarity_all[i], data.iloc[a,1], data.iloc[a,0]]]
            # 使用加载生成的模型预测新样本
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
 #       print("Count Done ！！！")
        print("Level {}, start calculating abundance\n".format(Similarity_all[i]))
        """计算 abundance"""
        matrix_flag = False
        max = 0.0
        #获得Q-matrix矩阵
        data = P_matrix2
        data=data.fillna(0)
        data1 = np.array(data)
        output_matrix = np.delete(data1, 0, axis=1)
        col_num = output_matrix.shape[1]
        m = output_matrix.shape[0]
    #    print(col_num)
        out_data = output_matrix[np.nonzero(output_matrix)]
        min_data = np.min(out_data)
        output_matrix[output_matrix == 0] = min_data
        #分别求行最大值及最大值所在索引
        arr=pd.DataFrame(output_matrix)
        arr['max_value']=arr.max(axis=1)
        arr['max_index']=np.argmax(output_matrix,axis=1)
        max_count = getlistnum(arr['max_index'],col_num)    #获取每列存在最大值的个数
        x_ori = []
        for k  in max_count:        #计算每列中存在每行Q-value最大值的占比作为目标函数初始值
            x = max_count[k]/m
            x_ori.append(x)
        #x_ori = [0.1,0.1,0.2,0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
        new_list = [w for w in x_ori if w > 0]
        for w in range(len(x_ori)):
            if x_ori[w] == 0:
                x_ori[w] = min(new_list)
        #设置目标函数
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
        # 获取函数约束条件
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
        # 设置边界值
        lines = []
        for w in range(0, n):
            lines.append((0, 1))
        bound = tuple(lines)
        #最小化目标函数
        cons = eval(cons2)
        x0 = np.array(x_ori)  # 设置初始值
        res = minimize(eval(f1), x0, method='SLSQP', jac=False, bounds=bound, constraints=cons)  # 最小化目标函数
        result_data = []
        for j in range(0,col_num):
            result_data.append(as_num(res.x[j]))
        re1 = map(result_data.index, heapq.nlargest(5, result_data))  # 求最大的五个索引    nsmallest与nlargest相反，求最小
        re2 = heapq.nlargest(5, result_data)  # 求最大的五个元素
        
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
        boot_file=D2
        boot_file.insert(0,'LowerConfidence', list(LC.values()))
        boot_file.insert(0,'UpperConfidence', list(UC.values()))
        boot_file.insert(0, 'Abundance', result_data)
        boot_file["Abundance"]=boot_file["Abundance"].astype('float64')
        boot_file=boot_file.round(4)
        boot_file = boot_file[["Abundance","LowerConfidence","UpperConfidence","P-value","Similarity"]]
        print(boot_file)
        boot_file.to_csv(output+'/_MIST_{}_measure.csv'.format(Similarity_all[i]))
        print("Level {}, done\n".format(Similarity_all[i]))
    print("All done \n")
