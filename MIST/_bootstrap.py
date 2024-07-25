import pandas as pd
import os,fnmatch

import warnings
from scipy import stats
from collections import defaultdict
warnings.simplefilter("ignore", category=DeprecationWarning)
import os, math, re, sys, datetime, heapq
from collections import defaultdict
import pandas as pd
import numpy as np
from scipy.optimize import minimize
from numpy import nan
import matplotlib.pyplot as plt
from sklearn import datasets, ensemble, naive_bayes
from sklearn.model_selection import train_test_split
from sklearn.externals import joblib

sys.setrecursionlimit(100000)
import os
import os.path
import click
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
import sys

sys.setrecursionlimit(100000)
import faulthandler

faulthandler.enable()


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

def strain(cluster_output, mismatch_matrix, read_length, output):
    if not os.path.exists(output + "/_MIST_strain/"):
        os.makedirs(output + "/_MIST_strain/")
    """Separate each cluster  """
    mismatch_result = pd.read_csv(mismatch_matrix, index_col=0)  # read data
    cluster = pd.read_csv(cluster_output, index_col=0)  # read data
    # ### panda:   这句是读入MIST_ref_cluster.csv文件，这个是需要保留的，是不是意味着，从这里开始，后面都是要保留的？
    Similarity_all = cluster.columns.values.tolist()

    for i in range(cluster.shape[1]):
        #       print(i)
        print("Level {}, start building prob matrix\n".format(Similarity_all[i]))
        ##############################################################################################
        dict_cluster = cluster.iloc[:, [i]].T.to_dict(orient='records')[0]
        for key in dict_cluster:
            dict_cluster[key] = str(dict_cluster[key])
        p_matrix = mismatch_result.copy()  ### panda:  这里再次出现了mismatch_result.copy，需要改一下
        count = len(p_matrix.index)  # ;print(count)
        p_matrix.rename(columns=dict_cluster, inplace=True)
        List1 = sorted(set(p_matrix.columns))  # print(List1)
        #         List1=list(map(int,List1))
        P_matrix1 = pd.DataFrame()
        for j in List1:
            if isinstance(p_matrix[j], pd.Series):
                df = p_matrix[j].to_frame()
            if isinstance(p_matrix[j], pd.DataFrame):
                df = p_matrix[j]
            P_matrix1[j] = df.min(axis=1)
        data1 = (0.05 ** P_matrix1.values) * ((1 - 0.05) ** (0 - P_matrix1.values))
        P_matrix2 = pd.DataFrame(data1)
        #         print(P_matrix1.columns)
        P_matrix2.columns = P_matrix1.columns
        P_matrix2.insert(0, 'reads', P_matrix1.index)
        #        print(" {} Done\n ".format(Similarity_all[i]))
        #    print(P_matrix1.head())
        """Count the minimum number of Mismatch, the number of unique minimum values, 
           and the total Mismatch and total similarity of each cluster"""
        DT = {}
        data = P_matrix1.copy()
        Total_Min = {};
        Unique_Min = {};
        Total_Mismatch = {};
        Total_Similarity = {}
        values = data.min(axis=1)
        Keys = data.columns
        # Diff_cutoff = read_length;  # cutoff for filtering false match reads
        Diff_cutoff = read_length / 10

        for k in Keys:
            Min = 0;
            Mismatch = 0;
            Similarity = 0
            temp = data[k]
            for j, z in zip(temp, values):
                # print(j)
                if j == z:  #
                    if z < Diff_cutoff:
                        Mismatch = Mismatch + z
                        Min = Min + 1
                        # if pair_1!=None and pair_2!=None and single_end==None:
                        #     Similarity = 1 - Mismatch / (2*read_length * Min)
                        # else:
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
            # print(min_count);print(Diff_cutoff)
            if min_count[0] < Diff_cutoff:
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
        data = Last_result.T
        #         print(data.sort_values(by="Total_Min" , ascending=False) )
        #####################################################################################################################
        if i == 0:
            b = list(data.index)
            result = {}
            result["Shared_best_reads"] = data['Total_Min'] - data['Unique_Min']
            result["Unique_best_reads"] = data['Unique_Min']
            result["Similarity"] = data["Similarity"]
            # print(result)
            Result = pd.DataFrame.from_dict(result, orient='index')
            D2 = Result.T.round(4)
            print("Level {}, start calculating abundance\n".format(Similarity_all[i]))
            """计算 abundance"""
            matrix_flag = False
            max = 0.0
            # Get Q-matrix matrix
            data = P_matrix2

            d = list(data.columns)
            d.remove('reads')
            #             print(d)
            w = list(set(d) - set(b))
            #             print(w)
            for x in w:
                del data[x]
            data = data.fillna(0)
            data1 = np.array(data)
            output_matrix = np.delete(data1, 0, axis=1)
            # print(output_matrix[:10])
            col_num = output_matrix.shape[1]
            m = output_matrix.shape[0]
            #    print(col_num)
            out_data = output_matrix[np.nonzero(output_matrix)]
            #         print(out_data)
            min_data = np.min(out_data)
            output_matrix[output_matrix == 0] = min_data
            # Find the index of the row maximum and the maximum value respectively
            arr = pd.DataFrame(output_matrix)
            arr['max_value'] = arr.max(axis=1)
            arr['max_index'] = np.argmax(output_matrix, axis=1)
            max_count = getlistnum(arr['max_index'], col_num)  # Get the number of maximum values in each column
            # print(max_count)
            x_ori = []
            for k in max_count:  # Calculate the proportion of the maximum value of Q-value in each column as the initial value of the objective function
                x = max_count[k] / m
                x_ori.append(x)
            # x_ori = [0.1,0.1,0.2,0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
            new_list = [w for w in x_ori if w > 0]
            for w in range(len(x_ori)):
                if x_ori[w] == 0:
                    x_ori[w] = min(new_list)
            # Set the objective function
            n = col_num
            # s = n
            a1 = ["x[" for x in range(0, n)]
            b1 = [str(w) for w in range(0, n)]
            d1 = ["]" for w in range(0, n)]
            c = list_add(a1, b1)
            final = list_add(c, d1)
            ar = np.array(final)
            matrix_list = output_matrix.tolist()
            signal = ["*" for w in range(0, n)]
            sign = [0 for w in range(0, n)]
            error_modify = len(matrix_list)
            name_list = [x for x in range(0, error_modify)]
            for j in range(0, error_modify):
                temp = list_add(matrix_list[j], signal)
                sign = list_add(temp, ar)
                temp1 = "+".join(sign)
                name_list[j] = "-math.log(" + temp1 + ")"
            func = "".join(name_list)
            f1 = "lambda x : " + func
            # Get function constraints
            func2 = "+".join(final) + "-1"
            var = [w for w in range(0, 2 * n)]
            for k in range(0, n):
                var[k] = "{'type': 'ineq', 'fun': %s x: x[%d]-0}," % ("lambda", k)
                var[k + n] = "{'type': 'ineq', 'fun': %s x: 1- x[%d]}," % ("lambda", k)
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
            for j in range(0, col_num):
                result_data.append(as_num(res.x[j]))
            boot_file = D2
            boot_file.insert(0, 'Abundance', result_data)
            boot_file["Abundance"] = boot_file["Abundance"].astype('float64')
            boot_file = boot_file.round(4)
            boot_file = boot_file[["Abundance", 'Unique_best_reads', 'Shared_best_reads', "Similarity"]]
            boot_file = boot_file[boot_file["Abundance"] > 0.05]
            boot_file.insert(0, "Cluster", boot_file.index)
            boot_file.reset_index(drop=True, inplace=True)
            boot_file['Unique_best_reads'] = boot_file['Unique_best_reads'].astype(int)
            boot_file['Shared_best_reads'] = boot_file['Shared_best_reads'].astype(int)
            #    print(boot_file)
            boot_file.to_csv(output + '/_MIST_strain/_MIST_{}_measure.csv'.format(Similarity_all[i]))
            print("Level {}, done\n".format(Similarity_all[i]))

        else:
            dict_ = {k: list(set(g[Similarity_all[i]].tolist())) for k, g in cluster.groupby(Similarity_all[i - 1])}
            all_ = set(cluster[Similarity_all[i - 1]].values)
            # print(all)
            BF = boot_file.sort_values(by="Abundance", ascending=False)
            top_ = list(map(int, BF["Cluster"].values.tolist()[:5]))
            #            print(top_)
            w = list(set(all_) - set(top_))
            #             data=data.drop(w,axis=1)
            for a in range(len(w)):
                del dict_[w[a]]  # Delete the cluster with a smaller percentage
            #             print(dict_)
            DF = pd.DataFrame();
            #         Filter_top=list()
            for key in dict_:
                #             Filter_top.append(dict_[key])
                #                 print(list(map(str,dict_[key])))
                df = data.ix[list(map(str, dict_[key])), :]
                if df.shape[0] >= 5:
                    df.sort_values(by="Total_Min", ascending=False, inplace=True)
                    df1 = df.iloc[:10]  # ;print(df1)
                    DF = DF.append(df1)
                else:
                    #                     print(df)
                    DF = DF.append(df)
            #             print(DF.sort_index())
            data = DF.sort_index()
            b = list(data.index)
            result = {}
            result["Shared_best_reads"] = data['Total_Min'] - data['Unique_Min']
            result["Unique_best_reads"] = data['Unique_Min']
            result["Similarity"] = data["Similarity"]
            # print(result)
            Result = pd.DataFrame.from_dict(result, orient='index')
            D2 = Result.T.round(4)
            # print(D2)
            #       print("Count Done ！！！")
            print("Level {}, start calculating abundance\n".format(Similarity_all[i]))
            """计算 abundance"""
            matrix_flag = False
            max = 0.0
            #  get Q-matrix
            data = P_matrix2

            ################################################################
            d = list(data.columns)
            d.remove('reads')
            #             print(d);print(b)
            w = list(set(d) - set(b))
            data = data.drop(w, axis=1)
            ##############################################################
            #             print(data.columns)
            data = data.fillna(0)
            data1 = np.array(data)
            output_matrix = np.delete(data1, 0, axis=1)
            # print(output_matrix[:10])
            col_num = output_matrix.shape[1]
            m = output_matrix.shape[0]
            #    print(col_num)
            out_data = output_matrix[np.nonzero(output_matrix)]
            #         print(out_data)
            min_data = np.min(out_data)
            output_matrix[output_matrix == 0] = min_data
            # Find the index of the row maximum and the maximum value respectively
            arr = pd.DataFrame(output_matrix)
            arr['max_value'] = arr.max(axis=1)
            arr['max_index'] = np.argmax(output_matrix, axis=1)
            max_count = getlistnum(arr['max_index'], col_num)  # Get the number of maximum values in each column
            #             print(arr)
            x_ori = []
            for k in max_count:  # Calculate the proportion of the maximum value of Q-value in each column as the initial value of the objective function
                x = max_count[k] / m
                x_ori.append(x)
            # x_ori = [0.1,0.1,0.2,0.8,0.1,0.1,0.1,0.1,0.1,0.1,0.1]
            new_list = [w for w in x_ori if w > 0]
            for w in range(len(x_ori)):
                if x_ori[w] == 0:
                    x_ori[w] = min(new_list)
            # Set the objective function
            n = col_num
            # s = n
            a1 = ["x[" for x in range(0, n)]
            b1 = [str(w) for w in range(0, n)]
            d1 = ["]" for w in range(0, n)]
            c = list_add(a1, b1)
            final = list_add(c, d1)
            ar = np.array(final)
            matrix_list = output_matrix.tolist()
            signal = ["*" for w in range(0, n)]
            sign = [0 for w in range(0, n)]
            error_modify = len(matrix_list)
            name_list = [x for x in range(0, error_modify)]
            for j in range(0, error_modify):
                temp = list_add(matrix_list[j], signal)
                sign = list_add(temp, ar)
                temp1 = "+".join(sign)
                name_list[j] = "-math.log(" + temp1 + ")"
            func = "".join(name_list)
            f1 = "lambda x : " + func
            # Get function constraints
            func2 = "+".join(final) + "-1"
            var = [w for w in range(0, 2 * n)]
            for k in range(0, n):
                var[k] = "{'type': 'ineq', 'fun': %s x: x[%d]-0}," % ("lambda", k)
                var[k + n] = "{'type': 'ineq', 'fun': %s x: 1- x[%d]}," % ("lambda", k)
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
            res = minimize(eval(f1), x0, method='SLSQP', jac=False, bounds=bound,
                           constraints=cons)  # Minimize the objective function
            result_data = []
            for j in range(0, col_num):
                result_data.append(as_num(res.x[j]))
            # D2=D2.dropna(axis=0,how='all')
            boot_file = D2
            # print(boot_file)
            # print(result_data)
            boot_file.insert(0, 'Abundance', result_data)
            boot_file["Abundance"] = boot_file["Abundance"].astype('float64')
            boot_file = boot_file.round(4)
            boot_file = boot_file[["Abundance", 'Unique_best_reads', 'Shared_best_reads', "Similarity"]]
            boot_file = boot_file[boot_file["Abundance"] > 0.05]
            boot_file.insert(0, "Cluster", boot_file.index)
            boot_file.reset_index(drop=True, inplace=True)
            boot_file['Unique_best_reads'] = boot_file['Unique_best_reads'].astype(int)
            boot_file['Shared_best_reads'] = boot_file['Shared_best_reads'].astype(int)
            # print(boot_file)
            boot_file.to_csv(output + '/_MIST_strain/_MIST_{}_measure.csv'.format(Similarity_all[i]))
            print("Level {}, done\n".format(Similarity_all[i]))
    print("All done \n")
    pass

def calculate_statistics(data):
    """计算每个菌株的丰度的平均值、标准差和95%置信区间"""
    results = {}

    for strain in data:
        abundances = data[strain]
        if len(abundances) == 0:
            continue
        abundances = np.array(abundances)

        # 计算平均值和标准差
        mean_abundance = np.mean(abundances)
        std_abundance = np.std(abundances, ddof=1)

        # 计算95%置信区间
        n = len(abundances)
        if n > 1:
            conf_interval = stats.t.interval(0.95, n - 1, loc=mean_abundance, scale=std_abundance / np.sqrt(n))
            ci_low, ci_high = conf_interval
        else:
            ci_low = ci_high = np.nan

        results[strain] = {
            'Average': mean_abundance,
            'Standard Error': std_abundance / np.sqrt(n),
            '95% CI Lower': ci_low,
            '95% CI Upper': ci_high
        }

    return results


def process_files(directory, suffix):
    """处理指定目录下的所有以suffix结尾的文件"""
    results = defaultdict(list)

    # 遍历目录中的所有文件夹
    for root, dirs, files in os.walk(directory):
       # print(directory)
        for file in files:
            if file.endswith(suffix):
                file_path = os.path.join(root, file)
                # print(f"Processing file: {file_path}")
                df = pd.read_csv(file_path, index_col=0)

                for strain in df.index:
                    abundance = df.loc[strain, 'Abundance'] if 'Abundance' in df.columns else 0
                    results[strain].append(abundance)
    #print(results)
    """将字典中每个值的长度填充为一致，较短的列表用0填充"""
    max_length = max(len(lst) for lst in results.values())

    for key in results:
        if len(results[key]) < max_length:
            results[key].extend([0] * (max_length - len(results[key])))
    return results


def bootstrap(cluster_output,mismatch, output_dir,read_length,bootstrap_numbers):
    
    # 获取文件名（不带路径和扩展名）
    base_file_name = os.path.splitext(os.path.basename(mismatch))[0]
    # 创建存放bootstrap样本的目录
    output_dir = output_dir +"/_MIST_bootstrap"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
    # 读取CSV文件并将第一列设为索引
    df = pd.read_csv(mismatch, index_col=0)
    # 获取总行数
    total_rows = len(df)
    # 计算抽样数量
    sample_size = int(total_rows / 2) + 1
    # 生成200个bootstrap样本并保存为新的CSV文件
    for i in range(1, bootstrap_numbers):
        # 进行bootstrap抽样
        bootstrap_sample = df.sample(n=sample_size, replace=True)
        # 添加标题行
        bootstrap_sample = pd.concat([df.head(0), bootstrap_sample])
        # 生成新的文件名
        new_file_name = f'{base_file_name}_{i}.csv'
        # 保存为新的CSV文件
        bootstrap_sample.to_csv(output_dir + "/" + new_file_name)

        print(f'Generated {new_file_name}')
        print(f'Run strain {new_file_name}')
        strain(cluster_output=cluster_output,
               mismatch_matrix=output_dir + "/" + new_file_name,
               read_length=read_length,
               output=output_dir + "/"+f'{base_file_name}_{i}')


    #suffixes = ['_0.995_measure.csv', '_0.9995_measure.csv', '_0.999_measure.csv', '_0.99_measure.csv']
    for root, dirs, files in os.walk(output_dir):
        suffixes = files
    for suffix in suffixes:
        data = process_files(output_dir, suffix)
        stats_results = calculate_statistics(data)

        # 生成对应的结果文件
        output_file = output_dir+f"/results_statistics{suffix.replace('.csv', '.csv')}"
        stats_df = pd.DataFrame(stats_results).T
         # 将索引重置为列，并重命名列名
        stats_df.reset_index(inplace=True)
        stats_df.rename(columns={'index': 'Cluster'}, inplace=True)
        print(stats_df)
        stats_df.to_csv(output_file)
        print(f"Results for {suffix} saved to {output_file}")
