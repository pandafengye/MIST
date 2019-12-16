import os,re,math
from collections import Counter
def map(indexpath,read_length,output,thread=8,single_end=None,pair_1=None,pair_2=None):
    if not os.path.exists(output + "/_MIST_map_alignment/"):
        os.mkdir(output+ "/_MIST_map_alignment/")
    """Map the fastq file to the reference genome using bowtie2"""
    index_list=os.listdir(indexpath)
    if pair_1==None and pair_2==None and single_end!=None:
        for i in index_list:
            os.system("bowtie2 --local  -p {} -x {} -U {} -S {}".format(thread,indexpath + i + "/" + i, single_end,
                                                                             output + "/_MIST_map_alignment/" + i + ".fq.sam"))
    if pair_1!=None and pair_2!=None and single_end==None:
        os.system("cat {} {} > {} ".format(pair_1, pair_2,output+"/_MIST_map_merge.fq"))
        for i in index_list:
            os.system("bowtie2 --local  -p {} -x {} -U {} -S {}".format(thread,indexpath + i + "/" + i, output+"/_MIST_map_merge.fq",
                                                                            output + "/_MIST_map_alignment/" + i + ".fq.sam"))
    os.system("rm -rf {}/_MIST_map_merge.fq".format(output))
    """Statistics Mismatch number"""
    DICT = {}
    alignment_filePath=output + "/_MIST_map_alignment/"
    files = os.listdir(alignment_filePath)  # 得到文件夹下的所有文件名称
    for file in files:  # 遍历文件夹
        if not os.path.isdir(file):  # 判断是否是文件夹，不是文件夹才打开
            f = open(alignment_filePath + "/" + file, 'r');  # 打开文件
            result = [];
            key = [];
            value = []
            reg = r'NM:i:(.{1})'
            wordreg = re.compile(reg)
            for line in f:
                # f.decode("utf8","ignore")
                li = line.strip()
                if not li.startswith("@"):  # 跳过@前三行
                    a = li.split()  # 以tab键划分
                    b = a[0:2];
                    c1 = a[5:6];  # 提取前两列和第五列
                    c2 = re.findall(wordreg, li)  # 查找NM:i:后面的数值
                    key.append(b);
                    value.append(c1 + c2)

            Mismatch_matrix = []  # 计算缺失值
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
            Dict = {}
            for i in Mismatch_matrix:
                if list(i.keys())[0] not in Dict:
                    #     print(DICT[list(b.keys())[0]])
                    Dict[list(i.keys())[0]] = list(i.values())[0]
                else:
                    Dict[list(i.keys())[0]] = list(i.values())[0] + Dict[list(i.keys())[0]]
            Dict = dict(sorted(Dict.items(), key=lambda d: d[0]))
        DICT['.'.join(file.split(".")[0:2])] = Dict
    import pandas as pd
    df = pd.DataFrame.from_dict(DICT, orient='index')
    df=df.T
    Count=0
    for i in range(len(df)):#一行一行的循环
        if any(df.iloc[i]<10):#寻找mismatch小于10的reads，即为完全比对上
            Count=Count+1#完全比对上的条数
    print("There are a total of {} sequences, a total of {} mapped!".format(len(df),Count))
    df.to_csv(output+"/_MIST_map_Mismatch_matrix.csv")
