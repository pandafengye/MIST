import os
import warnings
warnings.filterwarnings("ignore")

def species(threads, pair_1, pair_2, database, output):
    if not os.path.exists(output + "/_MIST_species/alignment/"):
        os.makedirs(output + "/_MIST_species/alignment/")
    for i, j, k in os.walk(database):
        # print(k[0].split('.')[0]plit('.')[0]
        index = k[0].split('.')[0]
    _ROOT = os.path.abspath(os.path.dirname(__file__))
    _ROOT = _ROOT.split('/')
    del (_ROOT[-1])
    str1 = '/'
    _ROOT = str1.join(_ROOT)
    bowtie2 = os.path.join(_ROOT, "bowtie2-2.4.1/bowtie2")
    os.system("{} -p {} -x {} -1 {} -2 {} --no-unal -S {}".format(bowtie2, threads, database + '/' + index, pair_1, pair_2,
                                                        output + "/_MIST_species/alignment/" + index + ".sam"))
   # os.system("sed -e '/*/d' {} > {}".format(output+"/_MIST_species/alignment/"+index + ".sam",output+"/_MIST_species/alignment/"+index + ".sam"))
    key = []
    alignment_filePath = output + "/_MIST_species/alignment/"
    files = os.listdir(alignment_filePath)  # 得到文件夹下的所有文件名称
    for file in files:  # 遍历文件
        f = open(alignment_filePath + "/" + file, 'r');  # 打开文件
        value = [];seq = {}
        number = 0
        for line in f:  # f.decode("utf8","ignore")
            li = line.strip()
            if not li.startswith("@"):  # 跳过@前三行
                a = li.split()  # 以tab键划分
                b = a[2:3];
                c = b[0];
                d = a[9:11];
                e = a[0:2];
                key.append(c)
                d = (d[0] + " " + d[1]).split("\n")
                e = (e[0] + "_" + e[1]).split()
                f = dict(zip(e, d))
                #print(f)
                if c not in seq:
                    seq[c] = f
                else:
                    seq[c] = {**seq[c], **f}
       # print(seq.keys())
        for k, v in seq.items():
            for k, v in seq.items():
                if k != "*":
        #            print(k)
                    fw = open(output + "/_MIST_species/" + "_MIST." + k + ".fq", 'w')
                    for i, j in seq[k].items():
                                # print(j)
                                # print(type(j))
                        m = j.split()
                        fw.write("@" + i + "\n")
                        fw.write(m[0] + "\n")
                        fw.write("+" + "\n")
                        fw.write(m[1])
                        fw.write("\n")
                    fw.close()
                
    DICT = dict.fromkeys(key, 0)
    #print(DICT)
    for i in key:
        DICT[i] += 1
   # DICT.pop("*")
    # print(DICT)
    # print(seq)
    import pandas as pd
    data = pd.DataFrame(list(DICT.items()), columns=['species', 'read_count'])
    # data=pd.DataFrame.from_dict(DICT,orient="index",columns=["read_count"])
    # data=data.reset_index().rename(columns={'index':'species'})
    print(data)
    data.to_csv(output + "/_MIST_species/species_count.txt")
    data1 = dict((k, v) for k, v in DICT.items() if v < 100)
    for k, v in DICT.items():
        if v < 100:
            os.remove(output + "/_MIST_species/" + "_MIST." + k + ".fq")
