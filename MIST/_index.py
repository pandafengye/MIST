import os
def index(refdir,output):
    ref_list=os.listdir(refdir)
    if not os.path.exists(output+"_MIST_index/"):
        os.mkdir(output+"_MIST_index/")
    for i in ref_list:
        os.mkdir(output+"_MIST_index/"+i)
        os.system("bowtie2-build {} {}".format(refdir+i,output+"_MIST_index/"+i+"/"+i))
