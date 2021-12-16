import os
def index(refdir,output):
#    if not os.path.exists(output):
#        os.makedirs(output)
#    _ROOT = os.path.abspath(os.path.dirname(__file__))
#    _ROOT=_ROOT.split('/')
#    del(_ROOT[-1])
#    str1='/'
#    _ROOT=str1.join(_ROOT)
   # bowtie2build=os.path.join(_ROOT,"bowtie2-2.4.1/bowtie2-build")
    ref_list=os.listdir(refdir)
#    if not os.path.exists(output+"_MIST_index/"):
#        os.makedirs(output+"_MIST_index/")
    for i in ref_list:
        os.makedirs(output+"_MIST_index/"+i)
        os.system("bowtie2-build {} {}".format(refdir+i,output+"_MIST_index/"+i+"/"+i))
