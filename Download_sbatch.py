import multiprocessing
import os
import time
from pathlib import Path

import pandas as pd

import utils_

def Download(x,_outdir,sra_dir):
    one_ = time.time()
    #print(
    #   "---------------------\n---------------------[ {} / {} ]---------------------\n".format(num + 1,
    #                                                                                           len(idlist)))
    #num += 1
    print("x = {}".format(x))
    # outdir__ = os.path.join(output, "out")
    outdir__ = os.path.join(_outdir, "Assembled")
    check_log = os.path.join(_outdir, "Analysischeck.log")
    final_dir = os.path.join(outdir__, "{}_contig.fa".format(x))
    sra_file=os.path.join(sra_dir,"{}/{}.sra".format(x,x))
    print(final_dir)
    print(sra_file)
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
    elif os.path.isfile(sra_file):
        print("was ran download ,sra is exist\n------------------------------\n\n")
    else:
        utils_.prefetch_srav2(x, sra_dir)
        print("Download {}\n.".format(x))
        #with open(Downloadcheck_log, "a+") as f:
        #    f.write("{}\n".format(x))
    dltime=time.time() - one_
    print('Done,total cost',dltime, 'secs')
    print("###########################################################")

def main():
    pool = multiprocessing.Pool(processes=50)
    start = time.time()
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    ## read SRAsetting.txt
    setting_path = os.path.join(current_path, "SRAsettings.txt")

    with open(setting_path, "r") as f:
        setList = f.readlines()

    print(setList)
    i = 0
    settings_dict = {}
    for line in setList:
        line = line.strip("\n")
        line_ = line.split("=")
        if line != "" and len(line_) == 2:
            print(line_)
            print("line{}. {}:{}\n".format(i, line_[0], line_[1]))
            settings_dict.update({line_[0]: line_[1]})
        i += 1
    print(settings_dict)
    # setting_df=pd.DataFrame(settings_dict)
    setting_df = pd.DataFrame.from_dict(settings_dict, orient='index').T
    print(setting_df.columns)

    start_date = str(setting_df['start_date'][0])
    expiry_date = str(setting_df['expiry_date'][0])
    thread = int(setting_df['cpu_thread'][0])
    cpu_process = int(setting_df['process'][0])
    gsize = str(setting_df['gsize'][0])
    outdir = str(setting_df['output_dir'][0])
    ref_dir = str(setting_df['Busco_ReferenceSequenceFileDir_Path'][0])
    buscoDB = str(setting_df['Busco_database'][0])
    buscoMode = str(setting_df['Busco_mode'][0])
    mlstS = str(setting_df['MLST_organism'][0])
    amrS = str(setting_df['AMR_organism'][0])
    shovill_RAM = str(setting_df['shovill_RAM'][0])
    # get (Date) to (Date)
    sd_Y = int(start_date.split("/")[0])
    sd_M = int(start_date.split("/")[1])
    sd_D = int(start_date.split("/")[2])
    ed_Y = int(expiry_date.split("/")[0])
    ed_M = int(expiry_date.split("/")[1])
    ed_D = int(expiry_date.split("/")[2])
    print(sd_Y, sd_M, sd_D)
    print(ed_Y, ed_M, ed_D)
    utils_.mkdir_join(outdir)

    sraList_txt=os.path.join(outdir,"sraList_test.txt")
    myfile = Path(sraList_txt)
    myfile.touch(exist_ok=True)
    with open(sraList_txt,"r")as f:
        lines=f.readlines()
    sralist = list(filter(lambda x: len(x.split(" ")) >= 4, lines))
    sra_run = list(map(lambda x: x.split(" ")[1], sralist))
    print(sra_run)
    pool_list=[]
    for x in sra_run:
        print("###################\n")
        print(x)
        sraid_outdir=os.path.join(outdir,x)
        utils_.mkdir_join(sraid_outdir)
        Download(x,outdir,sraid_outdir)
        pool_list.append(pool.apply_async(Download, (x,outdir,sraid_outdir,)))
    pool.close()
    print("pool.close()")
    pool.join()
    print("pool.join()")
if __name__ == '__main__':
    main()