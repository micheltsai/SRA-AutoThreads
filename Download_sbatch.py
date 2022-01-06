import datetime
import multiprocessing
import os
import sys
import time
import traceback
from pathlib import Path

import numpy
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
    outdir__ = os.path.join(sra_dir, "Assembled")
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
        sra_dir = os.path.join(sra_dir, "sra")
        utils_.prefetch_srav2(x, sra_dir)
        print("Download {}\n.".format(x))
        #with open(Downloadcheck_log, "a+") as f:
        #    f.write("{}\n".format(x))
    dltime=time.time() - one_
    print('Done,total cost',dltime, 'secs')
    print("###########################################################")

def sbatch_job(outdir,pdat,need_list,start):
    print("outdir:{}\nnew_outdir:{}\n".format(outdir, outdir))
    job_dir = os.path.join(outdir, "job")
    utils_.mkdir_join(job_dir)
    utils_.mkdir_join(outdir)
    job_file = os.path.join(job_dir, "{}.sh".format("test"))
    job_out = os.path.join(outdir, "JOBoutput")
    utils_.mkdir_join(job_out)

    check_log = os.path.join(outdir, "Analysischeck.log")
    check_file = Path(check_log)
    check_file.touch(exist_ok=True)
    sraList = os.path.join(outdir, "sraList_test.txt")
    needList = os.path.join(outdir, "need_run.txt")
    need_file = Path(needList)
    need_file.touch(exist_ok=True)

    #with open(sraList, "r") as f:
    #    run_list = f.readlines()
   # print("run_list: ", run_list)

    # run_list = run_list[0].split("\n")
    # print("run_list: ",run_list)
    #run_list = [rr.strip() for rr in run_list if rr.strip() != '']


    # print(run_list)
    f = open(check_log, 'r')
    line_Analysis = f.readlines()
    print("check log :{}\n".format(line_Analysis))
    f.close()
    for s in line_Analysis:
        print("{}\n".format(s))
    finish_Analysis = list(filter(lambda x: len(x.split(" ")) >= 4, line_Analysis))
    finish_Analysis_run = list(map(lambda x: x.split(" ")[1], finish_Analysis))
    need_run = list(filter(lambda x: x not in finish_Analysis_run, need_list))
    print("finish: {}\nfinish_run: {}\nneed_run: {}".format(finish_Analysis, finish_Analysis_run, need_run))
    print(
        "finish length: {}\nfinish_run length: {}\nneed_run length: {}".format(len(finish_Analysis), len(finish_Analysis_run),
                                                                               len(need_run)))
    sra_num_ = len(need_run) + len(finish_Analysis_run)
    finish_num = 0
    print("len(finish_run)+len(need_run) = {}".format(sra_num_))
    num = len(finish_Analysis_run)
    progress_list = []
    prog_num = 0

    finish_num = len(finish_Analysis_run)
    finish_num_ = len(finish_Analysis_run)
    print("finish_num = {}".format(finish_num))
    pool_list = []
    with open(needList, "w+") as f:
        f.write("")
    try:
        for k in need_run:
            k.strip("\n")
            print("########### hello %d ############\n" % prog_num)
            print(k)
            print(need_run.index(k))
            print("########## {}/{} ###########".format(finish_num+1, sra_num_))
            # utils_.run_cmd("sbatch -A MST109178 -J Job_test -p ngs48G -c 14 --mem=46g -o ./out/{}_array_out.log -e ./out/{}_array_out.log "
            #               "--mail-user=sj985517@gmail.com --mail-type=BEGIN,END --wrap='/home/linsslab01/miniconda/bin/python3 one_Analysis.py'--array=1-4")
            with open(needList, "a+") as f:
                f.write("{}\n".format(k))

            prog_num += 1
            finish_num += 1
            time.sleep(1)
        job_out_o = os.path.join(job_out, "array_test_%a.log")
        job_out_err = os.path.join(job_out, "array_err_%a.log")
        with open(job_file, "w+") as f:
            f.write("#!/usr/bin/sh\n")
            f.write("#SBATCH -A MST109178\n")
            f.write("#SBATCH -J {}_job\n".format("test"))
            f.write("#SBATCH -p ngs7G\n")
            f.write("#SBATCH -c 2\n")
            f.write("#SBATCH --mem 7g\n")
            f.write("#SBATCH --array=1-{}\n".format(len((need_run))))
            f.write("#SBATCH -o {}\n".format(job_out_o))
            f.write("#SBATCH -e {}\n".format(job_out_err))
            f.write("#SBATCH --mail-user=sj985517@gmail.com\n")
            f.write("#SBATCH --mail-type=BEGIN,END\n")
            f.write("echo $SLURM_ARRAY_TASK_ID\n")
            f.write("/home/linsslab01/miniconda3/bin/python3 one_Analysis_new.py {} {} {}\n".format(pdat,
                                                                                                sra_num_,
                                                                                                "$SLURM_ARRAY_TASK_ID"))
        ###
        print("sbatch before du-sh\n")
        print(utils_.run_cmd("du ./SRAtest -sh"))
        time.sleep(1)
        ###
        print("sbatch {}".format(job_file))
        utils_.run_cmd("sbatch {}".format(job_file))
    except KeyboardInterrupt:
        print("Catch keyboardinterdinterupterror\n")
        print("srart : {}\n".format(start))
        print("Download all ", 'Done,total cost', time.time() - start, 'secs')
        pid = os.getgid()
        with open("./SRA_run_error.txt", "a+") as f:
            f.write("Catch keyboardinterdinterupterror : {}/{}/{}\n".format())
        # with open("./Automate_check.log", "a+") as f:
        #    f.write("keyboardinterupter")
        #    f.write("{}:{}:{}\n".format(date, time.time() - ds, time.time() - start))
        # sys.exit("Catch keyboardinterdinterupterror")
        os.popen("taskkill.exe /f /pid:%d" % pid)
    except Exception as e:
        error_class = e.__class__.__name__  # 取得錯誤類型
        detail = e.args[0]  # 取得詳細內容
        cl, exc, tb = sys.exc_info()  # 取得Call Stack
        lastCallStack = traceback.extract_tb(tb)[-1]  # 取得Call Stack的最後一筆資料
        fileName = lastCallStack[0]  # 取得發生的檔案名稱
        lineNum = lastCallStack[1]  # 取得發生的行號
        funcName = lastCallStack[2]  # 取得發生的函數名稱
        errMsg = "File \"{}\", line {}, in {}: [{}] {}".format(fileName, lineNum, funcName, error_class,
                                                               detail)
        print(errMsg)
        ###
        # process = psutil.Process(os.getpid())
        # print(str(datetime.datetime.now()), process.memory_info().rss)
        utils_.run_cmd("free -h > ./checkmem.txt")
        ####
        with open("./SRA_run_error.txt", "a+") as f:
            f.write("{}\n".format(errMsg))
        sys.exit(errMsg)




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
    limit_num=int(str(setting_df['limit_number'][0]))
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
    pdat=str(datetime.datetime.now())

    print("pdat: {}".format(pdat))
    sraList_txt=os.path.join(outdir,"sraList_test.txt")
    myfile = Path(sraList_txt)
    myfile.touch(exist_ok=True)
    with open(sraList_txt,"r")as f:
        lines=f.readlines()
    sralist = list(filter(lambda x: len(x.split(" ")) >= 4, lines))
    sra_run = list(map(lambda x: x.split(" ")[1], sralist))
    pdat_run=list(map(lambda x: x.split(" ")[0].split(":")[0], sralist))
    print(sra_run)
    pool_list=[]
    pdat=""

    limit_list=list(range(0, len(sra_run), limit_num))



    for ll in limit_list:
        need_list=sra_run[ll:ll+limit_num]
        print("###############\n{} -> {}\n".format(ll,ll+limit_num))
        for x in need_list:
            print("###################\n")
            print(x)
            new_outdir = os.path.join(outdir, "output")
            sraid_outdir=os.path.join(new_outdir,x)
            utils_.mkdir_join(sraid_outdir)
            #Download(x,outdir,sraid_outdir)
            pool_list.append(pool.apply_async(Download, (x,outdir,sraid_outdir,)))
            pdat=pdat_run[sra_run.index(x)]
        sbatch_job(outdir,pdat,need_list,start)



    pool.close()
    print("pool.close()")
    pool.join()
    print("pool.join()")

    print(str(datetime.datetime.now()), ' Done,current total cost', time.time() - start, 'secs\n')
if __name__ == '__main__':
    main()