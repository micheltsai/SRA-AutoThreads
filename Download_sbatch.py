import datetime
import multiprocessing
import os
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path

import numpy
import pandas as pd

import utils_
def run_cmd2(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True, check=True)
    return p.stdout.decode().strip("\n")

def getProgramTime():
    cmd="squeue -u linsslab01"
    str = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True,
                         check=True).stdout.decode().strip("\n")
    str=str.split("\n")
    try:
        print(str)
        print(len(str))
        for x in range(1,len(str)):
            s=str[x].strip().split("   ")
            s=[x for x in s if x!='']
            print(s)
            programID=s[0]
            print(programID)
            time=s[2]
            time=time.split(":")
            print(time)
            if len(time)>=3:
                print("len(time)>2")
                if int(time[0])>=2:
                    print("int(time[0])>1")
                    print("time: {}\n".format(time))
                    print("time[0]: {}\n".format(time[0]))
                    print("scancel -b {}\n".format(programID))
                    run_cmd2("scancel -b {}".format(programID))
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
    return 0
def Download(x,_outdir,sra_dir):
    one_ = time.time()
    #print(
    #   "---------------------\n---------------------[ {} / {} ]---------------------\n".format(num + 1,
    #                                                                                           len(idlist)))
    #num += 1
    print("x = {}".format(x))
    # outdir__ = os.path.join(output, "out")
    outdir__ = os.path.join(sra_dir, "Assembled")

    check_dir = os.path.join(_outdir, "check")
    utils_.mkdir_join(check_dir)
    check_log = os.path.join(check_dir, "Analysischeck.log")

    final_dir = os.path.join(outdir__, "{}_contig.fa".format(x))
    utils_.mkdir_join(sra_dir)
    sra_file=os.path.join(sra_dir,"sra/{}/{}.sra".format(x,x))
    print("final_dir:",final_dir)
    print("sra_file:",sra_file)
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
    elif os.path.isfile(sra_file):
        print("was ran download ,sra is exist\n------------------------------\n\n")
    else:

        sra_dir = os.path.join(sra_dir, "sra")
        utils_.mkdir_join(sra_dir)
        print("sra_dir:",sra_dir)
        utils_.prefetch_srav2(x, sra_dir)
        print("Download {}\n.".format(x))
        #with open(Downloadcheck_log, "a+") as f:
        #    f.write("{}\n".format(x))
    dltime=time.time() - one_
    print('Done,total cost',dltime, 'secs')
    print("###########################################################")

def sbatch_job(outdir,pdat,need_list,ll,limit_number,sra_num_,finish_,start):
    print("outdir:{}\nnew_outdir:{}\n".format(outdir, outdir))

    check_dir = os.path.join(outdir, "check")
    utils_.mkdir_join(check_dir)

    job_dir = os.path.join(check_dir, "job")
    utils_.mkdir_join(job_dir)
    job_file = os.path.join(job_dir, "{}.sh".format("pdat"))

    job_out = os.path.join(check_dir, "JOBoutput")
    utils_.mkdir_join(job_out)

    check_log = os.path.join(check_dir, "Analysischeck.log")
    check_file = Path(check_log)
    check_file.touch(exist_ok=True)

    sraList = os.path.join(outdir, "sraList_test.txt")

    needList = os.path.join(check_dir, "need_run_{}.txt".format(ll+len(need_list)))
    need_file = Path(needList)
    need_file.touch(exist_ok=True)
    with open(needList,"w")as f:
        f.write("")
    #with open(sraList, "r") as f:
    #    run_list = f.readlines()
   # print("run_list: ", run_list)

    # run_list = run_list[0].split("\n")
    # print("run_list: ",run_list)
    #run_list = [rr.strip() for rr in run_list if rr.strip() != '']


    # print(run_list)
    f = open(check_log, 'r')
    line_Analysis = f.readlines()
    #print("check log :{}\n".format(line_Analysis))
    f.close()
    #for s in line_Analysis:
        #print("{}\n".format(s))
    #finish_Analysis = list(filter(lambda x: len(x.split(" ")) >= 4, line_Analysis))
    #finish_Analysis_run = list(map(lambda x: x.split(" ")[1], finish_Analysis))
    #need_run = list(filter(lambda x: x not in finish_Analysis_run, need_list))
    #print("finish: {}\nfinish_run: {}\nneed_run: {}".format(finish_Analysis, finish_Analysis_run, need_run))
    #print(
    #    "finish length: {}\nfinish_run length: {}\nneed_run length: {}".format(len(finish_Analysis), len(finish_Analysis_run),
    #                                                                           len(need_run)))
    #sra_num_ = len(need_run) + len(finish_Analysis_run)
    finish_num = finish_
    #print("len(finish_run)+len(need_run) = {}".format(sra_num_))
    #num = len(finish_Analysis_run)
    progress_list = []
    prog_num = 0

    #finish_num = len(finish_Analysis_run)
    #finish_num_ = len(finish_Analysis_run)
    print("finish_num = {}\n".format(finish_num))
    print("need_num={}\n".format(len(need_list)-finish_num))
    pool_list = []
    SRA_run_error = os.path.join(check_dir, "SRA_run_error.txt")
    try:
        for k in need_list:
            k.strip("\n")
            print("########### hello %d ############\n" % prog_num)
            print(k)
            print(need_list.index(k))
            print("########## {}/{} ###########".format(finish_num+1+ll, sra_num_))
            # utils_.run_cmd("sbatch -A MST109178 -J Job_test -p ngs48G -c 14 --mem=46g -o ./out/{}_array_out.log -e ./out/{}_array_out.log "
            #               "--mail-user=sj985517@gmail.com --mail-type=BEGIN,END --wrap='/home/linsslab01/miniconda/bin/python3 one_Analysis.py'--array=1-4")
            with open(needList, "a+") as f:
                f.write("{}\n".format(k))

            prog_num += 1
            finish_num += 1
            time.sleep(1)
        job_out_o = os.path.join(job_out, "%A_test_%a.log")
        job_out_err = os.path.join(job_out, "%A_err_%a.log")
        with open(job_file, "w+") as f:
            f.write("#!/usr/bin/sh\n")
            f.write("#SBATCH -A MST109178\n")
            f.write("#SBATCH -J {}_job\n".format("test"))
            f.write("#SBATCH -p ngs7G\n")
            f.write("#SBATCH -c 2\n")
            f.write("#SBATCH --mem 7g\n")
            f.write("#SBATCH --array=1-{}\n".format(len((need_list))))
            f.write("#SBATCH -o {}\n".format(job_out_o))
            f.write("#SBATCH -e {}\n".format(job_out_err))
            f.write("#SBATCH --mail-user=sj985517@gmail.com\n")
            f.write("#SBATCH --mail-type=BEGIN,END\n")
            f.write("echo $SLURM_ARRAY_TASK_ID\n")
            f.write("/home/linsslab01/miniconda3/bin/python3 one_Analysis_new.py {} {} {} {}\n".format(pdat,
                                                                                                sra_num_,
                                                                                                "$SLURM_ARRAY_TASK_ID",
                                                                                                       ll+len(need_list)))

        print("sbatch {}".format(job_file))
        utils_.run_cmd("sbatch {}".format(job_file))
        time.sleep(2)
        utils_.run_cmd("squeue -u linsslab01")

    except KeyboardInterrupt:
        print("Catch keyboardinterdinterupterror\n")
        print("srart : {}\n".format(start))
        print("Download all ", 'Done,total cost', time.time() - start, 'secs')
        pid = os.getgid()

        with open(SRA_run_error, "a+") as f:
            f.write("Catch keyboardinterdinterupterror\n")
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

        with open(SRA_run_error, "a+") as f:
            f.write("{}\n".format(errMsg))
        sys.exit(errMsg)

def pool_append(need_list,outdir):
    pool = multiprocessing.Pool(processes=50)
    pool_list=[]
    for x in need_list:
        print("###################\n")
        print(x)
        new_outdir = os.path.join(outdir, "output")
        utils_.mkdir_join(new_outdir)
        print(new_outdir)
        sraid_outdir = os.path.join(new_outdir, x)
        utils_.mkdir_join(sraid_outdir)
        print(sraid_outdir)
        # sraid_outdir = os.path.join(sraid_outdir, "sra")
        # utils_.mkdir_join(sraid_outdir)
        print(sraid_outdir)
        # Download(x,outdir,sraid_outdir)
        pool_list.append(pool.apply_async(Download, (x, outdir, sraid_outdir,)))
        #pdat = pdat_run[need_run.index(x)]
    pool.close()
    print("pool.close()")
    pool.join()
    print("pool.join()")


def main():
    pool = multiprocessing.Pool(processes=50)
    start = time.time()
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    ##########
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
    today1=str(datetime.datetime.now())
    ##################
    print("pdat: {}".format(today1))
    sraList_txt=os.path.join(outdir,"sraList_test.txt")
    myfile = Path(sraList_txt)
    myfile.touch(exist_ok=True)
    with open(sraList_txt,"r")as f:
        lines=f.readlines()

    #######<DATE>:Run <sraID> is ok.
    #sralist = list(filter(lambda x: len(x.split(" ")) >= 4, lines))
    #sra_run = list(map(lambda x: x.split(" ")[1], sralist))

    ###<DATE>:<sraID>
    sralist = list(filter(lambda x: len(x.split(":")) >= 2, lines))
    sra_run = list(map(lambda x: x.split(":")[1].strip("\n"), sralist))
    #######
    pdat_run=list(map(lambda x: x.split(" ")[0].split(":")[0], sralist))
    Year=pdat_run[0].split("/")[0]
    Month = pdat_run[0].split("/")[1]
    #print(sra_run)

    ##build analysis_final
    check_dir = os.path.join(outdir, "check")
    utils_.mkdir_join(check_dir)
    final_log= os.path.join(check_dir, "analysis_final.csv")
    if os.path.isfile(final_log):
        print("{} is exist.\n".format(final_log))
    else:
        myfile_final = Path(final_log)
        myfile_final.touch(exist_ok=True)
        print("build {}.\n".format(final_log))
        with open(final_log,"r+") as f:
            df_clowns = ""+","+"Accession"+","+"MLST"+","+"AMR"+","+"Point"+","+"Serotype"+","+"IncType"
            f.write(df_clowns+"\n")
            print(f.readlines())
        print("build {} Done.\n".format(final_log))
    ######

    check_log = os.path.join(check_dir, "Analysischeck.log")
    myfile2 = Path(check_log)
    myfile2.touch(exist_ok=True)
    with open(check_log,"r")as f:
        line=f.readlines()

    finish = list(filter(lambda x: len(x.split(" ")) >= 4, line))
    finish_run = list(map(lambda x: x.split(" ")[1], finish))
    need_run = list(filter(lambda x: x not in finish_run, sra_run))

    #####sra_id file not fill ANI>95 in nofillQC.txt
    QCcheck_log = os.path.join(check_dir, "nofillQC.txt")
    QC_file = Path(QCcheck_log)
    QC_file.touch(exist_ok=True)
    with open(QCcheck_log,"r")as f:
        QCline=f.readlines()
    QC_run=list(filter(lambda x: len(x.split(":")) >= 2, QCline))
    QCfinish_run = list(map(lambda x: x.split(" ")[1], QC_run))
    need_run = list(filter(lambda x: x not in QCfinish_run, need_run))

    print("sra_list length: {}\n".format(len(sra_run)))
    print("need_run length: {}\n".format(len(need_run)))
    print("finish_list length: {}\n".format(len(finish_run)))

    pool_list=[]
    pdat=""

    limit_list=list(range(0, len(need_run), limit_num))




    for ll in limit_list:


        need_list=need_run[ll:ll+limit_num]
        print("###############\n{} -> {}\n".format(ll+len(finish_run),len(need_list)+len(finish_run)))
        print("Download\n")
        pool_append(need_list,outdir)
        print("Download End\n")

        # git date for store gz file name
        pdat = pdat_run[need_run.index(need_list[0])]
        print("pdat:{}".format(pdat))


        print("################\nsbatch_job\n")
        sbatch_job(outdir,pdat,need_list,ll,limit_num,len(sra_run),len(finish_run),start)

        print("sbatch_job {}->{}\n".format(ll,ll+limit_num))
        #print(run_cmd2("squeue -u linsslab01 |wc -l"))
        #print(type(run_cmd2("csqueue -u linsslab01 |wc -l")))

        utils_.run_cmd2("squeue -u linsslab01")
        num=int(run_cmd2("squeue -u linsslab01 |wc -l"))
        # pd_start=time.time()
        # while num == 2:
        #     print("progresses status is PD, or one progress is running\n")
        #     try:
        #         #utils_.run_cmd2("squeue -u linsslab01")
        #         #time.sleep(2)
        #         tmp= run_cmd2("squeue -u linsslab01 |wc -l")
        #         num = int(tmp)
        #         time.sleep(60)
        #     except Exception as e:
        #         print(tmp)
        #         #print("again run 'squeue -u linsslab01'\n")
        #         #utils_.run_cmd2("squeue -u linsslab01")
        #         #time.sleep(2)
        #         ####
        #         print("print 'squeue -u linsslab01' error:\n")
        #         error_class = e.__class__.__name__  # 取得錯誤類型
        #         detail = e.args[0]  # 取得詳細內容
        #         cl, exc, tb = sys.exc_info()  # 取得Call Stack
        #         lastCallStack = traceback.extract_tb(tb)[-1]  # 取得Call Stack的最後一筆資料
        #         fileName = lastCallStack[0]  # 取得發生的檔案名稱
        #         lineNum = lastCallStack[1]  # 取得發生的行號
        #         funcName = lastCallStack[2]  # 取得發生的函數名稱
        #         errMsg = "File \"{}\", line {}, in {}: [{}] {}".format(fileName, lineNum, funcName, error_class,
        #                                                                detail)
        #         print(errMsg)
        #         time.sleep(5)
        #         pass
        #
        #
        # print(str(datetime.datetime.now()), 'PD Done,current total cost', time.time() - pd_start, 'secs\n')

        running_start=time.time()
        while num != 1:
            print("progresses is running\n")
            try:
                #utils_.run_cmd2("squeue -u linsslab01")
                #time.sleep(2)
                tmp=run_cmd2("squeue -u linsslab01 |wc -l")
                num=int(tmp)
                print("Quantity of running progress  = {}\n".format(num - 1))
                print(str(datetime.datetime.now()), 'Running,current total cost', time.time() - running_start, 'secs\n')
                if time.time() - running_start >= 7200:
                    print("getProgramTime()\n")
                    getProgramTime()
                if num == 1:
                    break
                time.sleep(60)
            except Exception as e:
                print(tmp)
                #print("again run 'squeue -u linsslab01'\n")
                #utils_.run_cmd2("squeue -u linsslab01")
                #time.sleep(2)
                ####
                print("print 'squeue -u linsslab01' error:\n")
                error_class = e.__class__.__name__  # 取得錯誤類型
                detail = e.args[0]  # 取得詳細內容
                cl, exc, tb = sys.exc_info()  # 取得Call Stack
                lastCallStack = traceback.extract_tb(tb)[-1]  # 取得Call Stack的最後一筆資料
                fileName = lastCallStack[0]  # 取得發生的檔案名稱
                lineNum = lastCallStack[1]  # 取得發生的行號
                funcName = lastCallStack[2]  # 取得發生的函數名稱
                errMsg = "File \"{}\", line {}, in {}: [{}]_ {}".format(fileName, lineNum, funcName, error_class,
                                                                       detail)
                print(errMsg)
                time.sleep(5)
                pass


        print(str(datetime.datetime.now()), 'sbatch Done,current total cost', time.time() - running_start, 'secs\n')

        output_dir = os.path.join(outdir, "output")

        joboutput_dir = os.path.join(check_dir, "JOBoutput")
        print("outputdir: {}\n".format(output_dir))
        print("JOBoutputdir: {}\n".format(joboutput_dir))
        ##################

        needList = os.path.join(check_dir, "need_run_{}.txt".format(ll))
        print("rm -rf {}\n".format(needList))
        utils_.run_cmd2("rm -rf {}".format(needList))
        #################
        mytarfile = os.path.join(outdir, "{}_{}.tar.gz".format(str(ed_M), ll + len(need_list) + len(finish_run)))
        print("mytarfile: {}\n".format(mytarfile))
        #tar.gz
        tar_start = time.time()
        print("tar zcvf {} {}\n".format(mytarfile,output_dir.replace(current_path,".")))

        utils_.run_cmd("tar zcvf {} {}".format(mytarfile,output_dir.replace(current_path,".")))
        print(str(datetime.datetime.now()), 'tar Done,current total cost', time.time() - tar_start, 'secs\n')
        time.sleep(1)
        ###################

        scp_start=time.time()
        print("scp -r {} root@140.112.165.124:/data/SRA_data/{}/{}/output/{}_{}.tar.gz\n".format(mytarfile,str(ed_Y),str(ed_M),str(ed_M),ll+len(need_list)+len(finish_run)))
        utils_.run_cmd("scp -r {} root@140.112.165.124:/data/SRA_data/{}/{}/output/{}_{}.tar.gz".format(mytarfile,str(ed_Y),str(ed_M),str(ed_M),ll+len(need_list)+len(finish_run)))
        print(str(datetime.datetime.now()), 'scp Done,current total cost', time.time() - scp_start, 'secs\n')
        time.sleep(1)
        ##################

        #output
        remove_start=time.time()
        print("rm -rf {}\n".format(output_dir))
        #utils_.run_cmd("rm -rf SRAtest/output")
        utils_.run_cmd("rm -rf {}".format(output_dir))
        #shutil.rmtree("./SRAtest/output")
        print("rm -rf {}\n".format(joboutput_dir))
        utils_.run_cmd("rm -rf {}".format(joboutput_dir))
        print("rm -rf {}\n".format(mytarfile))
        utils_.run_cmd("rm -rf {}".format(mytarfile))

        print(str(datetime.datetime.now()), 'remove Done,current total cost', time.time() - remove_start, 'secs\n')
        ##################
        time.sleep(1)

    print("analysislog: {}\n".format(check_log))

    #######
    print("scp {} root@140.112.165.124:/data/SRA_data/\n".format(final_log,pdat_run[0]))


    # print("scp {} root@140.112.165.124:/data/SRA_data/{}/{}\n".format(final_log, str(ed_Y), str(ed_M)))
    # utils_.run_cmd("scp {} root@140.112.165.124:/data/SRA_data/{}/{}".format(final_log, str(ed_Y), str(ed_M)))
    # print("scp {} root@140.112.165.124:/data/SRA_data/{}/{}\n".format(check_log,str(ed_Y),str(ed_M)))
    # utils_.run_cmd("scp {} root@140.112.165.124:/data/SRA_data/{}/{}".format(check_log,str(ed_Y),str(ed_M)))
    # print("scp {} root@140.112.165.124:/data/SRA_data/{}\n".format(sraList_txt,str(ed_Y),str(ed_M)))
    # utils_.run_cmd("scp {} root@140.112.165.124:/data/SRA_data/{}/{}".format(sraList_txt,str(ed_Y),str(ed_M)))
    # print("scp {} root@140.112.165.124:/data/SRA_data/{}\n".format(QCcheck_log, str(ed_Y), str(ed_M)))
    # utils_.run_cmd("scp {} root@140.112.165.124:/data/SRA_data/{}/{}".format(QCcheck_log, str(ed_Y), str(ed_M)))
    print("Progerss end\n")
    print(str(datetime.datetime.now()), ' Done,current total cost', time.time() - start, 'secs\n')
if __name__ == '__main__':
    main()