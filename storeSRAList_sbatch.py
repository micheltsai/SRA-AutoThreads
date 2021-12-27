import datetime
import multiprocessing
import os
import sys
import time
import traceback
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

def sra_stat_old(sra_id,outdir,sra_dir,isfinal):

    QC_error = os.path.join(outdir, "nofillQC.txt")

    print("SequenceReadArchive\n")
    sra = utils_.SequenceReadArchivev3(sra_id)
    _base_ = sra.base_percentage() * 100
    print("base percentage: ", _base_, "\n")
    #######Q30 base>=80%
    if _base_ < 80:
        # shutil.rmtree(outdir)
        with open(QC_error, "a+") as f:
            f.write("{}: Reads quality is too low\n".format(sra_id))
        sys.exit('Reads quality is too low.\n')
    ###### layout = 2

    if sra.layout != '2':
        with open(QC_error, "a+") as f:
            f.write("{}: File layout is not pair-end\n".format(sra_id))
        sys.exit(f'File layout is not pair-end\n')

    print("layout=2\n")

    # if sra_layout==2 continue
    Download(sra_id, outdir, sra_dir)
    sraList = os.path.join(outdir, "sraList.txt")
    with open(sraList, "a+") as f:
        f.write(sra_id)
        if isfinal == False:
            f.write(",")
    with open(sraList, "r") as f:
        print(f.readlines())


def sra_stat(sra_id, outdir, sra_dir, sraNUM, needNUM, date):
    QC_error = os.path.join(outdir, "nofillQC.txt")
    global fin_number
    print("SequenceReadArchive\n")
    sra = utils_.SequenceReadArchivev3(sra_id)
    _base_ = sra.base_percentage() * 100
    print("base percentage: ", _base_, "\n")
    #######Q30 base>=80%
    if _base_ < 80:
        # shutil.rmtree(outdir)
        with open(QC_error, "a+") as f:
            f.write("{}: Reads quality is too low\n".format(sra_id))
        fin_number+=1
        sys.exit('Reads quality is too low.\n')
    else:
        ###### layout = 2
        if sra.layout != '2':
            with open(QC_error, "a+") as f:
                f.write("{}: File layout is not pair-end\n".format(sra_id))
            fin_number+=1
            sys.exit(f'File layout is not pair-end\n')
        else:
            print("layout=2\n")
            # if sra_layout==2 continue
            Download(sra_id, outdir, sra_dir)
            sraList = os.path.join(outdir, "sraList.txt")
            with open(sraList, "a+") as f:
                f.write("Run {} is ok.\n".format(sra_id))
            with open(sraList, "r") as f:
                sraL = f.readlines()
                print("{}/{}: {}\n".format(fin_number, sraNUM, sraL))

            if needNUM == fin_number:
                print("{} store SRAList End.\n".format(date))
                # sys.exit("{} store SRAList End.\n".format(date))
            fin_number += 1

def sbatch_job(outdir,pdat):
    new_outdir = os.path.join(outdir, pdat)
    print("outdir:{}\nnew_outdir:{}\n".format(outdir, new_outdir))
    job_dir = os.path.join(outdir, "job")
    utils_.mkdir_join(job_dir)
    utils_.mkdir_join(new_outdir)
    job_file = os.path.join(job_dir, "{}.sh".format(pdat))
    job_out = os.path.join(outdir, "JOBoutput")
    utils_.mkdir_join(job_out)
    job_out = os.path.join(outdir, pdat)
    utils_.mkdir_join(job_out)

    check_log = os.path.join(new_outdir, "Analysischeck.log")

    sraList = os.path.join(new_outdir, "sraList.txt")
    needList = os.path.join(new_outdir, "need_run.txt")

    with open(sraList, "r") as f:
        run_list = f.readlines()
    print("run_list: ", run_list)

    # run_list = run_list[0].split("\n")
    # print("run_list: ",run_list)
    run_list = [rr.strip() for rr in run_list if rr.strip() != '']

    # print(run_list)
    f = open(check_log, 'r')
    line = f.readlines()
    print("check log :{}\n".format(line))
    f.close()
    for s in line:
        print("{}\n".format(s))
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, line))
    finish_run = list(map(lambda x: x, finish))
    need_run = list(filter(lambda x: x not in finish_run, run_list))
    print("finish: {}\nfinish_run: {}\nneed_run: {}".format(finish, finish_run, need_run))
    print(
        "finish length: {}\nfinish_run length: {}\nneed_run length: {}".format(len(finish), len(finish_run),
                                                                               len(need_run)))
    sra_num_ = len(need_run) + len(finish_run)
    finish_num = 0
    print("len(finish_run)+len(need_run) = {}".format(sra_num_))
    num = len(finish_run)
    progress_list = []
    prog_num = 0

    finish_num = len(finish_run)
    finish_num_ = len(finish_run)
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
            print("########## {}/{} ###########".format(finish_num, sra_num_))
            # utils_.run_cmd("sbatch -A MST109178 -J Job_test -p ngs48G -c 14 --mem=46g -o ./out/{}_array_out.log -e ./out/{}_array_out.log "
            #               "--mail-user=sj985517@gmail.com --mail-type=BEGIN,END --wrap='/home/linsslab01/miniconda/bin/python3 one_Analysis.py'--array=1-4")
            with open(needList, "a+") as f:
                f.write("{}\n".format(k))

            prog_num += 1
            finish_num += 1
            time.sleep(1)
        job_out_o = os.path.join(job_out, "test_array_%a.log")
        job_out_err = os.path.join(job_out, "array_err_%a.log")
        with open(job_file, "w+") as f:
            f.write("#!/usr/bin/sh\n")
            f.write("#SBATCH -A MST109178\n")
            f.write("#SBATCH -J {}_job\n".format(pdat))
            f.write("#SBATCH -p ngs448core\n")
            f.write("#SBATCH -c 56\n")
            f.write("#SBATCH -n 8\n")
            f.write("#SBATCH --array=1-{}\n".format(len((need_run))))
            f.write("#SBATCH -o {}\n".format(job_out_o))
            f.write("#SBATCH -e {}\n".format(job_out_err))
            f.write("#SBATCH --mail-user=sj985517@gmail.com\n")
            f.write("#SBATCH --mail-type=BEGIN,END\n")
            f.write("echo $SLURM_ARRAY_TASK_ID\n")
            f.write("/home/linsslab01/miniconda3/bin/python3 one_Analysis.py {} {} {}\n".format(pdat,
                                                                                                sra_num_,
                                                                                                "$SLURM_ARRAY_TASK_ID"))
        print("sbatch {}".format(job_file))
        #utils_.run_cmd("sbatch {}".format(job_file))
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
            f.write("{} :\n{}\n".format(pdat, errMsg))
        sys.exit(errMsg)
    with open("./Automate_check.log", "a+") as f:
        f.write("{}\n".format(pdat))

def main(yy,mon,d):
    pool=multiprocessing.Pool(processes=50)
    global fin_num
    fin_num=0
    pattern = "salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]"

    date = datetime.date(yy, mon, d).strftime("%Y/%m/%d")
    # temp="{}/{}/{}".format(str(2020),str(mon+1),str(d))
    ######
    pdat = date.replace("/", "")
    new_outdir = os.path.join(outdir, pdat)
    utils_.mkdir_join(new_outdir)
    print("output: {}\n".format(new_outdir))
    sra_dir = os.path.join(new_outdir, "sra")  # .sra file
    utils_.mkdir_join(sra_dir)

    pattern, count = utils_.count_egquery(pattern, date, date)
    print("pattern: {}\ncount: {}\n".format(pattern, count))

    idlist = utils_.IdList_esearch(pattern, 'sra', count)

    print(idlist)

    runinfo = utils_.Get_RunInfo(idlist)
    run_list = list(runinfo['Run'])  # get SRAfile nameList stored in run_list
    print("runinfo: {}\n run_list: {}\n".format(runinfo, run_list))

    sraList = os.path.join(new_outdir, "sraList.txt")
    myfile2 = Path(sraList)
    myfile2.touch(exist_ok=True)
    f = open(sraList, 'r')
    line = f.readlines()
    print("check log :{}\n".format(line))
    f.close()


    for s in line:
        print("{}\n".format(s))
    finish = list(filter(lambda x: len(x.split(",")) >= 4, line))
    finish_run = list(map(lambda x: x.split(" ")[1], finish))
    need_run = list(filter(lambda x: x not in finish_run, run_list))
    print("finish: {}\nfinish_run: {}\nneed_run".format(finish, finish_run, need_run))
    print(
        "finish length: {}\nfinish_run length: {}\nneed_run length: {}".format(len(finish), len(finish_run),
                                                                               len(need_run)))

    sra_num_=len(need_run)+len(finish_run)
    run_num=1
    ######
    check_log = os.path.join(new_outdir, "Analysischeck.log")
    global fin_number
    fin_number=len(finish)
    for aa in need_run:
        try:
            if fin_num>0 and fin_num%200==0: ##if download 200 to sbatch ,and next sbatch has n-200
                print("fin_num%200==0\n")
                print("fin_num={}\n".format(fin_num))
                sbatch_job(outdir,pdat)
            else:
                print("#########################\nhello {}\n".format(aa))
                pool_list.append(pool.apply_async(sra_stat, (aa, new_outdir, sra_dir,sra_num_,len(need_run),date)))
                # pool.apply_async(test, (k,new_outdir,))
                # sra_stat(aa, new_outdir, sra_dir)
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
            with open("./SRA_run_error.txt", "a+") as f:
                f.write("{} :\n{}\n".format(date, errMsg))

    # if fin_run<200
    print("need_run end. fin_run<200\n")
    print("fin_num={}\n".format(fin_num))
    sbatch_job(outdir, pdat)

    pool.close()
    print("pool.close()\n")
    pool.join()
    print("pool.join()\n")

if __name__ == '__main__':
    start = time.time()
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    ## read SRAsetting.txt
    utils_.progress_bar("read SRAsetting.txt")
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




    for yy in range(sd_Y, ed_Y + 1):
        Month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
        ########
        if (yy % 4) == 0:
            if (yy % 100) == 0:
                if (yy % 400) == 0:
                    Month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
            else:
                Month = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]

        ########
        if sd_Y != ed_Y:
            if sd_Y == yy:
                sM = sd_M
                eM = 12
                sD = sd_D
                print("a")
            elif yy == ed_Y:
                sM = 1
                eM = ed_M
                sD = 1
                print("b")
        else:
            sM = sd_M
            eM = ed_M
            sD = 1

        ########
        for mon in range(sM, eM + 1):
            ###########
            if yy != ed_Y:
                eD = Month[mon - 1]
                sD = 1
            else:
                if mon != eM:
                    eD = Month[mon - 1]
                    sD = 1
                else:
                    eD = ed_D
                    sD = sd_D
            ########
            pool_list = []
            for d in range(sD, eD + 1):
                ######
                date = datetime.date(yy, mon, d).strftime("%Y/%m/%d")
                pdat = date.replace("/", "")
                new_outdir = os.path.join(outdir, pdat)
                utils_.mkdir_join(new_outdir)
                check_log = os.path.join(new_outdir, "Analysischeck.log")
                myfile2 = Path(check_log)
                myfile2.touch(exist_ok=True)
                with open(check_log, "a+") as f:
                    f.write(str(datetime.datetime.now()).split(".")[0])
                    f.write("\n")
                ##########
                main(yy,mon,d)



    print("Program Done\n")
    print(str(datetime.datetime.now()), ' Done,current total cost', time.time() - start, 'secs\n')
