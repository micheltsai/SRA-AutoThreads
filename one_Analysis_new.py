from __future__ import print_function

import csv
import datetime
import errno
import glob
import logging
import multiprocessing
import os
import random
import re
import shlex
import shutil
import signal
import subprocess
import sys
import time
import traceback
from pathlib import Path
import pandas as pd

import utils_


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
shovill_RAM=str(setting_df['shovill_RAM'][0])
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

def run_cmd(cmd):
    cmd=shlex.split(cmd)
    print(cmd)
    p = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    #print (cmd)
    print("--------------------------------------\n{}\noutput:\n".format(cmd))
    while p.poll() is None:
        #progress_bar("sub excuting")
        line = p.stdout.readline()
        line = line.strip()
        #if line:
        #    line_=line.decode().split("\n")
        #    for s in line_:
        #        print (str("{}\n".format(s)))
        #    sys.stdout.flush()
        #    sys.stderr.flush()
    if p.returncode ==0:
        print ("Subprogram sucess")
    else:
        print ("Subprogram failed")

    print ("-------------------------\n")
    return p

def progress_bar(Category):
    for i in range(1, 101):
        print("\r{}: ".format(Category),end="")
        print("[{}] {}%".format("*" * (i // 2), i), end="")
        sys.stdout.flush()
        time.sleep(0.02)
    print ("\n")

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
        #utils_.prefetch_sra(x, sra_dir)
        #print("no Download {}\n.".format(x))
        #with open(Downloadcheck_log, "a+") as f:
        #    f.write("{}\n".format(x))
        sys.exit("no Download {}\n{} is no exist.\n".format(sra_file))
    dltime=time.time() - one_
    print('Done,total cost',dltime, 'secs')
    print("###########################################################")



def Assembled(x,_outdir,sra_dir,ass_dir,assemble_dir,fastq_dir,thread,gsize,start):
    ###
    print("Assemble du-sh\n")
    utils_.run_cmd("du ./SRAtest -sh")
    time.sleep(1)
    ###
    final_dir = os.path.join(ass_dir, "{}_contig.fa".format(x))
    check_log = os.path.join(_outdir, "Asembledcheck.log")
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
    else:
        # one_run_prefetch = time.time()
        # utils_.prefetch_sra(x,sra_dir)
        one_run_ass = time.time()
        sra_file_dir = os.path.join(sra_dir, x)
        print("sra_file_dir:",sra_file_dir)
        sra_file = os.path.join(sra_dir, "{}/{}.sra".format(x, x))
        print("sra_file:",sra_file)
        if os.path.isfile(sra_file):
            print("have {}.sra file.\n".format(x))
        else:
            utils_.prefetch_sra(x, sra_dir)
            print("not found {}.sra, Download now\n".format(x))

        utils_.run_for_114v3(x, sra_dir, fastq_dir, assemble_dir, _outdir, thread, gsize, start, check_log, shovill_RAM)
        ### unnecessary ERR file
        #ERR_path = os.path.join(os.path.abspath(os.getcwd()), x)
        #print("suERR_path: ", ERR_path, "\n")
        # print ("shutil.rmtree({})\n".format(current_path))
        #utils_.run_cmd2("rm -rf {}".format(current_path))
        #print("remove {}\n".format(current_path))
        #rmsra_cmd="rm -rf {}".format(sra_file_dir)
        #print(rmsra_cmd)
       # print("remove {}.sra".format(x))
        #run_cmd(rmsra_cmd)


def QualityCheck(sra_id,_outdir,ori_outdir,genome_Path,thread,gsize,start):
    ###
    #print("Qualitycheck_start du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(1)
    ###
    print("#####################  QualityCheck  #####################\n")

    Assem_path = os.path.join(_outdir, "Assembled/")
    BUSCOresult = os.path.join(_outdir, "BUSCOresult.txt")
    check = os.path.join(_outdir, "QCcheck.log")
    outdir = os.path.join(_outdir, "QualityCheck")
    utils_.mkdir_join(outdir)

    # outdir = utils_.mkdir_join(outdir, str(current_time))
    print("outdir: \n", _outdir)
    print("check: \n", check)
    print("BUSCOresult= {}".format(BUSCOresult))

    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    # genome_Path = utils_.getGenomeListPath(args.genome, outdir)



    gID = genome_Path.replace(Assem_path, "")
    # gID=os.path.basename(genome_Path)
    gID = gID.split(".")[0]
    print("gID: {}\n".format(gID))
    refDIR=os.path.join(_outdir,"{}_ref".format(gID))
    utils_.mkdir_join(refDIR)
    refPath = utils_.getRefListPath(ref_dir, refDIR)
    # refPath=args.ref

    # outdir_ani = os.path.join(outdir, 'fastani')
    outdir_ani = os.path.join(outdir, 'fastani')
    utils_.mkdir_join(outdir_ani)
    print("outdir_ani: {}\n".format(outdir_ani))
    # outdir_ani=os.path.join(outdir, 'fastani')

    outfile_ = '{}_ani.txt'.format(gID)

    outfile = os.path.join(outdir_ani, outfile_)  # stroed fastANI output in out.txt

    info_txt = os.path.join(outdir_ani, '{}_info.txt'.format(gID))  # stroed fastANI output in out.txt
    db = buscoDB
    mode = buscoMode


    ##fastANI-------
    fastANI_time = time.time()
    current_path = os.path.abspath(os.getcwd())
    print("current_path: ", current_path, "\n")
    replace_path = outdir.replace(current_path, ".")
    fastani_outdir = os.path.join(replace_path, '{}_ani.txt'.format(gID))

    print("-------------------------------fastANI start.-------------------------------")
    print("reseq: {}\n qen: {}\n outdir: {}\nout_txt: {}\n{}\n".format(refPath, genome_Path, outdir, outfile,
                                                                       os.path.join(outdir_ani, outfile_)))
    #utils_.progress_bar("fastANI excuting")
    #fastani_ = "/data/usrhome/LabSSLin/user30/Desktop/FastANI/fastANI -t {} --rl {} -q {} -o {}".format(thread,refPath, genome_Path, outfile)
    fastani_ = "/home/linsslab01/FastANI/fastANI -t {} --rl {} -q {} -o {}".format(thread, refPath,genome_Path,outfile)
    print(fastani_ + "\n")
    if os.path.isfile(outfile):
        print(outfile, " is exist.\n")
        print("fastANI was done.\n")
    else:
        utils_.run_cmd(fastani_)
        print("fastANI done.\n")


    # ANI>=95------
    print(
        "-------------------------------fastANI end.-------------------------------\ncompare and calculate ANI\nget ANIoutPath\n")

    rm_ref_cmd="rm -rf {}".format(refDIR)
    print(rm_ref_cmd)
    utils_.run_cmd(rm_ref_cmd)

    print("open fastANI output txt\n")
    # open fastANI output
    f = open(outfile, 'r')
    AverageANI = 0.0
    num = 0  # quantity of ANI>=95
    not_num = 0  # quantity of ANI<95
    ANI_total = 0.0  # total of all ANI value
    ANI_ = f.readlines()  # read file stored in ANI_
    for x in ANI_:
        tmp = float(x.split("\t")[2])  # temporarily ANI value
        if (tmp >= 95.0):
            num += 1
            ANI_total += tmp
        else:
            not_num += 1
    # num =0
    QC_error=os.path.join(ori_outdir,"nofillQC.txt")
    print("QC_errfile:", QC_error)
    QC_file = Path(QC_error)
    QC_file.touch(exist_ok=True)
    if num==0:
        with open(QC_error, "a+") as f:
            f.write("{}: all ANI value < 95\n".format(sra_id))
        sys.exit("all ANI value <95\n")
    AverageANI = ANI_total / num  # if ANI>=95 ,calulate average value of ANI
    print("Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num + not_num, num,
                                                                                     ANI_[0].split("\t")[2]))
    targetPath = ANI_[0].split("\t")[1]

    # save data info
    with open(info_txt, "w+") as f2:
        f2.write(
            "Average ANI: {}\ntotal number: {}\n>= quantity: {}\nmax ANI: {}\n".format(AverageANI, num + not_num, num,
                                                                                       ANI_[0].split("\t")[2]))
        f2.write("targetPath: {}\n".format(targetPath))

    f.close()

    if num == 0:
        print("All ANI value doesn't exceed 95, next genome run\n")
        with open(check, "a+") as f:
            f.write("{} is ANI<95.\n".format(gID))
        return 0


    # BUSCO------
    ###
    #print("busco du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(1)
    ###
    print("-------------------------------ANI>=95 continue, BUSCO start-------------------------------\n")
    # use conda enterring the busco VM(vm name is "busco")
    busco_time = time.time()
    outdir_bus = os.path.join(outdir, 'busco_db')
    busco_db="/work/linsslab01/busco_db"
    #busco_db = utils_.mkdir_join(outdir_bus)

    # -f overwrite

    # genome is "one excuting"
    # busco -i /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/RefSeq/GCF_000335875.2.fa -o cofig --out_path /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/QualityCheck -l enterobacterales_odb10 -m geno --download_path /data1/usrhome/LabSSLin/linss01/Desktop/SRA-AutoAnalysis/SRA-AutomatedAnalysis/QualityCheck/QualityCheck/busco_db -f
    #cmd_bus = 'bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco && busco -c {} -i {} -o {} --out_path {} -l {} -m {} --download_path {} -f"'.format(
    cmd_bus = 'bash -c "source /home/linsslab01/miniconda3/etc/profile.d/conda.sh && conda activate busco && busco -c {} -i {} -o {} --out_path {} -l {} -m {} --download_path {} -f --offline"'.format(
        thread,targetPath, gID, outdir, db, mode, busco_db)
    print(cmd_bus, "\n")
    utils_.run_cmd(cmd_bus)
    ###
    #print("busco du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(1)
    ###
    # get BUSCO complete>=95 & duplicate>=3 ,or exit
    buscopath = os.path.join(outdir, "{}".format(gID))
    buscopath = os.path.join(buscopath, "run_{}".format(db))
    buscopath = glob.glob(buscopath + "/*.txt")
    print(buscopath)
    buscopath = os.path.abspath(buscopath[0])
    print(buscopath)
    with open(buscopath, "r") as f:
        b = f.readlines()
    #print(b, "\n")
    b = b[8].strip("\t")
    #print(b)
    #print([float(s) for s in re.findall(r'-?\d+\.?\d*', b)])

    bC, bS, bD, bF, bM, bn = [float(s) for s in re.findall(r'-?\d+\.?\d*', b)]
    #print("c:{}%, d:{}%\n".format(bC, bD))

    with open(BUSCOresult, "a+") as f:
        print("stored BUSCO result\n")
        f.write("{}: {}\n".format(gID, [float(s) for s in re.findall(r'-?\d+\.?\d*', b)]))

    if bC < 95.0 and bD > 3.0:
        with open(check, "a+") as f:
            f.write("{} is C<95 or D>3.\n".format(gID))
        return 0

    # continue
    targettxt = os.path.join(_outdir, "target.txt")
    print("target.txt path: {}".format(targettxt))
    # continue need target path
    with open(targettxt, "a+") as f:
        f.write("{}:{}\n".format(genome_Path, targetPath))
    # check
    with open(check, "a+") as f:
        f.write("{} is ok.\n".format(gID))
        print("commit on check \n")

    print('Done,total cost', time.time() - start, 'secs\n')

    return targetPath

def Analysis(sra_id,input,target_ref,anoutdir,_outdir,thread,gsize,start):
    ###
    print("Analysis du-sh\n")
    utils_.run_cmd("du ./SRAtest -sh")
    time.sleep(1)
    ###
    print("#####################  Analysis  #####################\n")
    mlst_organism = mlstS
    amr_organism = amrS
    utils_.mkdir_join(anoutdir)

    # get input id
    inlist = input.split("/")
    inId = inlist[len(inlist) - 1]
    print("input Id: {}\n".format(inId))
    inId = inId.split(".")[0]
    print("input Id: {}\n".format(inId))


    # workdir
    current_path = os.path.abspath(os.getcwd())
    current_path2 = current_path.replace("/SRA-AutoThreads", "")
    print("current_path: ", current_path, "\n")
    relative_input_ = input.replace(current_path2, ".")
    relative_input = input.replace(current_path, ".")
    print("relative input: {}\n".format(relative_input))
    print("relative input_: {}\n".format(relative_input_))
    origin_outdir = _outdir
    allinfopath = os.path.join(origin_outdir, "info.txt")
    check = os.path.join(origin_outdir, "Anacheck.log")

    # add outpath "analysis"
    utils_.mkdir_join(_outdir)
    anoutdir_ = os.path.join(_outdir, "analysis")
    utils_.mkdir_join(anoutdir_)
    print("analysis outdir: {}\n".format(anoutdir_))

    # set {genomoe}_log_output
    logpath = os.path.join(anoutdir_, "analysis_log.txt")

    # get relative output dir path
    outdir_list = anoutdir_.split("/")
    relative_path_o2 = anoutdir_.replace(current_path, ".")
    relative_path2 = anoutdir_.replace(current_path2, ".")
    print("relative2: {}\n".format(relative_path2))
    print("relative_path: {}".format(relative_path2))
    #relative_path_o2 = os.path.join(relative_path_o2, inId)
    #relative_path2 = os.path.join(relative_path2, inId)
    utils_.mkdir_join(relative_path_o2)
    utils_.mkdir_join(relative_path2)
    print("relative_path: {}".format(relative_path2))

    # load log.txt read running statedat
    step = 0
    filename = Path(logpath)
    filename.touch(exist_ok=True)

    with open(logpath, "r") as f:
        line = f.readlines()
        print(line)
        for x in line:
            ana = x.split(" ")[0]
            if ana == "mlst":
                step = 1
            elif ana == "plasmidfinder":
                step = 2
            elif ana == "amr":
                step = 3
            elif ana == "sistr":
                step = 4
            print("ana: {}, step: {}\n".format(ana, step))

    # run MLST
    if step < 1:
        step1_time = time.time()
        print("STEP{}\n".format(step + 1))

        print("********** Now MLST analysis running. **********\n")
        #MLST_DB = "/home/linsslab01/mlst_db"
        mlst_outdir = os.path.join(anoutdir_, "mlst")
        utils_.mkdir_join(mlst_outdir)
        mlst_tmp=os.path.join(mlst_outdir,"tmp")

        mlst_datajson = os.path.join(mlst_outdir, "data.json")
        #mlst_cmd = "docker run --rm -it \-v {}:/databases \-v {}:/workdir \mlst -i {} -o {} -s {}".format(MLST_DB,current_path,relative_input,mlst_outdir,mlst_organism)

        with open(mlst_datajson,"w+") as f:
            f.close()

        #mlst_cmd = "singularity exec --containall --bind /work/linsslab01/:/home/linsslab01/ /work/linsslab01/mlst.sif python3 /home/linsslab01/mlst/mlst.py -i {} -o {} -s {}".format(relative_input_.replace("work","home"),mlst_outdir.replace("work", "home"),mlst_organism)
        mlst_cmd="/home/linsslab01/miniconda3/bin/python3 /work/linsslab01/mlst/mlst.py -i {} -o {} -s {} -t {}".format(relative_input,mlst_outdir,mlst_organism,mlst_outdir)
        print(mlst_cmd, "\n")
        try:
            time.sleep(random.randint(0,30))
            mlst, err = utils_.run_cmd3(mlst_cmd)

            with open(logpath, "a+") as f:
                if mlst.returncode != 0:
                    # print(mlst.stdout.readline())
                    # print(err)
                    f.write(str(err))
                    f.write("\n")
                    sys.exit(err)
                else:
                    f.write("mlst is ok\n")
            step += 1
            ###
            #print("mlst du-sh\n")
            #utils_.run_cmd("du ./SRAtest -sh")
            #time.sleep(1)
            ###

            rmTMP_cmd="rm -rf {}\n".format(mlst_tmp)
            print(rmTMP_cmd)
            run_cmd(rmTMP_cmd)
            time.sleep(1)
            ###
            #print("rm_mlstTMP du-sh\n")
            #utils_.run_cmd("du ./SRAtest -sh")
            ###
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
            print("mlst err\n")
            with open("./err_need_again.txt","a+") as f1:
                f1.write(sra_id)
                f1.close()
            utils_.run_cmd(mlst_cmd)
            pass




    else:
        print("**********       mlst was running.      **********\n next step\n")
    ###
    #print("mlst du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(2)
    ###
    # run plasmidfinder
    print("**********       mlst end.      **********\n next step\n")
    if step < 2:
        step2_time = time.time()
        print("STEP{}\n".format(step + 1))
        ###
        #print("plas du-sh\n")
        #utils_.run_cmd("du ./SRAtest -sh")
        #time.sleep(2)
        ###
        print("********** Now plasmidfinder analysis running. **********\n")
        #PLASMID_DB = "/home/linsslab01/plasmidfinder_db"
        utils_.mkdir_join(anoutdir_)
        plas_outdir = os.path.join(anoutdir_, "plasmidfinder")
        utils_.mkdir_join(plas_outdir)
        plas_outdir=plas_outdir.replace("work","home")
        #plas_cmd = "docker run --rm -it \-v {}:/databases \-v {}:/workdir \plasmidfinder -i {} -o {}".format(PLASMID_DB, current_path, relative_input, plas_outdir)
        plas_cmd = "singularity exec --containall --bind /work/linsslab01/:/home/linsslab01/ /work/linsslab01/plasmidfinder.sif python3 /home/linsslab01/plasmidfinder/plasmidfinder.py -i {} -o {}".format(relative_input_.replace("work","home"),plas_outdir)

        print(plas_cmd, "\n")
        plas = utils_.run_cmd(plas_cmd)

        with open(logpath, "a+") as f:
            if plas.returncode != 0:
                # print(mlst.stdout.readline())

                sys.exit()
            else:
                f.write("plasmidfinder is ok\n")
        step += 1
        # time
    else:
        print("********** plasmidfinder was running. **********\n next step\n")
    time.sleep(1)
    # run amrfinder
    print("**********       plasmidfinder end.      **********\n next step\n")
    ###
    #print("plas du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(2)
    ###

    if step < 3:
        step3_time = time.time()
        print("STEP{}\n".format(step + 1))
        ###
        #print("amr du-sh\n")
        #utils_.run_cmd("du ./SRAtest -sh")
        #time.sleep(2)
        ###
        print("********** Now amrfinder analysis running. **********\n")
        amr_outdir = os.path.join(relative_path_o2, "amrfinder")
        utils_.mkdir_join(amr_outdir)
        amr_outdir = os.path.join(amr_outdir, "amrout.tsv")
        amr_cmd = "amrfinder -n {} -o {} -O {}".format(input, amr_outdir, amr_organism)
        print(amr_cmd, "\n")
        amr = run_cmd(amr_cmd)
        with open(logpath, "a+") as f:
            if amr.returncode != 0:
                sys.exit()
            else:
                f.write("amr is ok\n")
        step += 1
        # time
    else:
        print("**********   amrfinder was running.   **********\n next step\n")
    ###
    #print("amr du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(2)
    ###
    print("**********       amrfinder end.      **********\n next step\n")
    # run sistr
    if step < 4:
        step4_time = time.time()
        print("STEP{}\n".format(step + 1))
        ###
        #print("sistr du-sh\n")
        #utils_.run_cmd("du ./SRAtest -sh")
        #time.sleep(2)
        ###
        print("********** Now sistr analysis running. **********")
        sistr_outdir = os.path.join(relative_path_o2, "sistr")
        utils_.mkdir_join(sistr_outdir)
        sistr_outdir = os.path.join(sistr_outdir, "sistr_out")
        input_list = input.split("/")
        input_name = input_list[len(input_list) - 1]
        print("name: ", input_name, "\n")
        sistr_cmd = "sistr --threads {} -i {} {} -f csv -o {} -m".format(thread,input, input_name, sistr_outdir)
        print(sistr_cmd, "\n")

        sistr = run_cmd(sistr_cmd)
        with open(logpath, "a+") as f:
            if sistr.returncode != 0:
                sys.exit()
            else:
                f.write("sistr is ok\n")
        step += 1
        # time
    else:
        print("********** sistr was running. **********\n next step\n")
    ###
    #print("sistr du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(2)
    ###

    print("**********       sistr end.      **********\n next step\n")
    ########################
    ########################
    # save data in analysis_final.csv and update DB

    # read mlst 'Sequence Type'
    mlst_file = os.path.join(relative_path_o2, "mlst/results.txt")
    print("mlst_file:",mlst_file)
    try:
        with open(mlst_file, "r") as f:
            data = f.readlines()
            print(data[6])
            # Sequence Type: 11
            sequenceType = data[6].split(" ")
            sequenceType = sequenceType[len(sequenceType) - 1].strip("\n")

            print(sequenceType)
    except Exception as e:
        time.sleep(3)
        print("not found mlst data.json\n")
        mlst_outdir = os.path.join(anoutdir_, "mlst")
        run_cmd("rm -rf {}".format(mlst_outdir))

        mlst_cmd = "singularity exec --containall --bind /work/linsslab01/:/home/linsslab01/ /work/linsslab01/mlst.sif python3 /home/linsslab01/mlst/mlst.py -i {} -o {} -s {}".format(
            relative_input_.replace("work","home"), mlst_outdir.replace("work", "home"), mlst_organism)
        print(mlst_cmd)
        mlst, err = utils_.run_cmd3(mlst_cmd)
        mlst_file = os.path.join(relative_path_o2, "mlst/results.txt")
        with open(mlst_file, "r") as f:
            data = f.readlines()
            print(data[6])
            # Sequence Type: 11
            sequenceType = data[6].split(" ")
            sequenceType = sequenceType[len(sequenceType) - 1].strip("\n")
        pass
    # read plasmidfinder 'gene'
    try:
        plas_file = os.path.join(relative_path_o2, "plasmidfinder/results_tab.tsv")
        pladf = pd.read_table(plas_file, sep='\t')
        # print(df)
        pladf = pd.DataFrame(pladf)
    except Exception as e:
        time.sleep(3)
        print("not found mlst data.json\n")
        plas = run_cmd(plas_cmd)
        plas_file = os.path.join(relative_path_o2, "plasmidfinder/results_tab.tsv")
        pladf = pd.read_table(plas_file, sep='\t')
        # print(df)
        pladf = pd.DataFrame(pladf)
        pass
    print(pladf)
    print(pladf.columns)
    print(pladf.Plasmid)
    plist = list(pladf.Plasmid)
    plas_format = ""

    for x in range(0, len(plist)):
        plas_format += plist[x]
        if x <= len(plist) - 1:
            plas_format += ","
    print(plas_format)

    # read amrfinder 'Gene symbol', 'subtype'
    ##add amr "Point"
    amr_file = os.path.join(relative_path_o2, "amrfinder/amrout.tsv")
    amrdf = pd.read_table(amr_file, sep='\t')
    # print(df)
    amrdf = pd.DataFrame(amrdf)
    print(amrdf)
    print(amrdf.columns)
    # replace ' ' as '_'
    amrdf.columns = ['Protein_identifier', 'Contig_id', 'Start', 'Stop', 'Strand',
                     'Gene_symbol', 'Sequence_name', 'Scope', 'Element_type',
                     'Element_subtype', 'Class', 'Subclass', 'Method', 'Target_length',
                     'Reference_sequence_length', 'Coverage_of_reference_sequence',
                     'Identity_to_reference_sequence', 'Alignment_length',
                     'Accession_of_closest sequence', 'Name_of_closest sequence', 'HMM_id',
                     'HMM_description']

    print(amrdf.columns)
    # if list is [], show NaN
    print(amrdf.Gene_symbol)
    print(amrdf.Element_subtype)

    amr_format = ""
    point_format = ""
    sym_list = list(amrdf.Gene_symbol)
    sub_list = list(amrdf.Element_subtype)
    for i in range(0, len(sym_list)):
        if sub_list[i] == "AMR":
            if len(amr_format) > 0 and i!=len(sym_list):
                amr_format += ","
            amr_format += sym_list[i]
        elif sub_list[i] == "POINT":
            if len(point_format) > 0 and i!=len(sym_list):
                point_format += ","
            point_format += sym_list[i]

    print(amr_format)
    print(point_format)

    sistr_file = os.path.join(relative_path_o2, "sistr")
    utils_.mkdir_join(sistr_file)
    sistr_file = os.path.join(sistr_file, "sistr_out.csv")

    sistrdf = pd.read_csv(sistr_file)
    # print(df)
    sistrdf = pd.DataFrame(sistrdf)
    print(sistrdf)
    print(sistrdf.columns)
    print(sistrdf.serovar)

    in_abspath = input.replace(".", current_path)

    dict = {'Accession': sra_id,
            'MLST': sequenceType,
            'AMR': amr_format,
            'Point': point_format,
            'Serotype': sistrdf.serovar,
            'Inc Type': plas_format
            }

    finaldf = pd.DataFrame(dict)
    print(finaldf)
    finalfile = os.path.join(outdir, "analysis_final.csv")
    ###
    #print("finalfile du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(1)
    ###
    print("finaldf.to_csv(finalfile, mode='a+', header=False)\n")
    finaldf.to_csv(finalfile, mode='a+', header=False)
    ###
    #print("finalfile du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(1)
    ###
    # after run all state, save ID in "Anackeck.log" and remove ./analysis
    with open(check, "a+") as f:
        f.write("Run {} is ok.\n".format(inId))
    return 0

def SRA_Analysis(sra_id,sra_dir,ass_dir,fastq_dir,assemble_dir,_outdir,thread,gsize,start,sra_num_,outdir):
    SRA_start=time.time()
    try:

        # if sra_layout==2 continue
        #Download(sra_id,_outdir,sra_dir)
        Assembled(sra_id,_outdir,sra_dir,ass_dir,assemble_dir,fastq_dir,thread,gsize,start)
        with open("./checkAssembled.txt","a+") as f:
            f.write("Run {} is ok.\n".format(sra_id))
        #####

        genome = os.path.join(ass_dir, "{}_contig.fa".format(sra_id))
        #sys.exit()
        targetPath=QualityCheck(sra_id,_outdir,outdir,genome,thread,gsize,start)
        time.sleep(1)
        with open("./checkQC.txt","a+") as f:
            f.write("Run {} is ok.\n".format(sra_id))
        print("targetPAth = {}\n######\n".format(targetPath.encode("utf-8").decode()))
        target_ = targetPath.replace(current_path, ".")
        print("target_= {}\n".format(target_))
        time.sleep(1)

        Analysis(sra_id,genome,target_,_outdir,_outdir,thread,gsize,start)
        with open("./checkAnalysis.txt","a+") as f:
            f.write("Run {} is ok.\n".format(sra_id))
        #global finish_num_
        #finish_num_ += 1
        check_log = os.path.join(_outdir, "Analysischeck.log")
        #with open(check_log, "a+") as f:
        #    f.write("Run {} is ok.\n".format(sra_id))

        ########
        utils_.run_cmd2("rm -rf {}".format(fastq_dir))
        print(fastq_dir)
        utils_.run_cmd2("rm -rf {}".format(sra_dir))
        print(sra_dir)
        utils_.run_cmd2("rm -rf {}".format(assemble_dir))
        print(assemble_dir)
        print("remove fastq dir, sra dir, assembled_result dir\n")

        print("Run {} is Done\n".format(sra_id))
        mycallback_write_Finish(sra_id)
        #######
        print(str(datetime.datetime.now()),' Done,current total cost', time.time() - start, 'secs\n')
    except Exception as e:
        error_class = e.__class__.__name__  # 取得錯誤類型
        detail = e.args[0]  # 取得詳細內容
        cl, exc, tb = sys.exc_info()  # 取得Call Stack
        lastCallStack = traceback.extract_tb(tb)[-1]  # 取得Call Stack的最後一筆資料
        fileName = lastCallStack[0]  # 取得發生的檔案名稱
        lineNum = lastCallStack[1]  # 取得發生的行號
        funcName = lastCallStack[2]  # 取得發生的函數名稱
        errMsg = "File \"{}\", line {}, in {}: [{}] {}\n".format(fileName, lineNum, funcName, error_class, detail)
        print(errMsg)
        ###
        #process = psutil.Process(os.getpid())
        #print(str(datetime.datetime.now()), process.memory_info().rss)
        #utils_.run_cmd("free -h")
        ####
        with open("./SRA_run_error.txt", "a+") as f:
            f.write("{} :\n{}\n".format(sra_id, errMsg))

        sys.exit(e)
    #sys.exit("subpreocess End\n")
    return sra_id

def mycallback_write_Finish(sra_id):
    print("mycallback_write\n")
    check_log = os.path.join("./SRAtest", "Analysischeck.log")
    with open(check_log, "a+") as f:
        f.write("Run {} is ok.\n".format(sra_id))
        check_lines = f.readlines()
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, check_lines))
    finish_num_s = list(map(lambda x: x.split(" ")[1], finish))
    print("finish num={}\n".format(len(finish_num_s)))
    print("{} / {}\n".format(len(finish_num_s), sra_num_))

if __name__ == '__main__':
    start = time.time()
    argvs=sys.argv

    print(argvs)

    pdat=argvs[1]
    #need_run=argvs[2]
    sra_num_ = argvs[2]
    sra_index=int(argvs[3])-1
    txt_index=argvs[4]
    #sra_id_test=argvs[5]


    needList = os.path.join(outdir, "need_run_{}.txt".format(txt_index))

    with open(needList,"r") as f:
        needlines=f.readlines()
    need_run=needlines
    need_run = [rr.strip() for rr in need_run if rr.strip()!='']
    print(need_run)
    sra_id=need_run[sra_index].strip("\n")
    new_outdir = os.path.join(outdir, "output")
    utils_.mkdir_join(new_outdir)
    new_outdir = os.path.join(new_outdir, sra_id)
    utils_.mkdir_join(new_outdir)
    print("output: {}\n".format(new_outdir))

    check_log = os.path.join(new_outdir, "Analysischeck.log")
    sra_dir = os.path.join(new_outdir, "sra")  # .sra file
    utils_.mkdir_join(sra_dir)
    ass_dir = os.path.join(new_outdir, "Assembled")
    utils_.mkdir_join(ass_dir)
    fastq_dir = os.path.join(new_outdir, 'fastq')
    utils_.mkdir_join(fastq_dir)
    assemble_dir = os.path.join(new_outdir, "assembly_result")
    utils_.mkdir_join(assemble_dir)

    print("sra_dir:{}\nass_dir={}\nfastq_dir={}\nassemble_dir={}\n".format(sra_dir, ass_dir, fastq_dir, assemble_dir))
    #print(need_run)
    #print(need_run.type)
    print(sra_index)
    print(sra_id)

    SRA_Analysis(sra_id, sra_dir, ass_dir, fastq_dir, assemble_dir, new_outdir, thread, gsize, start, sra_num_, outdir)

    #shutil.rmtree(sra_dir)
   #utils_.mkdir_join(fastq_dir)
    #shutil.rmtree(fastq_dir)
    #shutil.rmtree(assemble_dir)

