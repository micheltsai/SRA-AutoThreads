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
        #utils_.mkdir_join(sra_dir)
        #utils_.prefetch_sra(x, sra_dir)
        #print("no Download {}\n.".format(x))
        #with open(Downloadcheck_log, "a+") as f:
        #    f.write("{}\n".format(x))
        sys.exit("no Download {}\n{} is no exist.\n".format(sra_file))
    dltime=time.time() - one_
    print('Done,total cost',dltime, 'secs')
    print("###########################################################")

def Assembled(x,_outdir,sra_dir,ass_dir,assemble_dir,fastq_dir,thread,gsize,start):
    final_dir = os.path.join(ass_dir, "{}_contig.fa".format(x))
    check_log = os.path.join(_outdir, "Asembledcheck.log")
    if os.path.isfile(final_dir):
        print("was ran assembly ,contig.fa is exist\n------------------------------\n\n")
    else:
        utils_.mkdir_join(ass_dir)
        utils_.mkdir_join(fastq_dir)
        utils_.mkdir_join(assemble_dir)
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

    anifile_new = os.path.join(outdir, outfile_)

    if os.path.isfile(anifile_new):
        print(anifile_new, " is exist.\n")
        print("fastANI was done.\n")
    else:
        utils_.run_cmd(fastani_)
        print("fastANI done.\n")


    ######move anifile to QCdir/

    mvani_cmd="cp {} {}".format(outfile,anifile_new)
    print(mvani_cmd)
    os.system(mvani_cmd)



    # ANI>=95------
    print(
        "-------------------------------fastANI end.-------------------------------\ncompare and calculate ANI\nget ANIoutPath\n")

    ######
    rm_ref_cmd="rm -rf {}".format(refDIR)
    print(rm_ref_cmd)
    utils_.run_cmd(rm_ref_cmd)
    ######
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
    check_dir = os.path.join(ori_outdir, "check")
    utils_.mkdir_join(check_dir)
    print("check_dir:{}\n".format(check_dir))

    QC_error=os.path.join(check_dir,"nofillQC.txt")
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
            f.write("{}:ANI<95.\n".format(gID))
        return 0

    ###################################
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
    # busco -i /data1/usrhome/LabSSLin/linss01/Desktop/web-AutoAnalysis/RefSeq/GCF_000335875.2.fa -o cofig --out_path /data1/usrhome/LabSSLin/linss01/Desktop/web-AutoAnalysis/web-AutomatedAnalysis/QualityCheck -l enterobacterales_odb10 -m geno --download_path /data1/usrhome/LabSSLin/linss01/Desktop/web-AutoAnalysis/web-AutomatedAnalysis/QualityCheck/QualityCheck/busco_db -f
    #cmd_bus = 'bash -c "source /data/usrhome/LabSSLin/user30/anaconda3/etc/profile.d/conda.sh && conda activate busco && busco -c {} -i {} -o {} --out_path {} -l {} -m {} --download_path {} -f"'.format(
    cmd_bus = 'bash -c "source /home/linsslab01/miniconda3/etc/profile.d/conda.sh && conda activate busco && busco -c {} -i {} -o {} --out_path {} -l {} -m {} --download_path {} -f --offline"'.format(
        thread,targetPath, gID, outdir, db, mode, busco_db)
    print(cmd_bus, "\n")

    buscofile_newpath = os.path.join(outdir, "busco_short_summary.txt")
    if os.path.isfile(buscofile_newpath):
        print(buscofile_newpath, " is exist.\n")
        print("Busco was done.\n")
    else:
        utils_.run_cmd(cmd_bus)
        print("Busco done.\n")
    ###
    #print("busco du-sh\n")
    #utils_.run_cmd("du ./SRAtest -sh")
    #time.sleep(1)
    ###
    # get BUSCO complete>=95 & duplicate>=3 ,or exit
    buscoDirpath = os.path.join(outdir, "{}".format(gID))
    buscopath = os.path.join(buscoDirpath, "run_{}".format(db))
    #buscopath = glob.glob(buscopath + "/*.txt")
    buscopath = os.path.join(buscopath , "short_summary.txt")


    mvbuscocmd="cp {} {}".format(buscopath,buscofile_newpath)
    print(mvbuscocmd)
    os.system(mvbuscocmd)

    print(buscofile_newpath)
    buscopath = os.path.abspath(buscofile_newpath)
    print(buscofile_newpath)
    with open(buscofile_newpath, "r") as f:
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
            f.write("{}:C<95 or D>3.\n".format(gID))
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

    #### remove
    rm_fastaniDir = "rm -rf {}".format(outdir_ani)
    print(rm_fastaniDir)
    os.system(rm_fastaniDir)

    print(rm_fastaniDir)
    rm_buscoDir="rm -rf {}".format(buscoDirpath)
    os.system(rm_buscoDir)

    print('Done,total cost', time.time() - start, 'secs\n')

    return targetPath

def Analysis(sra_id,input,target_ref,outdir,thread,gsize,start):
    print("#####################  Analysis  #####################\n")
    mlst_organism = mlstS
    amr_organism = amrS
    anoutdir=outdir
    utils_.mkdir_join(anoutdir)
    print("anoutdir:{}\n".format(anoutdir))
    print("_outdir:{}\n".format(outdir))

    # get input id
    inlist = input.split("/")
    inId = inlist[len(inlist) - 1]
    print("input Id: {}\n".format(inId))
    inId = inId.split(".")[0]
    print("input Id: {}\n".format(inId))


    # workdir
    current_path = os.path.abspath(os.getcwd())
    current_path2 = current_path.replace("/web-AutoThreads-test", "")
    print("current_path: ", current_path, "\n")
    print("current_path2: ", current_path2, "\n")
    relative_input_ = input.replace(current_path2, ".")
    relative_input = input.replace(current_path, ".")
    print("relative input: {}\n".format(relative_input))
    print("relative input_: {}\n".format(relative_input_))
    origin_outdir = outdir
    print("origin_outdir:".format(origin_outdir))
    check = os.path.join(origin_outdir, "Anacheck.log")
    print("check: {}\n".format(check))
    # add outpath "analysis"
    utils_.mkdir_join(outdir)
    anoutdir_ = os.path.join(outdir, "analysis")
    utils_.mkdir_join(anoutdir_)
    print("anoutdir_: {}\n".format(anoutdir_))



    # set {genomoe}_log_output
    logpath = os.path.join(anoutdir_, "analysis_log.txt")
    print("logpath: {}\n".format(logpath))
    # get relative output dir path
    relative_path_o2 = anoutdir_.replace(current_path, ".")
    #relative_path2 = anoutdir_.replace(current_path2, ".")
    #print("relative2: {}\n".format(relative_path2))
    #print("relative_path: {}".format(relative_path2))
    #relative_path_o2 = os.path.join(relative_path_o2, inId)
    #relative_path2 = os.path.join(relative_path2, inId)
    utils_.mkdir_join(relative_path_o2)
   # utils_.mkdir_join(relative_path2)
    #qprint("relative_path: {}".format(relative_path2))

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

    # after run all state, save ID in "Anackeck.log"
    with open(check, "a+") as f:
        f.write("Run {} is ok.\n".format(inId))
    return 0

def getBenga2(sra_id,outdir):
    Benga_start=time.time()
    outdir=os.path.join(outdir,"output")
    utils_.mkdir_join(outdir)
    input=os.path.join(outdir,"{}/Assembled/{}_contig.fa".format(sra_id,sra_id))
    output = os.path.join(outdir, sra_id)
    utils_.mkdir_join(output)
    output=os.path.join(output,"cgMLST")
    utils_.mkdir_join(output)
    output = os.path.join(output, "{}.tsv".format(sra_id))
    cmd="/home/linsslab01/miniconda3/bin/python3 ./20210210cgMLST/Benga-2/Benga-2/profiling.py -i {} -o {} " \
        "--scheme ./20210210cgMLST/Benga-2/scheme.faa --prodigaltf ./20210210cgMLST/Benga-2/prodigaltf.trn".format(input,output)
    print(cmd)
    try:
        os.system(cmd)
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
        check_dir = os.path.join(outdir, "check")
        utils_.mkdir_join(check_dir)
        print("check_dir:{}\n".format(check_dir))
        err_txt=os.path.join(check_dir,"Benga_err.txt")
        with open(err_txt, "a+") as f:
            f.write("{}:{}\n".format(sra_id,e))
        sys.exit(e)
    print('getBenga2 Done,{} total cost'.format(sra_id), time.time() - Benga_start, 'secs\n')

def SRA_Analysis(sra_id,_outdir,thread,gsize,start,sra_num_,outdir):
    SRA_start=time.time()
    sra_dir = os.path.join(new_outdir, "sra")  # .sra file
    ass_dir = os.path.join(new_outdir, "Assembled")
    fastq_dir = os.path.join(new_outdir, 'fastq')
    assemble_dir = os.path.join(new_outdir, "assembly_result")
    print("sra_dir:{}\nass_dir={}\nfastq_dir={}\nassemble_dir={}\n".format(sra_dir, ass_dir, fastq_dir, assemble_dir))

    try:
        # if sra_layout==2 continue
        #Download(sra_id,_outdir,sra_dir)
        Assembled(sra_id,_outdir,sra_dir,ass_dir,assemble_dir,fastq_dir,thread,gsize,start)
        check_dir = os.path.join(outdir, "check")
        utils_.mkdir_join(check_dir)
        check_Assemble=os.path.join(check_dir,"checkAssembled.txt")
        print("check_dir:{}\n".format(check_dir))

        with open(check_Assemble,"a+") as f:
            f.write("Run {} is ok.\n".format(sra_id))
        #####
        genome = os.path.join(ass_dir, "{}_contig.fa".format(sra_id))

        print("getBenga2\n")
        getBenga2(sra_id, outdir)
        print("getBenga2 Done\n")

        check_Benga = os.path.join(check_dir, "checkBenga.txt")
        with open(check_Benga,"a+") as f:
            f.write("Run {} is ok.\n".format(sra_id))


        targetPath=QualityCheck(sra_id,_outdir,outdir,genome,thread,gsize,start)
        time.sleep(1)

        check_QC = os.path.join(check_dir, "checkQC.txt")
        with open(check_QC,"a+") as f:
            f.write("Run {} is ok.\n".format(sra_id))

        print("targetPAth = {}\n######\n".format(targetPath.encode("utf-8").decode()))
        target_ = targetPath.replace(current_path, ".")
        print("target_= {}\n".format(target_))
        time.sleep(1)


        Analysis(sra_id,genome,target_,_outdir,thread,gsize,start)

        check_Analysis = os.path.join(check_dir, "checkAnalysis.txt")
        with open(check_Analysis,"a+") as f:
            f.write("Run {} is ok.\n".format(sra_id))
        #global finish_num_
        #finish_num_ += 1

        #check_log = os.path.join(check_dir, "Analysischeck.log")
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

        mycallback_write_Finish(sra_id)
        print("Run {} is Done\n".format(sra_id))
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
        SRA_run_error= os.path.join(check_dir, "SRA_run_error.txt")
        with open(SRA_run_error, "a+") as f:
            f.write("{}:{}\n".format(sra_id, errMsg))

        sys.exit(e)
    #sys.exit("subpreocess End\n")

    return sra_id

def mycallback_write_Finish(sra_id):
    print("mycallback_write\n")
    check_dir = os.path.join(outdir, "check")
    utils_.mkdir_join(check_dir)
    check_log = os.path.join(check_dir,"Analysischeck.log")
    print("check_log: {}\n".format(check_log))
    with open(check_log, "a+") as f:
        f.write("Run {} is ok.\n".format(sra_id))
        check_lines = f.readlines()
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, check_lines))
    finish_num_s = list(map(lambda x: x.split(" ")[1], finish))
    print("finish num={}\n".format(len(finish_num_s)))
    print("{} / {}\n".format(len(finish_num_s)+1, sra_num_))

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

    check_dir = os.path.join(outdir, "check")
    utils_.mkdir_join(check_dir)
    needList = os.path.join(check_dir, "need_run_{}.txt".format(txt_index))

    with open(needList,"r") as f:
        needlines=f.readlines()
    need_run=needlines
    need_run = [rr.strip() for rr in need_run if rr.strip()!='']
    print(need_run)
    sra_id=need_run[sra_index].strip("\n")
    check_log = os.path.join(check_dir, "Analysischeck.log")

    new_outdir = os.path.join(outdir, "output")
    utils_.mkdir_join(new_outdir)
    new_outdir = os.path.join(new_outdir, sra_id)
    utils_.mkdir_join(new_outdir)
    print("output: {}\n".format(new_outdir))

    #print(need_run)
    #print(need_run.type)
    print(sra_index)
    print(sra_id)

    SRA_Analysis(sra_id,  new_outdir, thread, gsize, start, sra_num_, outdir)

    #shutil.rmtree(sra_dir)
    #utils_.mkdir_join(fastq_dir)
    #shutil.rmtree(fastq_dir)
    #shutil.rmtree(assemble_dir)

