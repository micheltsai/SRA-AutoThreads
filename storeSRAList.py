import datetime
import os
import sys
import time

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
        utils_.prefetch_sra(x, sra_dir)
        print("Download {}\n.".format(x))
        #with open(Downloadcheck_log, "a+") as f:
        #    f.write("{}\n".format(x))
    dltime=time.time() - one_
    print('Done,total cost',dltime, 'secs')
    print("###########################################################")

def sra_stat(sra_id,outdir,sra_dir):

    QC_error = os.path.join(outdir, "nofillQC.txt")

    print("SequenceReadArchive\n")
    sra = utils_.SequenceReadArchivev2(sra_id)
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
        f.write(",")
    with open(sraList, "r") as f:
        print(f.readlines())

def main():
    start=time.time()
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
    thread = 4

    for yy in range(sd_Y,ed_Y+1):
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
        for mon in range(sM, eM+1):
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
            for d in range(sD,eD+1):


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
                f = open(sraList, 'r')
                line = f.readlines()
                print("check log :{}\n".format(line))
                f.close()

                for s in line:
                    print("{}\n".format(s))
                finish = list(filter(lambda x: len(x.split(",")) >= 4, line))
                finish_run = list(map(lambda x: x, finish))
                need_run = list(filter(lambda x: x not in finish_run, run_list))
                print("finish: {}\nfinish_run: {}\nneed_run".format(finish, finish_run, need_run))
                print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run),
                                                                                           len(need_run)))


                for aa in need_run:
                    sra_stat(aa,new_outdir,sra_dir)
                print('Done,total cost', time.time() - start, 'secs')



if __name__ == '__main__':
    main()