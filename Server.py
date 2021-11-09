import os
import time
import datetime
from pathlib import Path
from socket import *
import socketserver

import numpy as np

import utils_
import pandas as pd

global current_path, start_date, expiry_date, thread, cpu_process, gsize, outdir, ref_dir, buscoDB, buscoMode, mlstS, amrS
### sd_Y/sd_M/sd_D - ed_Y/ed_M/ed_d
global sd_Y, sd_M, sd_D, ed_Y, ed_M, ed_D
##############
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
######################


class MyTCPHandler(socketserver.BaseRequestHandler):
    def handle(self):
        self.data = self.request.recv(1024).strip()
        print("{} wrote:".format(self.client_address[0]))
        print(self.data)
        self.request.sendall(self.data.upper())

class SplitSRAIdList(socketserver.BaseRequestHandler):
    def handle(self):
        self.data = self.request.recv(1024).strip()
        print("{} wrote:".format(self.client_address[0]))
        print(self.data)
        self.request.sendall(self.data.upper())
        return self.data

def getSRAIdList(yy,mon,d):
    pattern = "salmonella enterica[ORGN] AND illumina[PLAT] AND wgs[STRA] AND genomic[SRC] AND paired[LAY]"
    ds = time.time()

    date = datetime.date(yy, mon, d).strftime("%Y/%m/%d")
    # temp="{}/{}/{}".format(str(2020),str(mon+1),str(d))
    ######
    pdat = date.replace("/", "")
    new_outdir = os.path.join(outdir, pdat)
    utils_.mkdir_join(new_outdir)
    print("output: {}\n".format(new_outdir))

    # commit
    check_log = os.path.join(new_outdir, "Analysischeck.log")

    #myfile2 = Path(check_log)
    #myfile2.touch(exist_ok=True)
    #with open(check_log, "a+") as f:
    #    f.write(str(datetime.datetime.now()).split(".")[0])
    #    f.write("\n")

    pattern, count = utils_.count_egquery(pattern, date, date)
    print("pattern: {}\ncount: {}\n".format(pattern, count))

    i_e_ = time.time()
    idlist = utils_.IdList_esearch(pattern, 'sra', count)

    print(idlist)

    runinfo = utils_.Get_RunInfo(idlist)
    run_list = list(runinfo['Run'])  # get SRAfile nameList stored in run_list
    print("runinfo: {}\n run_list: {}\n".format(runinfo, run_list))

    sra_dir = os.path.join(new_outdir, "sra")  # .sra file
    utils_.mkdir_join(sra_dir)
    ass_dir = os.path.join(new_outdir, "Assembled")
    utils_.mkdir_join(ass_dir)
    fastq_dir = os.path.join(new_outdir, 'fastq')
    utils_.mkdir_join(fastq_dir)
    assemble_dir = os.path.join(new_outdir, "assembly_result")
    utils_.mkdir_join(assemble_dir)
    print("sra_dir:{}\nass_dir={}\nfastq_dir={}\nassemble_dir={}\n".format(sra_dir, ass_dir, fastq_dir, assemble_dir))

    f = open(check_log, 'r')
    line = f.readlines()
    print("check log :{}\n".format(line))
    f.close()

    for s in line:
        print("{}\n".format(s))
    finish = list(filter(lambda x: len(x.split(" ")) >= 4, line))
    finish_run = list(map(lambda x: x.split(" ")[1], finish))
    need_run = list(filter(lambda x: x not in finish_run, run_list))
    print("finish: {}\nfinish_run: {}\nneed_run".format(finish, finish_run, need_run))
    print("finish length: {}\nfinish_run length: {}\nneed_run length: ".format(len(finish), len(finish_run),
                                                                               len(need_run)))
    print("Toal", len(need_run), "sra runs need to downlaod.")
    need_run_=np.array_split(need_run,2)
    need_run_txt=os.path.join(new_outdir,"run_1.txt")
    need_run_txt2 = os.path.join(new_outdir, "run_2.txt")
    with open(need_run_txt,"a+") as f:
        f.write(date)
        for n in need_run_[0]:
            f.write(n)
    with open(need_run_txt2,"a+") as f:
        f.write(date)
        for n in need_run_[1]:
            f.write(n)

if __name__ == "__main__":
    HOST, PORT = '140.112.165.124', 8088
    # Create the server, binding to localhost on port 9999

    getSRAIdList(2020,8,1)
    server = socketserver.TCPServer((HOST, PORT), MyTCPHandler)

    # Activate the server; this will keep running until you
    # interrupt the program with Ctrl-C
    server.serve_forever()