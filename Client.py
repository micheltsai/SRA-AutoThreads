import os
import socket
import sys

import pandas as pd

import utils_
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
outdir = str(setting_df['output_dir'][0])

if __name__ == "__main__":
    run_txt_path = os.path.join(outdir,"run_2.txt")
    while True:
        if os.path.isfile(run_txt_path):
            print("get run_2.txt")

        else:
            os.sleep(10)

