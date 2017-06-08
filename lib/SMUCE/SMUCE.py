
# run_max.py
import subprocess
import numpy as np
# import matplotlib.pyplot as plt
import csv


def SMUCE(seq):
    # Define command and arguments
    # write data
    with open('./lib/SMUCE/TS.csv', "w") as file_:
        writer = csv.writer(file_, delimiter='\n')
        writer.writerow(seq)

    command = 'Rscript'
    path2script = './lib/SMUCE/SMUCE.Rd'

    # Build subprocess command
    cmd = [command, path2script]  # + args

    # check_output will run the command and store to result
    x = subprocess.check_output(cmd, universal_newlines=True, stderr=subprocess.STDOUT)
    #print x
    cp = []
    with open('./lib/SMUCE/CP.csv', "r") as file_c:
        reader = csv.reader(file_c, delimiter=' ')
        i = 0
        for row in reader:
            if i > 1:
                cp.append(int(row[0]))
            i = i + 1

    return cp

if __name__ == "__main__":
    #seq = [3.,5.,4.,2.,.2,3.,4.,5.,1.,4.,3.,2.,2.,3.,3.,9.,8.,9.,8.,8.9,8.,9.,8.,8.,7.,7.]
    # plt.plot(seq)
    # plt.show()
    print PELT(seq)
