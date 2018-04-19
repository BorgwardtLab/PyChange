import numpy as np
from PyChange import PyChange
from matplotlib import pyplot as plt
import csv
import matplotlib
import datetime


# plt.style.use('ggplot')
# plt.style.use('fivethirtyeight')
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)


if __name__ == "__main__":
    r = np.random.RandomState(42)

    seq = []
    dates = []
    with open('ftse100.csv', "r") as file_c:
        reader = csv.reader(file_c, delimiter=',')
        i = 0
        for row in reader:
            if i != 0:
                seq.append(float(row[1]))
                dates.append(row[0])
            i = i + 1
    #dates = matplotlib.dates.datestr2num(dates)
    dates = map(datetime.datetime.strptime, dates, len(dates) * ["%Y-%m-%d"])

    plt.xlabel('year', fontsize=18)
    plt.ylabel('daily returns FTSE100 [\%]', fontsize=18)
    plt.plot(dates, seq, color='k')
    plt.tight_layout()
    plt.savefig('FTSE100_raw.pdf')
    plt.clf()

    f, ((ax1, ax2, ax7), (ax3, ax4, ax8), (ax5, ax6, ax9)) = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(10, 10))

    axis = [ax1, ax2, ax7, ax3, ax4, ax8, ax5, ax6, ax9]
    methods = ['PELT', 'PELT', 'PELT', 'CUSUM', 'CUSUM', 'CUSUM', 'BCP', 'BCP', 'BCP']  # ['MaChaMP', 'PELT', 'WBS', 'SMUCE', 'E-Divise', 'BCP', 'Lepage', 'EWMA', 'Fpop']
    color = ['darkblue', 'darkred', 'darkorange', 'olive', 'gold', 'teal', 'salmon', 'steelblue', 'rosybrown']  # ['darkblue', 'darkred', 'darkorange', 'olive', 'gold', 'teal', 'salmon', 'steelblue', 'rosybrown']
    transformation = ['std', 'diff', 'logdiff', 'std', 'diff', 'logdiff', 'std', 'diff', 'logdiff']

    for a, m, c, t in zip(axis, methods, color, transformation):
        a.plot(seq, color='black', alpha=0.5)

        if m != 'none':
            loc = PyChange(seq, transform=t, method=m)
        else:
            loc = []
        print m, loc
        if len(loc) >= len(seq) - 5:
            ax_b = a.twinx()
            ax_b.plot(loc, color=c, linewidth=1)
            if t == 'logdiff':
                ax_b.set_ylabel('posterior probability', fontsize=20)
            ax_b.tick_params('y')
        else:
            for l in loc:
                a.axvline(l, color=c, linewidth=1)
        #a.set_title(r"$\texttt{" + m + "}$", fontsize=20, color=c)
        a.set_title(m + ' ' + t, fontsize=20, color=c)
    ax6.set_xlabel('Time', fontsize=20)
    ax3.set_ylabel('Values', fontsize=20)
    #ax1.set_xticks([100, 200, 300])
    ax1.set_xlim((0, len(seq)))
    plt.tight_layout()
    plt.savefig('FTSE100.pdf')
