import numpy as np
from PyChange import PyChange
from matplotlib import pyplot as plt

# plt.style.use('ggplot')
# plt.style.use('fivethirtyeight')
# plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=14)


if __name__ == "__main__":
    r = np.random.RandomState(42)
    seq = [1120, 1160, 963, 1210, 1160, 1160, 813, 1230, 1370, 1140, 995, 935, 1110, 994, 1020, 960, 1180, 799, 958, 1140, 1100, 1210, 1150, 1250, 1260, 1220, 1030, 1100, 774, 840, 874, 694, 940, 833, 701, 916, 692, 1020, 1050, 969, 831, 726, 456, 824, 702, 1120, 1100, 832, 764, 821, 768, 845, 864, 862, 698, 845, 744, 796, 1040, 759, 781, 865, 845, 944, 984, 897, 822, 1010, 771, 676, 649, 846, 812, 742, 801, 1040, 860, 874, 848, 890, 744, 749, 838, 1050, 918, 986, 797, 923, 975, 815, 1020, 906, 901, 1170, 912, 746, 919, 718, 714, 740]

    plt.xlabel('year', fontsize=18)
    plt.ylabel('Nile flow volume [$m^3/y$]', fontsize=18)
    ax = plt.gca()
    ax.plot(range(1871, 1971), seq)
    ax.annotate('Is this a changepoint?', xy=(1900, 1200), xytext=(1920, 1300), arrowprops=dict(facecolor='black', shrink=0.05))
    plt.tight_layout()
    plt.savefig('Nile_raw.pdf')
    plt.clf()

    f, ((ax1, ax2, ax7), (ax3, ax4, ax8), (ax5, ax6, ax9)) = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(10, 10))

    axis = [ax1, ax2, ax7, ax4, ax3, ax8, ax6, ax9, ax5]
    methods = ['PELT', 'MaChaMP', 'FPOP', 'CUSUM', 'EWMA', 'QChart', 'E-Divise', 'SMUCE', 'BCP']  # ['MaChaMP', 'PELT', 'WBS', 'SMUCE', 'E-Divise', 'BCP', 'Lepage', 'EWMA', 'Fpop']
    color = ['darkblue', 'darkred', 'darkorange', 'olive', 'gold', 'teal', 'salmon', 'steelblue', 'rosybrown']  # ['darkblue', 'darkred', 'darkorange', 'olive', 'gold', 'teal', 'salmon', 'steelblue', 'rosybrown']

    for a, m, c in zip(axis, methods, color):
        a.plot(seq, color='black', alpha=0.5)
        if m != 'none':
            loc = PyChange(seq, method=m)
        else:
            loc = []
        if len(loc) >= len(seq) - 5:
            ax_b = a.twinx()
            ax_b.plot(loc, color=c, linewidth=2)
            ax_b.set_ylabel('posterior probability', color=c, fontsize=20)
            ax_b.tick_params('y', colors=c)
        else:
            for l in loc:
                a.axvline(l, color=c, linewidth=3)
        #a.set_title(r"$\texttt{" + m + "}$", fontsize=20, color=c)
        a.set_title(m, fontsize=20, color=c)
    ax6.set_xlabel('Time', fontsize=20)
    ax3.set_ylabel('Values', fontsize=20)
    #ax1.set_xticks([100, 200, 300])
    ax1.set_xlim((0, len(seq)))
    plt.tight_layout()
    plt.savefig('Nile.pdf')
