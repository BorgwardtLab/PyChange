import numpy as np
from PyChange import PyChange
from matplotlib import pyplot as plt

# plt.style.use('ggplot')
# plt.style.use('fivethirtyeight')
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=16)


if __name__ == "__main__":
    r = np.random.RandomState(42)
    seq = list(r.randn(100)) + list(r.randn(100) + 3) + list(r.randn(100)) + list(r.randn(100) + 3)
    f, ((ax1, ax2, ax7), (ax3, ax4, ax8), (ax5, ax6, ax9)) = plt.subplots(3, 3, sharex=True, sharey=True, figsize=(10, 10))

    axis = [ax1, ax2, ax3, ax4, ax5, ax9, ax6, ax8, ax7]
    methods = ['MaChaMP', 'PELT', 'WBS', 'SMUCE', 'E-Divise', 'BCP', 'Lepage', 'Segmentor3IsBack', 'Fpop']
    color = ['cornflowerblue', 'darkred', 'darkorange', 'olive', 'gold', 'teal', 'salmon', 'steelblue', 'rosybrown']

    for a, m, c in zip(axis, methods, color):
        a.plot(seq, color='black', alpha=0.5)
        loc = PyChange(seq, method=m)
        print m, loc
        if len(loc) >= len(seq) - 5:
            ax_b = a.twinx()
            ax_b.plot(loc, color=c, linewidth=2)
            ax_b.set_ylabel('posterior probability', color=c, fontsize=20)
            ax_b.tick_params('y', colors=c)
        else:
            for l in loc:
                a.axvline(l, color=c, linewidth=3)
        a.set_title(r"$\texttt{" + m + "}$", fontsize=20, color=c)
    ax6.set_xlabel('Time', fontsize=20)
    ax3.set_ylabel('Values', fontsize=20)
    ax1.set_xticks([100, 200, 300])
    ax1.set_xlim((0, 400))
    plt.tight_layout()
    plt.savefig('Example.pdf')
