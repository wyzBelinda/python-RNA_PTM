import csv

import numpy as np
import pandas
import pandas as pd
from matplotlib import pyplot as plt


def plot_spectrum_matching(spectrum, rna, matched_peaks, root, title="Spectrum Matching"):
    # # 提取匹配的荷质比和强度
    # matched_mz = spectrum[matched_peaks, 0]
    # matched_intensity = spectrum[matched_peaks, 1]

    if len(matched_peaks) < 2:
        return

    # 绘制荷质比-强度柱形图
    fig, ax = plt.subplots()
    ax.bar(spectrum[:, 0], spectrum[:, 1], width=2, alpha=0.5, color='k')
    ax.bar(rna[:, 1], -3000000, width=2, alpha=0.5, color='b')
    ax.bar(matched_peaks[:, 0], matched_peaks[:, 1], width=2, alpha=0.5, color='r')

    # # 在谱峰上绘制标签
    # for i, mz in enumerate(spectrum[:, 0]):
    #     ax.text(mz, spectrum[i, 1],
    #             mz,
    #             # matched_peaks[i, 2] + "(" + matched_peaks[i, 3] + ")",
    #             ha='center',
    #             va='bottom')

    # # 在RNA计算好的碎片离子峰上绘制标签
    # for i, mz in enumerate(rna[:, 1]):
    #     ax.text(mz, 3000000,
    #             mz,
    #             # matched_peaks[i, 2] + "(" + matched_peaks[i, 3] + ")",
    #             ha='center',
    #             va='bottom')

    # 在匹配的峰上绘制标签
    for i, mz in enumerate(matched_peaks[:, 0]):
        ax.text(mz, matched_peaks[i, 1],
                "$" + matched_peaks[i, 2] + "$",
                # matched_peaks[i, 2] + "(" + matched_peaks[i, 3] + ")",
                ha='center',
                va='bottom')

    # 设置标题和坐标轴标签
    plt.title(title)
    plt.xlabel("m/z")
    plt.ylabel("Intensity")

    # 设置上边和右边无边框
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    # 设置x轴画在y=0
    ax.spines['bottom'].set_position(('data', 0))

    plt.savefig(root + title + ".png", dpi=800)
    # plt.show()

    plt.close(fig)


if __name__ == '__main__':
    root = "tests/results/02_AGUC/"
    # root = args.root
    # 读取csv文件
    re = csv.reader(open(root + "beat_matches.csv", "r"))

    for row in re:
        if row[0] == 'matched_score':
            continue
        print(row)
        try:
            s = pd.read_csv(root + "spectrum_" + str(row[0]) + ".csv")
            r = pd.read_csv(root + "oligo_ions_" + str(row[1]) + ".csv")
            m = pd.read_csv(root + "matched_" + str(row[0]) + "___" + str(row[1]) + ".csv")
        except pandas.errors.EmptyDataError:
            continue
        spectrum = np.array(s.loc[[0, 1]])
        rna = np.array(r.loc[[0, 1]])
        matched_peaks = np.array(m)
        plot_spectrum_matching(spectrum, rna, matched_peaks, root,
                               title=str(row[0]) + "___" + str(row[1]))
