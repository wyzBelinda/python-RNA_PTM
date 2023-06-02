import csv
import re

import numpy as np
import pandas
import pandas as pd
from matplotlib import pyplot as plt


#
# def gen_intensity_matrix(oligo_seq):
#     matrix = [[ -1 for i in range(n)] for j in range(m)]

def draw_peak_match(ax, title, spectrum, oligo_ions_mz, match_infos_list):
    ax.bar(spectrum[:, 0], spectrum[:, 1], width=2, alpha=0.5, color='k')
    ax.bar(oligo_ions_mz, -3000000, width=2, alpha=0.5, color='b')

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

    for x, y, info in match_infos_list:
        ax.bar(x, y, width=2, alpha=0.5, color='r')
        ax.text(x, y, "$" + info + "$", ha='center', va='bottom')

    # 设置标题和坐标轴标签
    ax.set_title(title, y=1.1)
    ax.set_xlabel("m/z")
    ax.set_ylabel("Intensity")

    # 设置上边和右边无边框
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')
    # 设置x轴画在y=0
    ax.spines['bottom'].set_position(('data', 0))

    return ax


def draw_match_ion_info(ax, sequence, match_inten_matrix, match_charge_matrix):
    # 绘制序列图
    im = ax.imshow(match_inten_matrix, cmap='jet')

    for _i in range(len(match_charge_matrix)):
        for _j in range(len(match_charge_matrix[0])):
            ax.text(_j, _i, match_charge_matrix[_i][_j],
                    ha="center", va="center", color="w")

    # print(seq.split(""))
    # 添加坐标轴和标题
    ax.set_xticks(np.arange(len(sequence)))
    # 等间距设置刻度标签
    ax.set_xticklabels([str(c) for c in sequence])
    ax.set_yticks(np.arange(9))
    ax.set_yticklabels(['a-B', 'a', 'b', 'c', 'd', 'w', 'x', 'y', 'z'])
    ax.set_xlabel('sequence', fontsize=14)
    ax.set_ylabel('ion_type', fontsize=14)
    ax.set_title('Ion Peak Matches Intensity', fontsize=14, y=1.05)

    # 添加颜色条
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Intensity', rotation=-90, va='bottom')

    return ax


def plot_spectrum_matching_test(spectrum, oligo_ions, matched_peaks, root, title="Spectrum Matching", sequence=""):
    if len(matched_peaks) < 2:
        return

    match_inten_matrix = [[-1 for i in range(len(sequence))] for j in range(9)]
    match_charge_matrix = [[0 for i in range(len(sequence))] for j in range(9)]

    match_infos_list = []
    for i, j, _ in matched_peaks:
        # 在匹配的峰上绘制标签
        i = int(i)
        j = int(j)
        match_ion_info = oligo_ions[j, 0]

        match_infos_list.append([spectrum[i][0], spectrum[i][1], match_ion_info])

        ion_type = re.findall(r'(\w)', match_ion_info)[0]
        is_B = re.findall(r'-\w', match_ion_info)
        ion_pos = int(re.findall(r'_(\d+)', match_ion_info)[0])
        charge = re.findall(r'\^(\d+)', match_ion_info)[0]

        if is_B:
            x = 0
            y = ion_pos
        elif ord(ion_type) < ord('e'):
            x = int(ord(ion_type) - ord('a')) + 1
            y = ion_pos
        else:
            x = int(ord(ion_type) - ord('w')) + 5
            y = len(sequence) - ion_pos
        match_inten_matrix[x][y] = spectrum[i][1]
        match_charge_matrix[x][y] = charge

    # 绘制荷质比-强度柱形图
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    draw_peak_match(ax1, title, spectrum, oligo_ions[:, 1], match_infos_list)
    draw_match_ion_info(ax2, sequence, match_inten_matrix, match_charge_matrix)

    # plt.savefig(root + title + ".png", dpi=800)
    plt.show(dpi=800)

    plt.close(fig)

def plot_spectrum_matching(spectrum, oligo_ions, matched_peaks, root, title="Spectrum Matching", sequence=""):
    # # 提取匹配的荷质比和强度
    # matched_mz = spectrum[matched_peaks, 0]
    # matched_intensity = spectrum[matched_peaks, 1]

    if len(matched_peaks) < 2:
        return

    match_inten_matrix = [[-1 for i in range(len(sequence))] for j in range(9)]
    match_charge_matrix = [[0 for i in range(len(sequence))] for j in range(9)]

    # 绘制荷质比-强度柱形图
    fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    ax1.bar(spectrum[:, 0], spectrum[:, 1], width=2, alpha=0.5, color='k')
    ax1.bar(oligo_ions[:, 1], -3000000, width=2, alpha=0.5, color='b')

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

    for i, j, _ in matched_peaks:
        # 在匹配的峰上绘制标签
        i = int(i)
        j = int(j)
        match_ion_info = oligo_ions[j, 0]
        ax1.bar(spectrum[i][0], spectrum[i][1], width=2, alpha=0.5, color='r')
        ax1.text(spectrum[i][0], spectrum[i][1],
                 "$" + match_ion_info + "$",
                 # matched_peaks[i, 2] + "(" + matched_peaks[i, 3] + ")",
                 ha='center',
                 va='bottom')

        ion_type = re.findall(r'(\w)', match_ion_info)[0]
        is_B = re.findall(r'-\w', match_ion_info)
        ion_pos = int(re.findall(r'_(\d+)', match_ion_info)[0])
        charge = re.findall(r'\^(\d+)', match_ion_info)[0]

        if is_B:
            x = 0
            y = ion_pos
        elif ord(ion_type) < ord('e'):
            x = int(ord(ion_type) - ord('a')) + 1
            y = ion_pos
        else:
            x = int(ord(ion_type) - ord('w')) + 5
            y = len(sequence) - ion_pos
        match_inten_matrix[x][y] = spectrum[i][1]
        match_charge_matrix[x][y] = charge

    # 设置标题和坐标轴标签
    ax1.set_title(title, y=1.1)
    ax1.set_xlabel("m/z")
    ax1.set_ylabel("Intensity")

    # 设置上边和右边无边框
    ax1.spines['right'].set_color('none')
    ax1.spines['top'].set_color('none')
    # 设置x轴画在y=0
    ax1.spines['bottom'].set_position(('data', 0))

    # ----绘制序列匹配示意图------
    ax = ax2

    # 绘制序列图
    im = ax.imshow(match_inten_matrix, cmap='jet')

    for _i in range(len(match_charge_matrix)):
        for _j in range(len(match_charge_matrix[0])):
            ax.text(_j, _i, match_charge_matrix[_i][_j],
                    ha="center", va="center", color="w")

    # print(seq.split(""))
    # 添加坐标轴和标题
    ax.set_xticks(np.arange(len(sequence)))
    # 等间距设置刻度标签
    ax.set_xticklabels([str(c) for c in sequence])
    ax.set_yticks(np.arange(9))
    ax.set_yticklabels(['a-B', 'a', 'b', 'c', 'd', 'w', 'x', 'y', 'z'])
    ax.set_xlabel('sequence', fontsize=14)
    ax.set_ylabel('ion_type', fontsize=14)
    ax.set_title('Ion Peak Matches Intensity', fontsize=14, y=1.05)

    # 添加颜色条
    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel('Intensity', rotation=-90, va='bottom')

    # plt.savefig(root + title + ".png", dpi=800)
    plt.show(dpi=800)

    plt.close(fig)


def best_matchings_plot(root, best_match_filename):
    # 读取csv文件
    re = csv.reader(open(root + best_match_filename, "r"))

    for row in re:
        if row[0] == 'matched_score':
            continue
        print(row)
        try:
            s = pd.read_csv(root + "spectrum_" + str(row[3]).replace(" ", "_").replace(",", "_").replace(".", "_")
                            .replace('+', '_').replace("/", "_") + ".csv")
            # 左右两端本应该自定义，这里设默认了
            chem5 = "O(1)H(1)"
            chem3 = "O(1)H(1)"
            o = pd.read_csv(root + "oligo_ions_" + chem5 + "_" + str(row[1]) + "_" + chem3 + ".csv")
            m = pd.read_csv(
                root + "peak_match_infos_" + str(row[3]).replace(" ", "_").replace(",", "_").replace(".", "_")
                .replace('+', '_').replace("/", "_") + "_" + str(row[1]) + ".csv")
        except FileNotFoundError:
            continue
        spectrum = s.values
        oligo_ions = o.values
        matched_peaks = np.array(m)
        # plot_spectrum_matching(spectrum, oligo_ions, matched_peaks, root,
        #                        title="plot_best_match_" + str(row[3]) + "___" + str(row[1]), sequence=row[1])
        plot_spectrum_matching_test(spectrum, oligo_ions, matched_peaks, root,
                                       title="plot_best_match_" + str(row[3]) + "___" + str(row[1]), sequence=row[1])


if __name__ == '__main__':
    root = "F:\\results\\python-RNA_PTM\\sun_zuo\\"  # 02_OH-AGUC-OH.379.-1.CID_10
    # root = args.root

    # best_matchings_plot(root, "best_matches.csv")
    best_matchings_plot(root, "best_matches_false.csv")
