import numpy as np
from matplotlib import pyplot as plt

param_grid = {'bm25_k': np.arange(0.3, 2.1, 0.1).tolist(),
              'bm25_t': np.arange(0.5, 5.6, 0.5).tolist()}

matrix = [[147, 156, 161, 163, 164, 168, 172, 170, 166, 141, 160, 150, 156, 160, 164, 165, 166, 169],
          [172, 165, 126, -1, 152, 157, 161, 164, 165, 170, 169, 171, 167, 144, 167, 153, 157, 161],
          [164, 165.5, 167, 172, 171, 168, 167, 168, 154, 158, 162, 164, 165, 169, 169, 171, 168, 167],
          [168, 154, 158, 162, 164, 167, 171, 169, 172, 168, 168, 168, 153, 159, 163, 164, 167, 170],
          [169, 172, 165, 168, 168, 154, 159, 163, 164, 168, 166, 169, 172, 163, 168, 169, 155, 160],
          [163, 164, 168, 167, 170.666, 172, 159, 168, 169, 155, 160, 163, 164, 168, 166, 170, 172, 155],
          [168, 168, 155, 160, 163, 164, 168, 167, 170, 172, 152, 168, 168, 155, 161, 163, 164, 168],
          [169, 170, 172, 153, 168, 168, 155, 161, 163, 164, 168, 168, 170, 172, 153, 168, 168, 155],
          [161, 163, 164, 168, 167, 170, 172, 157, 168, 168, 156, 161, 163, 164, 168, 167, 170, 172]]
# 绘制热图展示网格搜索结果
fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(matrix, cmap='jet')

seq="AZSXDCFVGBHJMKL"
# print(seq.split(""))
# 添加坐标轴和标题
ax.set_xticks(np.arange(len(seq)))
# 等间距设置刻度标签
ax.set_xticklabels([str(c) for c in seq])
ax.set_yticks(np.arange(9))
ax.set_yticklabels(['a-B', 'a', 'b', 'c', 'd', 'w', 'x', 'y', 'z'])
ax.set_xlabel('sequence', fontsize=14)
ax.set_ylabel('ion_type', fontsize=14)
ax.set_title('Ion Peak Matches Intensity', fontsize=16)

# 添加颜色条
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel('Intensity', rotation=-90, va='bottom')

# # 在热图上标注最佳参数组合的位置
# best_bm25_k, best_bm25_t = best_params
# a = np.where(np.array(param_grid['bm25_k']) == best_bm25_k)
# b = np.where(np.array(param_grid['bm25_t']) == best_bm25_t)
# best_i, best_j = a[0][0], b[0][0]
# ax.scatter(best_j, best_i, marker='*', color='white', s=100, linewidths=1, edgecolors='red', label='Best Parameters')
ax.legend()

plt.show()
