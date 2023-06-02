import itertools

# 定义待优化的参数空间
import numpy as np
from matplotlib import pyplot as plt

from src.RNA_PTM import flows_continuing

param_grid = {'bm25_k': np.arange(0, 1.1, 0.1).tolist(),
              'bm25_t': np.arange(0, 1.1, 0.1).tolist()}


def my_func(x, y):
    # return x*y
    return flows_continuing(x, y)


# 定义目标函数
def objective_func(param1, param2):
    # 调用手写方法并返回结果
    result = my_func(param1, param2)
    return result


# 定义所有可能的参数组合
all_params = list(itertools.product(*param_grid.values()))

# 遍历每个参数组合，并记录最佳组合和最佳得分
best_params = None
best_score = float('-inf')
scores = []
for params in all_params:
    # 调用目标函数计算当前参数组合的得分
    score = objective_func(*params)
    print("score:", score)
    scores.append(score)

    # 如果当前得分比之前得到的最佳得分更高，则更新最佳得分和最佳参数组合
    if score > best_score:
        best_score = score
        best_params = params

# 输出最佳参数组合和最佳得分
print("Best parameters: ", best_params)
print("Best score: ", best_score)

# 将得分列表转换为矩阵
score_matrix = np.array(scores).reshape(len(param_grid['bm25_k']), len(param_grid['bm25_t']))
print(score_matrix)

# 绘制热图展示网格搜索结果
fig, ax = plt.subplots(figsize=(10, 8))
im = ax.imshow(score_matrix, cmap='jet')

# 添加坐标轴和标题
ax.set_xticks(np.arange(len(param_grid['bm25_k']))[::10])
# 等间距设置刻度标签
ax.set_xticklabels(param_grid['bm25_k'][::10])
ax.set_yticks(np.arange(len(param_grid['bm25_t']))[::10])
ax.set_yticklabels(param_grid['bm25_t'][::10])
ax.set_xlabel('bm25_k', fontsize=14)
ax.set_ylabel('bm25_t', fontsize=14)
ax.set_title('Grid Search Scores', fontsize=16)

# 添加颜色条
cbar = ax.figure.colorbar(im, ax=ax)
cbar.ax.set_ylabel('Score', rotation=-90, va='bottom')

# 在热图上标注最佳参数组合的位置
best_bm25_k, best_bm25_t = best_params
a = np.where(np.array(param_grid['bm25_k']) == best_bm25_k)
b = np.where(np.array(param_grid['bm25_t']) == best_bm25_t)
best_i, best_j = a[0][0], b[0][0]
ax.scatter(best_j, best_i, marker='*', color='white', s=100, linewidths=1, edgecolors='red', label='Best Parameters')
ax.legend()

plt.show()
