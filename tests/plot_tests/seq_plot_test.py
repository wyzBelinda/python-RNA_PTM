import matplotlib.pyplot as plt

# 序列及其离子位置信息
sequence = "MKALLEDPQIAAQAR"
ion_locations = {"d2": 2, "b2": 2, "c2": 2,"a2": 2, 'a2-B': 2, "w3": 11, "x3": 11, "y3": 11, "z3": 11}

# 绘制序列图
height_center = 2
fig, ax = plt.subplots(figsize=(len(sequence) * 0.5, 1))
ax.set_xlim([0, len(sequence)])
ax.set_ylim([0, height_center * 2+1])
ax.axis('off')

# 绘制序列
for i in range(len(sequence)):
    ax.text(i + 0.5, height_center, sequence[i], ha='center', va='center', fontsize=12)

# 绘制离子标记
for ion_name, loc in ion_locations.items():
    if ion_name[0] == 'a' and '-'in ion_name:
        ax.plot([loc - 1, loc - 1], [height_center - 0.4, height_center + 0.4], color='blue')
        ax.plot([loc - 1 - 0.2, loc - 1], [height_center + 0.4, height_center + 0.4], color='blue')
        ax.text((loc - 1 + loc) / 2 - 1, height_center + 0.4, ion_name, ha='center', va='bottom', color='blue',
                fontsize=8)
    elif ion_name[0] == 'a':
        ax.plot([loc - 1, loc - 1], [height_center - 0.4, height_center + 1], color='blue')
        ax.plot([loc - 1 - 0.2, loc - 1], [height_center + 1, height_center + 1], color='blue')
        ax.text((loc - 1 + loc) / 2 - 1, height_center + 1, ion_name, ha='center', va='bottom', color='blue',
                fontsize=8)
    elif ion_name[0] == 'b':
        ax.plot([loc - 1, loc - 1], [height_center - 0.4, height_center +1.6], color='blue')
        ax.plot([loc - 1 - 0.2, loc - 1], [height_center + 1.6, height_center + 1.6], color='blue')
        ax.text((loc - 1 + loc) / 2 - 1, height_center + 1.6, ion_name, ha='center', va='bottom', color='blue',
                fontsize=8)
    elif ion_name[0] == 'c':
        ax.plot([loc - 1, loc - 1], [height_center - 0.4, height_center + 2.2], color='blue')
        ax.plot([loc - 1 - 0.2, loc - 1], [height_center + 2.2, height_center + 2.2], color='blue')
        ax.text((loc - 1 + loc) / 2 - 1, height_center + 2.2, ion_name, ha='center', va='bottom', color='blue',
                fontsize=8)
    elif ion_name[0] == 'd':
        ax.plot([loc - 1, loc - 1], [height_center - 0.4, height_center + 2.8], color='blue')
        ax.plot([loc - 1 - 0.2, loc - 1], [height_center + 2.8, height_center + 2.8], color='blue')
        ax.text((loc - 1 + loc) / 2 - 1, height_center + 2.8, ion_name, ha='center', va='bottom', color='blue',
                fontsize=8)
    elif ion_name[0] == 'w':
        ax.plot([len(sequence) - loc, len(sequence) - loc], [height_center - 2, height_center + 0.5], color='red')
        ax.plot([len(sequence) - loc, len(sequence) - loc + 0.2], [height_center-0.5, height_center-0.5], color='red')
        ax.text((len(sequence) - loc + len(sequence) - loc + 1) / 2, height_center - 0.5 - 0.6, ion_name, ha='center',
                va='bottom',
                color='red',  fontsize=8)
    elif ion_name[0] == 'x':
        ax.plot([len(sequence) - loc, len(sequence) - loc], [height_center - 2, height_center + 0.5], color='red')
        ax.plot([len(sequence) - loc, len(sequence) - loc + 0.2], [height_center-0.5- 0.6, height_center-0.5- 0.6], color='red')
        ax.text((len(sequence) - loc + len(sequence) - loc + 1) / 2, height_center - 0.5 - 1.2, ion_name, ha='center',
                va='bottom',
                color='red', fontsize=8)
    elif ion_name[0] == 'y':
        ax.plot([len(sequence) - loc, len(sequence) - loc], [height_center - 2, height_center + 0.5], color='red')
        ax.plot([len(sequence) - loc, len(sequence) - loc + 0.2], [height_center-0.5- 1.2, height_center-0.5- 1.2], color='red')
        ax.text((len(sequence) - loc + len(sequence) - loc + 1) / 2, height_center - 0.5 - 1.8, ion_name, ha='center',
                va='bottom',
                color='red',  fontsize=8)
    else:
        ax.plot([len(sequence) - loc, len(sequence) - loc], [height_center-2.4, height_center+0.5], color='red')
        ax.plot([len(sequence) - loc, len(sequence) - loc + 0.2], [height_center-0.5- 1.8, height_center-0.5- 1.8], color='red')
        ax.text((len(sequence) - loc + len(sequence) - loc + 1) / 2, height_center-0.5-2.4, ion_name, ha='center', va='bottom',
                color='red', fontsize=8)

plt.show()

plt.close()

# import matplotlib.pyplot as plt
#
# # 在同一画布上创建两个子图
# fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)
#
# # 设置第一个子图的属性（即您已经创建的柱形图）
# ax1.bar([1, 2, 3], [4, 5, 6])
# ax1.set_xlabel('X Label')
# ax1.set_ylabel('Y Label')
# ax1.set_title('Bar Chart')
#
# # 设置第二个子图的属性（即您要展示的图表）
# ax2.text(0.5, 0.5, 'Your Other Chart', ha='center', va='center', fontsize=18)
#
# # 调整第二个子图的位置和大小
# pos = ax1.get_position()
# # pos.width = pos.width*2
# ax2.set_position([pos.x0 + pos.width * 0.7, pos.y0 + pos.height * 0.7, 0.2, 0.2])
#
# plt.show()
#
# #
# import matplotlib.pyplot as plt
# from Bio.SeqUtils.ProtParam import ProteinAnalysis
#
# # 序列及其离子位置信息
# sequence = "MKALLEDPQIAAQAR"
# ion_locations = {"b2":2, "y3":11}
#
# # 创建ProteinAnalysis对象，计算序列的各种属性
# protein = ProteinAnalysis(sequence)
# charge = protein.charge_at_pH(7.0)
# isoelectric_point = protein.isoelectric_point()
#
# # 设置画布和轴属性
# fig, ax = plt.subplots(figsize=(len(sequence)*0.5, 1))
# ax.set_xlim([0, len(sequence)])
# ax.set_ylim([0, 1])
# ax.spines['top'].set_visible(False)
# ax.spines['right'].set_visible(False)
# ax.spines['bottom'].set_visible(False)
# ax.spines['left'].set_visible(False)
# ax.tick_params(axis='both', length=0)
#
# # 绘制序列
# for i in range(len(sequence)):
#     if sequence[i] in "AGVIL":
#         color = 'darkgray'
#     elif sequence[i] in "FYW":
#         color = 'blue'
#     else:
#         color = 'black'
#     ax.text(i+0.5, 0.5, sequence[i], ha='center', va='center', fontsize=12, color=color)
#
# # 绘制离子标记
# for ion_name, loc in ion_locations.items():
#     if ion_name[0] == 'b':
#         ax.plot([loc-1, loc], [0.2, 0.8], color='blue', linewidth=2)
#         ax.text((loc-1+loc)/2, 0.9, ion_name, ha='center', va='bottom', color='blue', fontsize=8)
#     else:
#         ax.plot([len(sequence)-loc, len(sequence)-loc+1], [0.2, 0.8], color='red', linewidth=2)
#         ax.text((len(sequence)-loc+len(sequence)-loc+1)/2, 0.9, ion_name, ha='center', va='bottom', color='red', fontsize=8)
#
# # 添加序列长度标签
# ax.text(len(sequence)+0.5, 0.5, str(len(sequence)), ha='center', va='center', fontsize=12, color='gray')
#
# # 添加氨基酸属性标签
# ax.text(0.5, 1.2, f'Charge: {charge:.2f}', ha='left', va='bottom', fontsize=8)
# ax.text(len(sequence)/2, 1.2, f'Isoelectric Point: {isoelectric_point:.2f}', ha='center', va='bottom', fontsize=8)
#
# plt.show()
