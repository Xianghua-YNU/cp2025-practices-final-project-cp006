"""
可视化函数模块：用于绘制模拟结果和分析数据。
"""
#芯片热传导问题的数值求解与方法对比研究
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['mathtext.fontset'] = 'cm'

def plot_temperature(T, x, y, dx, method_name):
    """绘制温度分布与热流场，增加温度棒"""
    fig = plt.figure(figsize=(14, 6))  # 增加画布宽度以容纳温度棒

    # 3D温度分布
    ax1 = fig.add_subplot(121, projection='3d')
    X, Y = np.meshgrid(x, y)
    surf = ax1.plot_surface(X, Y, T, cmap=cm.jet, alpha=0.8)
    
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_zlabel('温度 (K)')
    ax1.set_title(f'芯片温度分布\n({method_name})')
    
    # 添加温度棒，设置位置和大小
    cbar_ax = fig.add_axes([0.91, 0.15, 0.02, 0.7])  # [left, bottom, width, height]
    fig.colorbar(surf, cax=cbar_ax, label='温度 (K)')

    # 温度等值线与热流场
    ax2 = fig.add_subplot(122)
    X, Y = np.meshgrid(x, y)
    # 计算热流场（负温度梯度）
    Ty, Tx = np.gradient(-T, dx)
    # 绘制温度等值线
    levels = np.linspace(T.min(), T.max(), 15)
    contour = ax2.contour(X, Y, T, levels=levels, colors='red', linestyles='solid')
    ax2.clabel(contour, inline=True, fontsize=6)
    # 绘制热流场流线
    ax2.streamplot(X, Y, Tx, Ty, density=1.5, color='blue', linewidth=1, arrowsize=1.5)
    # 标记热源位置
    heat_source = 100 * np.exp(-((X - 0.5) ** 2 + (Y - 0.7) ** 2) / (2 * 0.1 ** 2))
    ax2.contour(X, Y, heat_source, levels=[5], colors='black', linestyles='dashed')
    ax2.text(0.5, 0.7, '热源', fontsize=9, ha='center')

    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    ax2.set_title(f'温度等值线与热流场\n({method_name})')
    ax2.set_aspect('equal')
    
    # 为2D图添加温度棒
    cbar_ax2 = fig.add_axes([0.91, 0.15, 0.02, 0.3])
    fig.colorbar(contour, cax=cbar_ax2, label='温度 (K)')

    plt.tight_layout()
    plt.show()
