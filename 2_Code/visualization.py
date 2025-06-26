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

波动方程
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# 参数设置
wave_speed = 1.0       
space_length = 1.0     
total_time = 2.0       
grid_points = 100      
time_steps = 400       

dx = space_length / (grid_points - 1)
dt = total_time / time_steps
stability_factor = wave_speed * dt / dx

if stability_factor > 1:
    raise ValueError("稳定性条件不满足！")

# 初始条件：高斯脉冲
x = np.linspace(0, space_length, grid_points)
initial_displacement = np.exp(-((x - 0.5)/0.15)**2)
initial_velocity = np.zeros_like(x)

# 计算前两层解
u_prev_prev = initial_displacement.copy()
u_prev = u_prev_prev + 0.5 * stability_factor**2 * (
    np.roll(u_prev_prev, -1) - 2*u_prev_prev + np.roll(u_prev_prev, 1)
)
u_prev[0], u_prev[-1] = 0, 0

# 预计算所有时间层
u = np.zeros((time_steps, grid_points))
u[0] = u_prev_prev
u[1] = u_prev

for n in range(1, time_steps-1):
    u_current = 2 * u[n] - u[n-1] + stability_factor**2 * (
        np.roll(u[n], -1, axis=-1) - 2*u[n] + np.roll(u[n], 1, axis=-1)
    )
    u_current[[0, -1]] = 0
    u[n+1] = u_current

# 动画配置
fig, ax = plt.subplots(figsize=(8, 4))
ax.set_xlim(0, space_length)
ax.set_ylim(-1.1, 1.1)
ax.set_xlabel('位置 x')
ax.set_ylabel('位移 u')
ax.set_title('波动方程动态解')
ax.grid(True, linestyle='--', alpha=0.7)

line, = ax.plot(x, u[0], 'r-', lw=2)

def animate(frame):
    line.set_ydata(u[frame])
    return line,

ani = FuncAnimation(fig, animate, frames=time_steps, interval=10, blit=True, repeat=True)
plt.show()
