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

# 前两层解计算（显式有限差分法）
u_prev = u_prev_prev + 0.5 * stability_factor**2 * (
    np.roll(u_prev_prev, -1) - 2*u_prev_prev + np.roll(u_prev_prev, 1)
)  # 半步长公式

# 递推公式（二阶中心差分）
u_current = 2 * u[n] - u[n-1] + stability_factor**2 * (
    np.roll(u[n], -1, axis=-1) - 2*u[n] + np.roll(u[n], 1, axis=-1)
)  # 波动方程差分形式

# 边界条件：两端固定为0
u_prev[0], u_prev[-1] = 0, 0
u_current[[0, -1]] = 0
2.输运方程隐式格式求解、边界条件处理及稳定性验证研究

import numpy as np
import matplotlib.pyplot as plt
from numerical_methods import solve_heat_equation, build_coefficients

def plot_concentration_evolution(x, u_history, Nt, dt):
    """绘制浓度分布随时间的演化"""
    plt.figure(figsize=(10, 6))
    for i in range(0, Nt + 1, max(Nt // 10, 1)):
        plt.plot(x, u_history[i], label=f't={i * dt:.2f} s')
    plt.xlabel('位置 x (m)')
    plt.ylabel('浓度 u (kg/m$^3$)')
    plt.title('浓度分布随时间的演化（Crank-Nicolson 方法）')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('concentration_evolution.png')
    plt.show()

def pure_diffusion_analytical(x, t, D, L=1.0, sigma0=0.1):
    """纯扩散场景的解析解"""
    sig = np.sqrt(4 * D * t)
    if t == 0:
        return np.exp(-((x - L / 2)** 2) / (2 * sigma0** 2))
    return np.exp(-((x - L / 2)** 2) / (2 * D * t)) / (sig * np.sqrt(np.pi))

def mse(u_num, u_exact):
    """计算均方误差"""
    return np.mean((u_num - u_exact)** 2)

def plot_error_analysis(x, u0, Nx, dx, dt, T, D):
    """绘制纯扩散场景下的误差分析"""
    # 重新初始化并求解纯扩散问题（v=0）
    u = u0.copy()
    u_history = [u.copy()]
    v_pure = 0.0
    a_band_pure = build_coefficients(Nx, dx, dt, D, v_pure)
    
    for n in range(1, int(T / dt) + 1):
        b = u.copy()
        b[0] = 1
        b[-1] = 0
        solution = solve_banded((1, 1), a_band_pure, b)
        u[:] = solution
        u_history.append(u.copy())
    u_history_pure = np.array(u_history)
    
    # 计算MSE
    mse_list = []
    times = np.arange(0, T + dt, dt)
    for i, t in enumerate(times):
        u_exact = pure_diffusion_analytical(x, t, D)
        mse_val = mse(u_history_pure[i], u_exact)
        mse_list.append(mse_val)
    
    # 绘制误差曲线
    plt.figure(figsize=(8, 5))
    plt.plot(times, mse_list, label='MSE', color='blue')
    plt.xlabel('时间 t (s)')
    plt.ylabel('均方误差 (MSE) (kg/m$^3$)$^2$')
    plt.title('纯扩散场景下 Crank-Nicolson 方法的误差分析')
    plt.yscale('log')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('error_analysis_pure_diffusion.png')
    plt.show()

def plot_cfl_verification(x, u0, Nx, dx, T, D, v):
    """验证CFL条件对数值解的影响"""
    dt_values = [0.1, 0.01, 0.001]  # 包括不稳定的dt=0.1
    colors = ['red', 'blue', 'green']
    labels = [f'dt={dt_val} s' for dt_val in dt_values]
    
    plt.figure(figsize=(8, 5))
    for dt_val, color, label in zip(dt_values, colors, labels):
        # 重新初始化并求解
        u = u0.copy()
        u_history = [u.copy()]
        Nt_val = int(T / dt_val)
        a_band = build_coefficients(Nx, dx, dt_val, D, v)
        
        for n in range(1, Nt_val + 1):
            b = u.copy()
            b[0] = 1
            b[-1] = 0
            solution = solve_banded((1, 1), a_band, b)
            u[:] = solution
            u_history.append(u.copy())
        
        u_final = u_history[-1]
        plt.plot(x, u_final, color=color, label=label)
    
    plt.xlabel('位置 x (m)')
    plt.ylabel('浓度 u (kg/m$^3$)')
    plt.title('不同时间步长下的稳态浓度分布')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('cfl_condition_verification.png')
    plt.show()

def plot_nonhomogeneous_boundary(x, u0, Nx, dx, dt, T, D, v):
    """分析非齐次边界条件的影响"""
    def linear_dirichlet(t):
        return t  # 左边界浓度随时间线性增加
    
    # 重新初始化并求解非齐次边界问题
    u = u0.copy()
    u_history = [u.copy()]
    a_band = build_coefficients(Nx, dx, dt, D, v)
    
    for n in range(1, int(T / dt) + 1):
        t = n * dt
        b = u.copy()
        b[0] = linear_dirichlet(t)  # 左边界u=t
        b[-1] = 0
        solution = solve_banded((1, 1), a_band, b)
        u[:] = solution
        u_history.append(u.copy())
    u_history = np.array(u_history)
    
    # 绘制浓度演化
    plt.figure(figsize=(10, 6))
    for i in range(0, int(T / dt) + 1, max(int(T / dt) // 10, 1)):
        plt.plot(x, u_history[i], label=f't={i * dt:.2f} s')
    plt.xlabel('位置 x (m)')
    plt.ylabel('浓度 u (kg/m$^3$)')
    plt.title('非齐次 Dirichlet 边界条件下浓度分布随时间的演化')
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig('nonhomogeneous_dirichlet_evolution.png')
    plt.show()
