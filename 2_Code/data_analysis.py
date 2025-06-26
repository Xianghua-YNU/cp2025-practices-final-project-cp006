"""
数据后处理与分析模块：用于对模拟结果进行分析。
"""
#芯片热传导问题的数值求解与方法对比研究
import numpy as np

def calculate_error(T1, T2):
    """计算两个温度场之间的误差"""
    return np.max(np.abs(T1 - T2))

def analyze_convergence(history, method_name):
    """分析收敛历史"""
    iterations = len(history)
    final_error = history[-1]
    
    print(f"\n{method_name}收敛分析:")
    print(f"总迭代次数: {iterations}")
    print(f"最终误差: {final_error:.6e}")
    
   


# 输运方程
# 构建系数矩阵
def build_coefficients(Nx, dx, dt, D, v):
    alpha = D * dt / dx ** 2
    beta = v * dt / (2 * dx)

    # 三对角矩阵的对角线元素
    main_diag = np.full(Nx, 1 + 2 * alpha, dtype=float)
    off_diag = np.full(Nx - 1, -alpha, dtype=float)

    # 下对角线（对流项）
    lower_diag = np.full(Nx - 1, -beta, dtype=float)
    # 上对角线（对流项）
    upper_diag = np.full(Nx - 1, beta, dtype=float)

    # 合并所有对角线
    a_band = np.zeros((3, Nx))
    a_band[0, 1:] = upper_diag  # 上对角线从第二个元素开始
    a_band[1, :] = main_diag  # 主对角线
    a_band[2, :-1] = lower_diag  # 下对角线到倒数第二个元素

    return a_band


# 构建系数矩阵
a_band = build_coefficients(Nx, dx, dt, D, v)

# 存储解用于可视化
u_history = []
u_history.append(u.copy())

# 时间迭代
for n in range(1, Nt + 1):
    # 右端向量 (基于前一时间步的解)
    b = u.copy()

    # 处理边界条件
    # Dirichlet 边界条件示例：左边界 u=1，右边界 u=0
    b[0] = 1  # u(0) = 1
    b[-1] = 0  # u(L) = 0

    # 构建三对角矩阵系统 Au^{n+1} = b
    # 使用 solve_banded 解决带状矩阵
    # solve_banded 需要带宽参数 (下带宽, 上带宽)，这里均为1
    solution = solve_banded((1, 1), a_band, b)
    u[:] = solution
    u_history.append(u.copy())

# 转换为numpy数组以便后续处理
u_history = np.array(u_history)

# 可视化浓度随时间演化
plt.figure(figsize=(10, 6))
for i in range(0, Nt + 1, max(Nt // 10, 1)):
    plt.plot(x, u_history[i], label=f't={i * dt:.2f}')
plt.xlabel('位置 x')
plt.ylabel('浓度 u')
plt.title('浓度分布随时间的演化（Crank-Nicolson 方法）')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('concentration_evolution.png')
plt.show()


# 稳定性验证：误差分析
# 纯扩散场景（v=0）与解析解对比
def pure_diffusion_analytical(x, t, D):
    sig = np.sqrt(4 * D * t)
    # 避免除以零的情况（当t=0时）
    if t == 0:
        return np.exp(-((x - L / 2) ** 2) / (2 * 0.1 ** 2))
    return np.exp(-((x - L / 2) ** 2) / (2 * D * t)) / (sig * np.sqrt(np.pi))


# 重新初始化
u = u0.copy()
u_history = [u.copy()]

# 设置对流速度为0
v_pure = 0.0
a_band_pure = build_coefficients(Nx, dx, dt, D, v_pure)

for n in range(1, Nt + 1):
    b = u.copy()
    # Dirichlet 边界条件
    b[0] = 1
    b[-1] = 0
    solution = solve_banded((1, 1), a_band_pure, b)
    u[:] = solution
    u_history.append(u.copy())

u_history_pure = np.array(u_history)


# 计算均方误差（MSE）
def mse(u_num, u_exact):
    return np.mean((u_num - u_exact) ** 2)


mse_list = []
times = np.arange(0, T + dt, dt)
for i, t in enumerate(times):
    u_exact = pure_diffusion_analytical(x, t, D)
    mse_val = mse(u_history_pure[i], u_exact)
    mse_list.append(mse_val)

# 绘制误差随时间变化
plt.figure(figsize=(8, 5))
plt.plot(times, mse_list, label='MSE')
plt.xlabel('时间 t')
plt.ylabel('均方误差 (MSE)')
plt.title('纯扩散场景下 Crank-Nicolson 方法的误差分析')
plt.yscale('log')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('error_analysis_pure_diffusion.png')
plt.show()

# CFL 条件验证
# 改变dt，观察稳定性
dt_values = [0.1, 0.01, 0.001]  # 包括不稳定的dt=0.1
colors = ['red', 'blue', 'green']
labels = ['dt=0.1', 'dt=0.01', 'dt=0.001']

plt.figure(figsize=(8, 5))
for dt_val, color, label in zip(dt_values, colors, labels):
    # 重新初始化
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

plt.xlabel('位置 x')
plt.ylabel('浓度 u')
plt.title('不同时间步长下的稳态浓度分布')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('cfl_condition_verification.png')
plt.show()


# 非齐次边界条件影响分析
# Dirichlet 边界条件随时间线性增加：u(0,t) = t, u(L,t)=0
def linear_dirichlet(t):
    return t


# 重新初始化
u = u0.copy()
u_history = [u.copy()]

for n in range(1, Nt + 1):
    t = n * dt
    b = u.copy()
    # Dirichlet 左边界：u= t, 右边界 u=0
    b[0] = linear_dirichlet(t)
    b[-1] = 0
    solution = solve_banded((1, 1), a_band, b)
    u[:] = solution
    u_history.append(u.copy())

u_history = np.array(u_history)

# 绘制浓度分布随时间演化（非齐次边界）
plt.figure(figsize=(10, 6))
for i in range(0, Nt + 1, max(Nt // 10, 1)):
    plt.plot(x, u_history[i], label=f't={i * dt:.2f}')
plt.xlabel('位置 x')
plt.ylabel('浓度 u')
plt.title('非齐次 Dirichlet 边界条件下浓度分布随时间的演化')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig('nonhomogeneous_dirichlet_evolution.png')
plt.show()
#输运方程求解

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve_banded

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
    raise ValueError("稳定性条件不满足！")  # 验证数值稳定性

# 初始条件：高斯脉冲
x = np.linspace(0, space_length, grid_points)
initial_displacement = np.exp(-((x - 0.5)/0.15)**2)  # 初始位移
initial_velocity = np.zeros_like(x)  # 初始速度为0

