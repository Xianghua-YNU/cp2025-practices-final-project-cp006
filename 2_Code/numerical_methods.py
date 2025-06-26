"""
核心数值算法模块：包含各种数值方法实现。
"""
#芯片热传导问题的数值求解与方法对比研究
import numpy as np

def solve_heat_jacobi(N, L=1.0, tol=1e-5):
    """使用Jacobi迭代法求解芯片热传导问题"""
    dx = L / (N - 1)
    x = np.linspace(0, L, N)
    y = np.linspace(0, L, N)
    X, Y = np.meshgrid(x, y)

    # 初始化温度场
    T = np.zeros((N, N)) + 300  # 初始温度设为底部恒温值

    iterations = 0
    max_iter = 10000
    convergence_history = []

    while iterations < max_iter:
        T_old = T.copy()

        # Jacobi迭代（内部点）
        for i in range(1, N - 1):
            for j in range(1, N - 1):
                # 计算热源
                q = 100 * np.exp(-((X[i, j] - 0.5) ** 2 + (Y[i, j] - 0.7) ** 2) / (2 * 0.1 ** 2))
                # 热传导方程离散形式
                T[i, j] = 0.25 * (T[i + 1, j] + T[i - 1, j] + T[i, j + 1] + T[i, j - 1] + dx ** 2 * q)

        # 处理非齐次边界条件
        T[0, :] = 300                          # 底部：恒温边界 T=300K
        T[:, -1] = 300 + 20 * y                # 右侧：T(1,y)=300+20y
        T[1, :] = T[2, :] - 10 * dx * np.sin(2 * np.pi * y)  # 左侧：热流边界
        T[-2, :] = (T[-3, :] + 10 * dx * 310) / (1 + 5 * dx)  # 顶部：对流边界

        # 计算收敛性
        max_change = np.max(np.abs(T - T_old))
        convergence_history.append(max_change)
        iterations += 1

        if max_change < tol:
            break

    return T, iterations, convergence_history, x, y, dx

def solve_heat_sor(N, L=1.0, omega=1.8, tol=1e-5):
    """使用SOR迭代法求解芯片热传导问题"""
    dx = L / (N - 1)
    x = np.linspace(0, L, N)
    y = np.linspace(0, L, N)
    X, Y = np.meshgrid(x, y)

    # 初始化温度场
    T = np.zeros((N, N)) + 300  # 初始温度设为底部恒温值

    iterations = 0
    max_iter = 10000
    convergence_history = []

    while iterations < max_iter:
        T_old = T.copy()

        # SOR迭代（内部点）
        for i in range(1, N - 1):
            for j in range(1, N - 1):
                # 计算热源
                q = 100 * np.exp(-((X[i, j] - 0.5) ** 2 + (Y[i, j] - 0.7) ** 2) / (2 * 0.1 ** 2))
                # SOR迭代公式
                r_ij = 0.25 * (T[i + 1, j] + T[i - 1, j] + T[i, j + 1] + T[i, j - 1] + dx ** 2 * q)
                T[i, j] = (1 - omega) * T[i, j] + omega * r_ij

        # 处理非齐次边界条件
        T[0, :] = 300                          # 底部：恒温边界 T=300K
        T[:, -1] = 300 + 20 * y                # 右侧：T(1,y)=300+20y
        T[1, :] = T[2, :] - 10 * dx * np.sin(2 * np.pi * y)  # 左侧：热流边界
        T[-2, :] = (T[-3, :] + 10 * dx * 310) / (1 + 5 * dx)  # 顶部：对流边界

        # 计算收敛性
        max_change = np.max(np.abs(T - T_old))
        convergence_history.append(max_change)
        iterations += 1

        if max_change < tol:
            break

    return T, iterations, convergence_history, x, y, dx

输运方程隐式格式求解、边界条件处理及稳定性验证研究

import numpy as np
from scipy.linalg import solve_banded

def build_coefficients(Nx, dx, dt, D, v):
    """构建热传导方程的三对角系数矩阵"""
    alpha = D * dt / dx** 2  # 扩散项系数
    beta = v * dt / (2 * dx)  # 对流项系数
    
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

def solve_heat_equation(Nx, dx, dt, T, D, v, u0):
    """使用Crank-Nicolson方法求解热传导方程"""
    Nt = int(T / dt)
    u = u0.copy()
    a_band = build_coefficients(Nx, dx, dt, D, v)
    u_history = [u.copy()]
    
    # 时间迭代
    for n in range(1, Nt + 1):
        # 右端向量 (基于前一时间步的解)
        b = u.copy()
        
        # 处理边界条件：Dirichlet边界（左边界u=1，右边界u=0）
        b[0] = 1  # u(0) = 1 kg/m³
        b[-1] = 0  # u(L) = 0 kg/m³
        
        # 求解三对角矩阵系统
        solution = solve_banded((1, 1), a_band, b)
        u[:] = solution
        u_history.append(u.copy())
    
    return np.array(u_history)

#波动方程

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
