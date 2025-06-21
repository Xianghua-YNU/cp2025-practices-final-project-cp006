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
