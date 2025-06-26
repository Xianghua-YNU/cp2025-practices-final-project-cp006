# 2_code/ -源代码目录

**目的:** 存放所有用于模拟、分析和可视化的代码。代码的质量是评分的重要组成部分。

**内容:**
README.md`（代码说明书）
# 芯片热传导数值模拟项目

## 项目结构
- `main_simulation.py`: 主程序，控制模拟流程
- `numerical_methods.py`: 数值求解算法（Jacobi和SOR迭代法）
- `visualization.py`: 可视化绘图函数
- `data_analysis.py`: 数据分析工具
- `utils.py`: 辅助工具函数
- `requirements.txt`: 项目依赖

## 环境配置
1. 安装Python 3.7+
2. 安装依赖库：
```bash
pip install -r requirements.txt
```

## 运行过程
运行过程说明
参数输出：
程序开始后会显示模拟参数（网格尺寸、芯片边长、收敛容差等）。
求解过程：
先运行 Jacobi 迭代法，显示迭代次数和耗时；
再运行 SOR 迭代法，同样显示迭代次数和耗时；
最后输出两种方法的性能对比（加速比）。
可视化窗口：
自动弹出两个图形窗口，分别展示 Jacobi 和 SOR 方法的温度分布：
左侧为 3D 温度曲面图（含温度棒）；
右侧为 2D 等值线与热流场图（含温度棒）。

## 功能说明
1. **主程序**：设置模拟参数，调用数值方法求解热传导方程，分析结果并可视化
2. **数值方法**：实现了Jacobi和SOR迭代法求解非齐次边界条件下的热传导方程
3. **可视化**：生成3D温度分布图和2D温度等值线+热流场图
4. **数据分析**：计算误差和收敛率等指标
5. **辅助工具**：提供计时装饰器等实用功能


### 模块化优势
1. **代码复用**：`numerical_methods.py`中的求解器可用于其他热传导问题
2. **可维护性**：功能分离，便于修改和扩展
3. **协作开发**：不同模块可由不同人员并行开发
4. **测试方便**：各模块可独立测试
### 芯片热传导项目完整代码预览（方便读者运行复现）
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import time
from matplotlib import cm  

plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['mathtext.fontset'] = 'cm'


def solve_heat_jacobi(N, L=1.0, tol=1e-5):
    dx = L / (N - 1)
    x = np.linspace(0, L, N)
    y = np.linspace(0, L, N)
    X, Y = np.meshgrid(x, y)

    T = np.zeros((N, N)) + 300  

    iterations = 0
    max_iter = 10000
    convergence_history = []

    while iterations < max_iter:
        T_old = T.copy()

        for i in range(1, N - 1):
            for j in range(1, N - 1):
                q = 100 * np.exp(-((X[i, j] - 0.5) ** 2 + (Y[i, j] - 0.7) ** 2) / (2 * 0.1 ** 2))
                T[i, j] = 0.25 * (T[i + 1, j] + T[i - 1, j] + T[i, j + 1] + T[i, j - 1] + dx ** 2 * q)

        T[0, :] = 300
        T[:, -1] = 300 + 20 * y
        T[1, :] = T[2, :] - 10 * dx * np.sin(2 * np.pi * y)
        T[-2, :] = (T[-3, :] + 10 * dx * 310) / (1 + 5 * dx)

        max_change = np.max(np.abs(T - T_old))
        convergence_history.append(max_change)
        iterations += 1

        if max_change < tol:
            break

    return T, iterations, convergence_history, x, y, dx  


def solve_heat_sor(N, L=1.0, omega=1.8, tol=1e-5):
    dx = L / (N - 1)
    x = np.linspace(0, L, N)
    y = np.linspace(0, L, N)
    X, Y = np.meshgrid(x, y)

    T = np.zeros((N, N)) + 300  

    iterations = 0
    max_iter = 10000
    convergence_history = []

    while iterations < max_iter:
        T_old = T.copy()

        for i in range(1, N - 1):
            for j in range(1, N - 1):
                q = 100 * np.exp(-((X[i, j] - 0.5) ** 2 + (Y[i, j] - 0.7) ** 2) / (2 * 0.1 ** 2))
                r_ij = 0.25 * (T[i + 1, j] + T[i - 1, j] + T[i, j + 1] + T[i, j - 1] + dx ** 2 * q)
                T[i, j] = (1 - omega) * T[i, j] + omega * r_ij

        T[0, :] = 300
        T[:, -1] = 300 + 20 * y
        T[1, :] = T[2, :] - 10 * dx * np.sin(2 * np.pi * y)
        T[-2, :] = (T[-3, :] + 10 * dx * 310) / (1 + 5 * dx)

        max_change = np.max(np.abs(T - T_old))
        convergence_history.append(max_change)
        iterations += 1

        if max_change < tol:
            break

    return T, iterations, convergence_history, x, y, dx  


def plot_temperature(T, x, y, dx, method_name):
    fig = plt.figure(figsize=(12, 5))

    ax1 = fig.add_subplot(121, projection='3d')
    X, Y = np.meshgrid(x, y)
    surf = ax1.plot_surface(X, Y, T, cmap=cm.jet, alpha=0.8)
    ax1.set_xlabel('X (m)')
    ax1.set_ylabel('Y (m)')
    ax1.set_zlabel('温度 (K)')
    ax1.set_title(f'芯片温度分布\n({method_name})')


    ax2 = fig.add_subplot(122)
    X, Y = np.meshgrid(x, y)
    Ty, Tx = np.gradient(-T, dx)  
    levels = np.linspace(T.min(), T.max(), 15)
    contour = ax2.contour(X, Y, T, levels=levels, colors='red', linestyles='solid')
    ax2.clabel(contour, inline=True, fontsize=6)
    ax2.streamplot(X, Y, Tx, Ty, density=1.5, color='blue', linewidth=1, arrowsize=1.5)
    heat_source = 100 * np.exp(-((X - 0.5) ** 2 + (Y - 0.7) ** 2) / (2 * 0.1 ** 2))
    ax2.contour(X, Y, heat_source, levels=[5], colors='black', linestyles='dashed')
    ax2.text(0.5, 0.7, '热源', fontsize=9, ha='center')

    ax2.set_xlabel('X (m)')
    ax2.set_ylabel('Y (m)')
    ax2.set_title(f'温度等值线与热流场\n({method_name})')
    ax2.set_aspect('equal')

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    N = 64  
    L = 1.0  
    tol = 1e-4  

    print("求解芯片热传导问题...")
    print(f"网格尺寸: {N}×{N}，芯片边长: {L}m")
    print(f"收敛容差: {tol}")

    print("\n1. Jacobi迭代法:")
    start_time = time.time()
    T_jacobi, iter_jacobi, conv_jacobi, x, y, dx_jacobi = solve_heat_jacobi(N, L, tol)
    time_jacobi = time.time() - start_time
    print(f"   收敛于 {iter_jacobi} 次迭代")
    print(f"   耗时: {time_jacobi:.3f} 秒")

    print("\n2. SOR迭代法 (ω=1.8):")
    start_time = time.time()
    T_sor, iter_sor, conv_sor, x, y, dx_sor = solve_heat_sor(N, L, 1.8, tol)
    time_sor = time.time() - start_time
    print(f"   收敛于 {iter_sor} 次迭代")
    print(f"   耗时: {time_sor:.3f} 秒")

    print("\n3. 方法性能对比:")
    print(f"   Jacobi: {iter_jacobi} 次迭代, {time_jacobi:.3f}s")
    print(f"   SOR:    {iter_sor} 次迭代, {time_sor:.3f}s")
    print(f"   加速比: 迭代次数 {iter_jacobi / iter_sor:.1f}倍, 时间 {time_jacobi / time_sor:.2f}倍")

    plot_temperature(T_jacobi, x, y, dx_jacobi, "Jacobi方法")
    plot_temperature(T_sor, x, y, dx_sor, "SOR方法")


    plt.tight_layout()
    plt.show()

