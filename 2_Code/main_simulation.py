"""
主程序入口：用于运行主要的物理模拟。
"""
#芯片热传导问题的数值求解与方法对比研究
import time
from numerical_methods import solve_heat_jacobi, solve_heat_sor
from visualization import plot_temperature

def main():
    # 模拟参数
    N = 64        # 网格数
    L = 1.0       # 芯片边长(m)
    tol = 1e-4    # 收敛容差

    print("求解芯片热传导问题...")
    print(f"网格尺寸: {N}×{N}，芯片边长: {L}m")
    print(f"收敛容差: {tol}")

    # 使用Jacobi方法求解
    print("\n1. Jacobi迭代法:")
    start_time = time.time()
    T_jacobi, iter_jacobi, conv_jacobi, x, y, dx = solve_heat_jacobi(N, L, tol)
    time_jacobi = time.time() - start_time
    print(f"   收敛于 {iter_jacobi} 次迭代")
    print(f"   耗时: {time_jacobi:.3f} 秒")

    # 使用SOR方法求解
    print("\n2. SOR迭代法 (ω=1.8):")
    start_time = time.time()
    T_sor, iter_sor, conv_sor, x, y, dx = solve_heat_sor(N, L, 1.8, tol)
    time_sor = time.time() - start_time
    print(f"   收敛于 {iter_sor} 次迭代")
    print(f"   耗时: {time_sor:.3f} 秒")

    # 性能对比
    print("\n3. 方法性能对比:")
    print(f"   Jacobi: {iter_jacobi} 次迭代, {time_jacobi:.3f}s")
    print(f"   SOR:    {iter_sor} 次迭代, {time_sor:.3f}s")
    print(f"   加速比: 迭代次数 {iter_jacobi / iter_sor:.1f}倍, 时间 {time_jacobi / time_sor:.2f}倍")

    # 绘制温度分布
    plot_temperature(T_jacobi, x, y, dx, "Jacobi方法")
    plot_temperature(T_sor, x, y, dx, "SOR方法")

if __name__ == "__main__":
    main()

#波动方程
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# 计算前两层解
u_prev_prev = initial_displacement.copy()
u_prev = u_prev_prev + 0.5 * stability_factor**2 * (
    np.roll(u_prev_prev, -1) - 2*u_prev_prev + np.roll(u_prev_prev, 1)
)
u_prev[0], u_prev[-1] = 0,  # 固定边界条件

# 预计算所有时间层（主循环）
u = np.zeros((time_steps, grid_points))
u[0] = u_prev_prev
u[1] = u_prev

for n in range(1, time_steps-1):
    u_current = 2 * u[n] - u[n-1] + stability_factor**2 * (
        np.roll(u[n], -1, axis=-1) - 2*u[n] + np.roll(u[n], 1, axis=-1)
    )
    u_current[[0, -1]] = 0  # 边界条件
    u[n+1] = u_current

# 调用可视化（动画部分）
# （注：主函数最后一步触发可视化，代码见下方“可视化”部分）


2.输运方程隐式格式求解、边界条件处理及稳定性验证研究

import numpy as np
from numerical_methods import solve_heat_equation, build_coefficients
from visualization import plot_concentration_evolution, plot_error_analysis
from visualization import plot_cfl_verification, plot_nonhomogeneous_boundary
from utils import setup_font

# 主模拟程序：协调参数设置、数值求解与可视化
if __name__ == "__main__":
    # 初始化字体设置
    setup_font()
    
    # 参数设置（添加单位注释）
    L = 1.0  # 空间长度（m）
    Nx = 50  # 空间离散点数
    dx = L / (Nx - 1)  # 空间步长（m）
    dt = 0.01  # 时间步长（s）
    T = 1.0  # 总时间（s）
    Nt = int(T / dt)  # 时间步数
    D = 0.01  # 扩散系数（m²/s）
    v = 0.1  # 对流速度（m/s）
    
    # 输出模拟参数
    print("===== 模拟参数 =====")
    print(f"空间长度: L = {L} m")
    print(f"离散点数: Nx = {Nx}")
    print(f"空间步长: dx = {dx:.4f} m")
    print(f"时间步长: dt = {dt} s")
    print(f"总时间: T = {T} s")
    print(f"扩散系数: D = {D} m²/s")
    print(f"对流速度: v = {v} m/s")
    print("====================\n")
    
    # 初始条件（高斯分布）
    x = np.linspace(0, L, Nx)
    u0 = np.exp(-((x - L / 2) **2) / (2 * 0.1** 2))  # 初始浓度分布（kg/m³）
    
    # 执行热传导方程求解
    u_history = solve_heat_equation(Nx, dx, dt, T, D, v, u0)
    
    # 执行可视化分析
    plot_concentration_evolution(x, u_history, Nt, dt)
    plot_error_analysis(x, u0, Nx, dx, dt, T, D)
    plot_cfl_verification(x, u0, Nx, dx, T, D, v)
    plot_nonhomogeneous_boundary(x, u0, Nx, dx, dt, T, D, v)
