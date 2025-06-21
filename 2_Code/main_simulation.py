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
