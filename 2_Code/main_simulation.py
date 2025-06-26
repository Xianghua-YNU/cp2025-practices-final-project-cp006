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
