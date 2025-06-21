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
    
    # 计算收敛率（最后10%迭代的平均下降率）
    if iterations > 10:
        last_10pct = int(iterations * 0.1)
        rate = np.mean(np.diff(np.log(history[-last_10pct:])))
        print(f"收敛率: {rate:.4f} (指数衰减率)")
    
    return {
        'iterations': iterations,
        'final_error': final_error,
        'convergence_rate': rate if iterations > 10 else None
    }


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
