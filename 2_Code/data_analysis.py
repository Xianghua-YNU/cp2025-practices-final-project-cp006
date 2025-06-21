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
