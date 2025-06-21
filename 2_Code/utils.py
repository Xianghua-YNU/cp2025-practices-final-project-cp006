"""
工具函数模块：包含一些通用的辅助函数。
"""
#芯片热传导问题的数值求解与方法对比研究
import time

def timer_decorator(func):
    """计时装饰器"""
    def wrapper(*args, **kwargs):
        start_time = time.time()
        result = func(*args, **kwargs)
        end_time = time.time()
        print(f"{func.__name__} 耗时: {end_time - start_time:.3f}秒")
        return result
    return wrapper
