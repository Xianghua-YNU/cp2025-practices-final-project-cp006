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

2.输运方程隐式格式求解、边界条件处理及稳定性验证研究


import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

def get_available_chinese_font():
    """检测系统中可用的中文字体"""
    available_fonts = {f.name for f in fm.fontManager.ttflist}
    chinese_fonts = [
        'SimHei',       # Windows常用黑体
        'Microsoft YaHei',  # Windows雅黑
        'SimSun',       # 宋体
        'KaiTi',        # 楷体
        'WenQuanYi Micro Hei',  # Linux常用
        'Heiti TC'       # macOS黑体
    ]
    for font in chinese_fonts:
        if font in available_fonts:
            return font
    return None

def setup_font():
    """设置 matplotlib 字体"""
    available_font = get_available_chinese_font()
    if available_font:
        plt.rcParams["font.family"] = available_font
        print(f"使用可用字体: {available_font}")
    else:
        plt.rcParams["font.family"] = "sans-serif"
        print("警告: 未找到中文字体，使用默认字体")
    plt.rcParams["axes.unicode_minus"] = False  # 正确显示负号
