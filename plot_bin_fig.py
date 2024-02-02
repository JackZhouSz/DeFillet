import matplotlib.pyplot as plt
import numpy as np

# 从文件中读取所有行的数据
def read_all_data(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
        data = [float(line.strip()) for line in lines]
    return data

# 绘制柱状分布图
def plot_histogram(data):
    plt.hist(data, bins=20, color='blue', edgecolor='black', alpha=0.7)
    plt.title('Histogram of Data')
    plt.xlabel('Values')
    plt.ylabel('Frequency')
    plt.show()

# 文件路径
file_path = './output.txt'  # 替换为实际的文件路径

# 读取所有数据
data = read_all_data(file_path)

# 绘制柱状分布图
plot_histogram(data)
