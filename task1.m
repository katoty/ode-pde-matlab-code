
% 使用改进的 Euler 方法一阶微分方程
clear all;
clc;
% 定义微分方程 dy/dx = f(x, y)，这里选择的是简单例子 y' = -5y
f = @(x, y) -5 * y;

% 初始条件
x0 = 0;
y0 = 1;

% 区间和步长
a = 0;
b = 1;
h = 0.05;

% 计算步数
N = (b - a) / h;

% 初始化结果数组
x = zeros(1, N+1);
y = zeros(1, N+1);

% 将初始值存入结果数组
x(1) = x0;
y(1) = y0;

% 设定 epsilon 值
epsilon = 0.0001; % 可以根据需要调整

% 使用改进的 Euler 方法进行迭代
for i = 1:N
    x(i+1) = x(i) + h; % 更新 x
    
    % 迭代
    iterations = 100; % 每个步骤上进行多次迭代
    
    temp=y(i);
  
    % 在每个步骤上进行多次迭代
    for j = 1:iterations
        % 计算 k1 和 k2
        k1 = h * f(x(i), y(i));
        k2 = h * f(x(i+1), temp);
        
        % 计算改进的 Euler 迭代结果
        y_new = y(i) + 0.5 * (k1 + k2); 
        
        % 检查当前迭代结果与上一次的差异
        if abs(y_new - temp) < epsilon
            y(i+1) = y_new; % 当差异低于阈值时，接受当前迭代结果
            break; % 停止迭代
        else
            temp = y_new; % 否则，继续迭代
        end
    end
end

% 绘制结果
plot(x, y);
xlabel('x');
ylabel('y');
title('改进的 Euler 方法求解常微分方程 y'' = -5y');

% 创建结果矩阵
results = [x', y'];

% % 将结果保存为 Excel 文件
% filename = 'C:\Users\lty\Desktop\task1_改进euler.xlsx';
% xlswrite(filename, results);
% disp(['Results saved to ' filename]);


