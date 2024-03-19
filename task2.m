

%euler法求解二阶微分方程
clear all;
clc;
% 定义初始条件和参数
u0 = 0;  % 初始条件 u(0)=0
v0 = 1;  % 初始条件 u'(0)=1
h = 0.05;  % 时间步长
N = 1/h;  % 

% 初始化数组
t=0:0.05:1
u = zeros(1, N + 1);
v = zeros(1, N + 1);
u(1) = u0;
v(1) = v0;

% Euler 方法迭代求解
for n = 1:N
    u(n+1) = u(n) + h * v(n);
    v(n+1) = v(n) - h * u(n);
end

results=[t',u'];

% 将结果保存为 Excel 文件
filename = 'xxx'; %xxx写保存路径
xlswrite(filename, results);
disp(['Results saved to ' filename]);