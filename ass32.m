
%改进euler法解决二阶微分方程
clear all;
clc;
% 定义初始条件和参数
u0 = 0;  % 初始条件 u(0)=0
v0 = 1;  % 初始条件 u'(0)=1
h = 0.05;  % 时间步长
N = round(1/h);  % 

% 初始化数组
t = (0:N)*h; 
u = zeros(1, N + 1);
v = zeros(1, N + 1);
u(1) = u0;
v(1) = v0;

% 改进的 Euler 方法迭代求解
for n = 1:N
    % 预测步
    u_pred = u(n) + h * v(n);
    v_pred = v(n) - h * u(n);
    
    % 校正步
    iterations = 100; % 设置迭代次数
    for i = 1:iterations
        u_new = u(n) + h/2 * (v(n) + v_pred);
        v_new = v(n) - h/2 * (u(n) + u_pred);
        
        % 检查收敛条件
        if abs(u_new - u(n+1)) < 1e-6 && abs(v_new - v(n+1)) < 1e-6
            u(n+1) = u_new;
            v(n+1) = v_new;
            break;
        else
            u(n+1) = u_new;
            v(n+1) = v_new;
        end
    end
end

results = [t', u']; % 将结果保存到 results 中

