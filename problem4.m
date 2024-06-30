clear;clc;

% 函数定义
solve_heat_equation_BD = @(h, tau) heat_equation_solver_backward_difference(h, tau);
solve_heat_equation_CN = @(h, tau) heat_equation_solver_crank_nicolson(h, tau);

% 方案1: h=1/40, τ=1/1600
[x1_BD, u1_BD] = solve_heat_equation_BD(1/40, 1/1600);
[x1_CN, u1_CN] = solve_heat_equation_CN(1/40, 1/1600);

% 方案2: h=1/80, τ=1/3200
[x2_BD, u2_BD] = solve_heat_equation_BD(1/80, 1/3200);
[x2_CN, u2_CN] = solve_heat_equation_CN(1/80, 1/3200);

% 解析解
x_eval = [1/5, 2/5, 3/5, 4/5, 1];
T = 1;
analytical = exp(-pi^2*T) * cos(pi * x_eval) + (1 - cos(T));

% 直接从计算结果中提取所需点的值
indices1 = round(x_eval / (1/40)) + 1;  % 方案1的索引
indices2 = round(x_eval / (1/80)) + 1;  % 方案2的索引

u1_BD_eval = u1_BD(indices1);
u2_BD_eval = u2_BD(indices2);
u1_CN_eval = u1_CN(indices1);
u2_CN_eval = u2_CN(indices2);

% 显示结果
fprintf('时间 t = %.1f 时的结果:\n', T);
fprintf('x\t\t方案1(BD)\t方案1(CN)\t方案2(BD)\t方案2(CN)\t解析解\n');
for i = 1:5
    fprintf('%.1f\t\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n', x_eval(i), u1_BD_eval(i), u1_CN_eval(i), u2_BD_eval(i), u2_CN_eval(i), analytical(i));
end

%% 向后差分
function [x, u] = heat_equation_solver_backward_difference(h, tau)
    L = 1;        
    T = 1;        
    a = 1;        % 热系数

    % 参数设置
    Nx = L / h;  
    Nt = T / tau; 

    
    r = a * tau / h^2;

    % 初始条件
    x = linspace(0, L, Nx+1)'; 
    t = linspace(0, T, Nt+1)';
    u = cos(pi * x); 

    % 系数矩阵A
    A = diag((1 + 2*r) * ones(Nx+1, 1)) + diag(-r * ones(Nx, 1), 1) + diag(-r * ones(Nx, 1), -1);

    % Neumann边界条件的调整
    A(1,2) = -2*r;
    A(end,end-1) = -2*r;

    % 右端函数
    f = @(t) sin(t);

    % 求解
    for n = 2:Nt+1
        b = u + tau * f(t(n));
        u = A \ b;
    end
end

%% CN法
function [x, u] = heat_equation_solver_crank_nicolson(h, tau)
    L = 1;        
    T = 1;        
    a = 1;        % 热系数

    % 参数设置
    Nx = L / h;   
    Nt = T / tau;
    r = a * tau / h^2;

    % 初始条件
    x = linspace(0, L, Nx+1)'; 
    t = linspace(0, T, Nt+1)';
    u = cos(pi * x);

    % 系数矩阵
    A = zeros(Nx+1, Nx+1);
    B = zeros(Nx+1, Nx+1);
    A(1,1:2) = [1+r, -r];
    B(1,1:2) = [1-r, r];
    for i = 2:Nx
        A(i,i-1:i+1) = [-r/2, 1+r, -r/2];
        B(i,i-1:i+1) = [r/2, 1-r, r/2];
    end
    A(Nx+1,Nx:Nx+1) = [-r, 1+r];
    B(Nx+1,Nx:Nx+1) = [r, 1-r];

    % 右端函数
    f = @(t) sin(t);

    % 求解
    for n = 2:Nt+1
        b = B * u + 0.5 * tau * (f(t(n)) + f(t(n-1)));
        u = A \ b;
    end
end
