clear;clc;
% 网格步长
h = 1/64;

% 计算节点数
N = 1/h + 1;

% 初始化参数
x = linspace(0, 1, N);  % x方向的坐标
y = linspace(0, 1, N);  % y方向的坐标

% 右端项的函数定义
f = @(x,y) 2*pi^2*exp(pi*(x+y))*(sin(pi*x)*cos(pi*y)+cos(pi*x)*sin(pi*y));

% 精确解函数定义
u_exact = @(x,y) exp(pi*(x+y))*sin(pi*x)*sin(pi*y);

% 内部网格点数
n = N - 2;  % 不包括边界点

% 初始化系数矩阵和右端项向量
A = zeros(n*n, n*n);
b = zeros(n*n, 1);

% 填充右端项向量
for i = 1:n
    for j = 1:n
        K = (j-1)*n + i;
        b(K) = f(x(i+1), y(j+1));
    end
end

% 填充系数矩阵
for i = 1:n
    for j = 1:n
        K = (j-1)*n + i;
        
        A(K, K) = -4;  % 中心点
        if i > 1
            A(K, K-1) = 1;  % 左边点
        end
        if i < n
            A(K, K+1) = 1;  % 右边点
        end
        if j > 1
            A(K, K-n) = 1;  % 下边点
        end
        if j < n
            A(K, K+n) = 1;  % 上边点
        end
    end
end

u_vector = gauss_seidel(A,b,1e-4,100);
% u_vector = A \ (h^2 * b);

% 将解向量转换回网格
u = zeros(N, N);
for i = 1:n
    for j = 1:n
        K = (j-1)*n + i;
        u(i+1, j+1) = u_vector(K);
    end
end

% 计算精确解在节点处的值
U_exact = zeros(N, N);
for i = 1:N
    for j = 1:N
        U_exact(i,j) = u_exact(x(i), y(j));
    end
end

% 绘制数值解图
[X, Y] = meshgrid(x, y);
figure;
surf(X, Y, u);
title('系数矩阵法求解泊松方程的数值解');
xlabel('x');
ylabel('y');
zlabel('u');

% 绘制精确解图
figure;
surf(X, Y, U_exact);
title('泊松方程的精确解');
xlabel('x');
ylabel('y');
zlabel('u\_exact');




