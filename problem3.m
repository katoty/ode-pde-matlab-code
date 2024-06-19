clear;clc;close all;
% 网格步长
H = [1/16,1/32,1/64];
for t =1:3
    h = H(t);

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
    % 因高斯赛德尔计算较慢，调试时用\计算。
    % u_vector = gauss_seidel(A,h^2*b,1e-4,1000);
    u_vector = A \ (h^2 * b);

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
    % 输出指定位置的数值解，解析解。
    fprintf("h= %s 时: \n",rats(h));
    fprintf("(1/4,1/4)数值解 %f \n",u(1/(4*h)+1,1/(4*h)+1));
    fprintf("(1/4,1/4)解析解 %f \n",U_exact(1/(4*h)+1,1/(4*h)+1));
    fprintf("(1/2,1/2)数值解 %f \n",u(1/(2*h)+1,1/(2*h)+1));
    fprintf("(1/2,1/2)解析解 %f \n",U_exact(1/(2*h)+1,1/(2*h)+1));
    fprintf("(3/4,3/4)数值解 %f \n",u(3/(4*h)+1,3/(4*h)+1));
    fprintf("(3/4,3/4)解析解 %f \n",U_exact(3/(4*h)+1,3/(4*h)+1));

end

% 绘制数值解图,h=1/32
figure;

% y =1/4 为第9列 = y/h +1(因为第1列为0,下面同理。)
plot(x,u(:,9));
hold on;

% y =1/2 
plot(x,u(:,17));

% y =3/4
plot(x,u(:,25));

legend('y=1/4','y=1/2','y=3/4');
xlabel('x');
ylabel('z');
title("差分解曲线");



[X, Y] = meshgrid(x, y);
figure;
surf(X, Y, u);
title('Gauss\_sidel求解泊松方程的数值解');
xlabel('x');
ylabel('y');
zlabel('u');

figure;
surf(X, Y, U_exact);
title('泊松方程的精确解');
xlabel('x');
ylabel('y');
zlabel('u\_exact');




