clear;clc;close all;

% 网格步长
H = [1/16,1/32,1/64];
for t =1:3
    h = H(t);

    % 计算节点数
    N = 1/h + 1;

    % 初始化参数
    x = linspace(0, 1, N); 
    y = linspace(0, 1, N);  

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
    g_s = @(A, b, tol, max_iter)  gauss_seidel(A, b, tol, max_iter);

    % 高斯赛德尔计算较慢，调试时用\计算。
    %u_vector = g_s(A,h^2*b,1e-4,1000);
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

    % 输出指定位置的数值解，解析解。设(x，y)位于矩阵（i,j)处，则
    % i = x/h +1 （+1因为第一列为0）。纵坐标同理
    fprintf("h= %s 时: \n",rats(h));
    fprintf("(1/4,1/4)数值解：%f",u(1/(4*h)+1,1/(4*h)+1));
    fprintf(" 解析解： %f \n",U_exact(1/(4*h)+1,1/(4*h)+1));
    fprintf("(1/2,1/2)数值解：%f",u(1/(2*h)+1,1/(2*h)+1));
    fprintf(" 解析解： %f \n",U_exact(1/(2*h)+1,1/(2*h)+1));
    fprintf("(3/4,3/4)数值解：%f",u(3/(4*h)+1,3/(4*h)+1));
    fprintf(" 解析解： %f \n",U_exact(3/(4*h)+1,3/(4*h)+1));

    % 绘制数值解图,h=1/32
    if h==1/32
    
        figure;

        % y =1/4
        plot(x,u(:,9));
        hold on;

        % y =1/2
        plot(x,u(:,17));

        % y =3/4
        plot(x,u(:,25));

        legend('y=1/4','y=1/2','y=3/4');
        xlabel('x');
        ylabel('z');
        title("h=1/32时的差分解曲线");
    end

end

%% 比较数值解 精确解
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

function [x, iter, error] = gauss_seidel(A, b, tol, max_iter)
    % 高斯赛德尔迭代法
    % A: 系数矩阵
    % b: 常数项向量
    % tol: 误差容限
    % max_iter: 最大迭代次数
    % 返回值：
    % x: 解向量
    % iter: 实际迭代次数
    % error: 最终误差
    
    % 输入检查
    if nargin < 4
        error('需要四个输入参数：A, b, tol, max_iter');
    end
    
    [m, n] = size(A);
    if m ~= n || length(b) ~= n
        error('A必须是方阵，且A和b的维度必须匹配');
    end
    
    % 初始化
    x = zeros(n, 1);  % 初始解向量
    x_old = x;        % 前一次迭代的解向量
    inv_diag = 1 ./ diag(A);  % 预先计算对角线元素的倒数
    
    for iter = 1:max_iter
        for i = 1:n
            % 计算当前行的新的 x 值
            sum1 = A(i, 1:i-1) * x(1:i-1);
            sum2 = A(i, i+1:n) * x_old(i+1:n);
            x(i) = (b(i) - sum1 - sum2) * inv_diag(i);
        end
        
        % 检查收敛性
        error = norm(x - x_old, inf);
        if error < tol
            fprintf('迭代收敛于第 %d 次迭代。\n', iter);
            return;
        end
        
        % 更新旧的解向量
        x_old = x;
    end
    
    warning('达到最大迭代次数 %d，未收敛。最终误差: %e', max_iter, error);
end



