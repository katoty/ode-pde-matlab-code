clear;clc;close all;

% 参数设置
L = 1;
T = 1;
h = 0.05;
tau = 0.04;
a = -1;
r = a*tau/h;
t_idx = [0, 0.12 , 0.2 , 0.8];

% 区间数
Nx = L/h ; 
Nt = T/tau;
x = linspace(0,L,Nx+1);
t = linspace(0,T,Nt+1);

% 初始条件
u = zeros(Nx+1,Nt+1);
for i = 1:Nx+1
    u(i,1) = (sin(pi*x(i)))^40;
end

% 求解 u(j,n+1)=(1+r)u(j,n)-ru(j+1,n)
for n = 2:Nt+1
    for i = 1:Nx
        u(i,n) = (1+r)*u(i,n-1)-r*u(i+1,n-1); 
    end
    u(end,n) = u(1,n);
end

% 精确解
u_exact = @(x,t) (sin(pi*(x+t))).^40;
U_exact = zeros(Nx+1, Nt+1);
for n = 1:Nt+1
    for i = 1:Nx+1
        U_exact(i,n) = (sin(pi*(x(i) + t(n))))^40;
    end
end

% 绘图
figure;

% 绘制第一组对比图 (t = 0)
plot(x, u(:, 1), 'b--', 'LineWidth', 2); 
hold on;
plot(x, u_exact(x, 0), 'b-', 'LineWidth', 2); 

% 绘制第二组对比图 (t = 0.12)
plot(x, u(:, 4), 'g--', 'LineWidth', 2); 
plot(x, u_exact(x, 0.12), 'g-', 'LineWidth', 2); 

% 绘制第三组对比图 (t = 0.2)
plot(x, u(:, 6), 'c--', 'LineWidth', 2); 
plot(x, u_exact(x, 0.2), 'c-', 'LineWidth', 2);

% 绘制第四组对比图 (t = 0.8)
plot(x, u(:, 21), 'y--', 'LineWidth', 2); 
plot(x, u_exact(x, 0.8), 'y-.', 'LineWidth', 2); 

legend('t = 0 数值解', 't = 0 精确解', 't = 0.12 数值解', 't = 0.12 精确解', 't = 0.20 数值解', 't = 0.20 精确解', 't = 0.80 数值解', 't = 0.80 精确解');
xlabel('x');
ylabel('u');



