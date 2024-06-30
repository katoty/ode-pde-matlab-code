clear;clc;close all;

% 参数设置
L = 1; 
T = 1;

h = 1/40; 
tau = 1/1600;

Nx = 1/h;
Ny = 1/h;
Nt = 1/tau;

x = linspace(0,L,Nx+1);
y = linspace(0,L,Ny+1);
t = linspace(0,T,Nt+1); 
r = tau/h^2;

% 结果矩阵
U = zeros(Nx+1,Ny+1,Nt+1);

% 初始条件
for i = 1:Nx+1
    for j = 1:Ny+1
        U(i,j,1) = sin(pi*x(i))*cos(pi*y(j));
    end
end

% 边界条件
U(1,:,:) = 0; % 左边界
U(Nx+1,:,:) = 0; % 右边界

% 系数矩阵1 (x方向)
A1 = zeros(Nx+1,Nx+1);
A1(1,1) = 1;
A1(Nx+1,Nx+1) = 1;
for i = 2:Nx
    A1(i,i-1:i+1) = [-r/32, 1+r/16, -r/32];
end

% 系数矩阵2 (y方向)
A2 = zeros(Ny+1,Ny+1);
A2(1,1:2) = [1+r/16, -r/16];
A2(Ny+1,Ny:Ny+1) = [-r/16, 1+r/16];
for i = 2:Ny
    A2(i,i-1:i+1) = [-r/32, 1+r/16, -r/32];
end

for n = 1:Nt
    % 时间对偶点n+1/2处值
    U_half = zeros(Nx+1,Ny+1);
    
    for j = 2:Ny
        b = zeros(Nx+1,1);
        for i = 2:Nx
            b(i) = (r/32)*(U(i,j+1,n) + U(i,j-1,n)) + (1-r/16)*U(i,j,n);
        end
        U_half(:,j) = A1\b;
    end
    
    % 处理y方向的边界条件
    U_half(:,1) = U_half(:,2);
    U_half(:,Ny+1) = U_half(:,Ny);
    
    % 计算n+1时刻
    for i = 2:Nx
        b = zeros(Ny+1,1);
        for j = 1:Ny+1
            b(j) = (r/32)*(U_half(i+1,j) + U_half(i-1,j)) + (1-r/16)*U_half(i,j);
        end
        U(i,:,n+1) = (A2\b)';
    end
end

% 对比图
figure;

t1 = 0.3/tau + 1;
t2 = 0.5/tau + 1;
t3 = 0.8/tau + 1;

subplot(3,2,1); mesh(U(:,:,t1)); title('数值解 t=0.3');
subplot(3,2,2); mesh(sin(pi*x)' * cos(pi*y) * exp(-pi^2*0.3/8)); title('精确解 t=0.3');

subplot(3,2,3); mesh(U(:,:,t2)); title('数值解 t=0.5');
subplot(3,2,4); mesh(sin(pi*x)' * cos(pi*y) * exp(-pi^2*0.5/8)); title('精确解 t=0.5');

subplot(3,2,5); mesh(U(:,:,t3)); title('数值解 t=0.8');
subplot(3,2,6); mesh(sin(pi*x)' * cos(pi*y) * exp(-pi^2*0.8/8)); title('精确解 t=0.8');

% 数值结果
disp('t = 1 时刻的数值结果：')
disp('  x   y   数值解')

for j = 1:3
    for k = 1:3
        x_idx = j*10 + 1;
        y_idx = k*10 + 1;
        val = U(x_idx, y_idx, end);
        fprintf('%.2f  %.2f  %.6f\n', j/4, k/4, val)
    end
end




