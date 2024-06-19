clear;clc;
% 初始化条件和参数
h1 = 0.1;
h2 = 0.05;

t1 = 0:h1:1;
t2 = 0:h2:1;

N1 = 1/h1;
N2 = 1/h2;

u1 = zeros(1,N1+1);
v1 = zeros(1,N1+1);
u2 = zeros(1,N2+1); 
v2 = zeros(1,N2+1);
u_exact = zeros(1,N2+1);

%% Euler
disp('Euler:');
% h1
u1(1) = 0;
v1(1) = 1;
for n=1:N1
    u1(n+1) = u1(n) + h1*v1(n);
    v1(n+1) = v1(n) - h1*u1(n);
end
result1 = [t1',u1'];
disp('h=0.1的结果')
disp(result1);

% h2
u2(1) = 0;
v2(1) = 1;
for n=1:N2
    u2(n+1) = u2(n) + h2*v2(n);
    v2(n+1) = v2(n) - h2*u2(n);
end
disp('h=0.05的结果')
result2 = [t2',u2'];
disp(result2);

%% improved euler
disp('improved Euler:');
% h1
for n=1:N1
    u1(n+1) = ((1-h1^2/4)*u1(n)+h1*v1(n))/(1+h1^2/4);
    v1(n+1) = ((1-h1^2/4)*v1(n)-h1*u1(n))/(1+h1^2/4);
end
result1 = [t1',u1'];
disp('h=0.1的结果')
disp(result1);
%h2
for n=1:N2
    u2(n+1) = ((1-h2^2/4)*u2(n)+h2*v2(n))/(1+h2^2/4);
    v2(n+1) = ((1-h2^2/4)*v2(n)-h2*u2(n))/(1+h2^2/4);
end
disp('h=0.05的结果')
result2 = [t2',u2'];
disp(result2);

%精确解
for n=1:N2+1
    u_exact(n) = sin(t2(n));
end
disp('精确解的结果')
result3 = [t2',u_exact'];
disp(result3);

%% 绘图
% Euler
figure(1);
plot(t1,u1,"r*");
hold on;
plot(t2,u2,"b^");
plot(t2,u_exact);
hold off;
grid on;
title('Euler法');
legend('h=0.1','h=0.05','精确解');
% Improved Euler
figure(2);
plot(t1,u1,"r*");
hold on;
plot(t2,u2,"b^");
plot(t2,u_exact);
hold off;
grid on;
title('改进Euler法');
legend('h=0.1','h=0.05','精确解');