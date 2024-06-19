clear;clc;

% 初始化条件和参数
h1 = 0.2;
h2 = 0.4;
h3 = 0.5;

t1 = 0:h1:2;
t2 = 0:h2:2;
t3 = 0:h3:2;

[N1, N2, N3] = deal(2/h1, 2/h2, 2/h3);
[u1, u2, u3, u_exact] = deal(zeros(1, N1+1), zeros(1, N2+1), zeros(1, N3+1), zeros(1, N1+1));

[u1(1), u2(1), u3(1), u_exact(1)] = deal(1);

% 定义函数 f(t, u)
f = @(t, u) 4 * t * sqrt(u);

% h1
disp('h=0.2的结果')
for n = 1:N1
    k1 = f(t1(n), u1(n));
    k2 = f(t1(n) + 0.5 * h1, u1(n) + 0.5 * h1 * k1);
    k3 = f(t1(n) + 0.5 * h1, u1(n) + 0.5 * h1 * k2);
    k4 = f(t1(n) + h1, u1(n) + h1 * k3);
    u1(n+1) = u1(n) + h1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
end
disp([t1', u1']);

% h2
disp('h=0.4的结果')
for n = 1:N2
    k1 = f(t2(n), u2(n));
    k2 = f(t2(n) + 0.5 * h2, u2(n) + 0.5 * h2 * k1);
    k3 = f(t2(n) + 0.5 * h2, u2(n) + 0.5 * h2 * k2);
    k4 = f(t2(n) + h2, u2(n) + h2 * k3);
    u2(n+1) = u2(n) + h2 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
end
disp([t2', u2']);

% h3
disp('h=0.5的结果')
for n = 1:N3
    k1 = f(t3(n), u3(n));
    k2 = f(t3(n) + 0.5 * h3, u3(n) + 0.5 * h3 * k1);
    k3 = f(t3(n) + 0.5 * h3, u3(n) + 0.5 * h3 * k2);
    k4 = f(t3(n) + h3, u3(n) + h3 * k3);
    u3(n+1) = u3(n) + h3 / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
end
disp([t3', u3']);

% 精确解
disp('精确解的结果')
for n = 1:N1+1
    u_exact(n) = (1 + t1(n)^2)^2;
end
disp([t1', u_exact']);
