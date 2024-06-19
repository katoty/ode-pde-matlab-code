function x = gauss_seidel(A, b, tol, max_iter)
    % 高斯赛德尔迭代法
    % A: 系数矩阵
    % b: 常数项向量
    % tol: 误差容限
    % max_iter: 最大迭代次数
    
    % 初始化
    n = length(b);
    x = zeros(n, 1);  % 初始解向量
    x_old = x;        % 前一次迭代的解向量
    
    for iter = 1:max_iter
        for i = 1:n
            % 计算当前行的新的 x 值
            sum1 = A(i, 1:i-1) * x(1:i-1);
            sum2 = A(i, i+1:n) * x_old(i+1:n);
            x(i) = (b(i) - sum1 - sum2) / A(i, i);
        end
        
        % 检查收敛性
        if norm(x - x_old, inf) < tol
            fprintf('迭代收敛于第 %d 次迭代。\n', iter);
            return;
        end
        
        % 更新旧的解向量
        x_old = x;
    end
    
    warning('达到最大迭代次数，未收敛。');
end