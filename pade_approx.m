function [val] = pade_approx(coeffs)
    % pade_approx: 使用 Wynn's epsilon 算法计算幂级数在 s=1 处的 Padé 逼近
    % 输入 coeffs: 某个节点的电压级数系数向量 [c0, c1, c2, ..., ck]
    
    N = length(coeffs);
    % 计算部分和序列 S_n = sum_{i=0}^n c_i
    S = cumsum(coeffs);
    
    % Wynn's epsilon 表格初始化
    % ep(i, j) 对应 epsilon_{j-2, i}
    % 第一列 (j=1) 是 epsilon_{-1, n} = 0
    % 第二列 (j=2) 是 epsilon_{0, n} = S_n
    ep = zeros(N, N + 1);
    ep(:, 2) = S(:);
    
    % 递归填充表格
    % 递推公式: ep(i, j) = ep(i+1, j-2) + 1 / (ep(i+1, j-1) - ep(i, j-1))
    for j = 3:N+1
        for i = 1:(N - (j - 2) - 1)
            diff = ep(i+1, j-1) - ep(i, j-1);
            
            if abs(diff) < 1e-18 % 防止除零
                ep(i, j) = ep(i+1, j-2); 
            else
                ep(i, j) = ep(i+1, j-2) + 1 / diff;
            end
        end
    end
    
    % 结果通常取偶数列 (j=2, 4, 6...) 的最后一个有效值
    % 偶数列对应于 [L/M] 型 Padé 逼近
    even_cols = 2:2:N+1;
    last_valid_col = even_cols(end);
    
    % 在表格对角线上寻找最高阶的结果
    val = ep(1, last_valid_col);
end