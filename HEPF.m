function [VV1, Vh, Qh, k, delS] = HEPF(Ybus, Sbus, V0, ref, pv, pq, tol)
    % 输入说明:
    % Vsp: PV节点的电压目标幅值向量 (长度应与pv一致)
    Vsp = V0(pv);
    n_bus = length(Sbus);
    max_k = 20; % 增加阶数以提高精度
    Vh = zeros(n_bus, max_k);
    W  = zeros(n_bus, max_k);
    Qh = zeros(n_bus, max_k);
    
    % 导纳矩阵拆分
    Ysh = sum(Ybus, 2);
    Ytr = Ybus - diag(Ysh);
    Gtr = real(Ytr); Btr = imag(Ytr);
    
    % Step 1: 0阶初始化 (Germ Solution)
    Vh(:, 1) = 1.0; 
    W(:, 1) = 1.0;
    k = 2;%幂级数的阶数
    delS = 1;

    while (delS > tol) && (k <= max_k)
        n = k - 1; % 当前阶数
        
        % --- 计算 PV 节点已知的实部 V_re[n] (公式 3) ---
        V_re_known = zeros(n_bus, 1);
        if n == 1
            V_re_known(pv) = (Vsp.^2 - 1) / 2;
            V_re_known(ref) = real(V0(ref)) - 1; % 公式 4
        else
            % 卷积项: -0.5 * sum_{m=1}^{n-1} V[m] * conj(V[n-m])
            for i = 1:length(pv)
                idx = pv(i);
                conv_val = 0;
                for m = 1:(n-1)
                    conv_val = conv_val + Vh(idx, m+1) * conj(Vh(idx, k-m));
                end
                V_re_known(idx) = -0.5 * real(conv_val);
            end
            V_re_known(ref) = 0;
        end
        
        % 预填实部到 Vh，方便后续移项计算
        Vh(pv, k) = V_re_known(pv);
        Vh(ref, k) = V_re_known(ref);

        % --- 计算 RHS (公式 1 & 2) ---
        rhs_complex = zeros(n_bus, 1);
        % PQ 节点
        rhs_complex(pq) = conj(Sbus(pq)) .* conj(W(pq, k-1)) - Ysh(pq) .* Vh(pq, k-1);
        
        % PV 节点 (处理无功级数卷积项)
        if ~isempty(pv)
            conv_qw = zeros(length(pv), 1);
            for i = 1:length(pv)
                idx = pv(i);
                for m = 1:(n-1)
                    conv_qw(i) = conv_qw(i) + Qh(idx, m+1) * conj(W(idx, k-m));
                end
            end
            rhs_complex(pv) = real(Sbus(pv)) .* conj(W(pv, k-1)) - Ysh(pv) .* Vh(pv, k-1) - 1j * conv_qw;
        end
        
        % --- 关键：移项补偿 ---
        % 将左侧已知的 V_re[n] 贡献移到右侧
        % rhs_final = rhs_complex - Ytr * V_re_known
        rhs_final = rhs_complex - Ytr(:, :) * Vh(:, k);

        % --- Step 2: 求解 2N 阶实数增广矩阵 ---
        [A_aug, B_aug] = construct_aug_system(Gtr, Btr, rhs_final, pq, pv, ref, W(:, 1));
        
        % 检查奇异性
        if cond(A_aug) > 1e15
            error('增广矩阵奇异，请检查网络拓扑或节点分类。');
        end
        
        x = A_aug \ B_aug;

        % --- Step 3: 回填解出的系数 ---
        n_pq = length(pq);
        n_pv = length(pv);
        % 回填 PQ
        Vh(pq, k) = x(1:n_pq) + 1j * x(n_pq+1 : 2*n_pq);
        % 回填 PV (实部已在上面预填，这里填虚部)
        Vh(pv, k) = Vh(pv, k) + 1j * x(2*n_pq+1 : 2*n_pq+n_pv);
        % 回填 Q 系数
        Qh(pv, k) = x(2*n_pq+n_pv+1 : end);

        % --- Step 4: 更新倒数级数 W (公式 V*W = 1) ---
        temp_w_sum = zeros(n_bus, 1);
        for m = 0:(n-1)
            temp_w_sum = temp_w_sum + W(:, m+1) .* Vh(:, k-m);
        end
        W(:, k) = -temp_w_sum ./ Vh(:, 1);

        % --- Step 5: 计算残差 ---
        VV1 = sum(Vh(:, 1:k), 2);
        Scal = VV1 .* conj(Ybus * VV1);
        delS = max(abs(Sbus([pq;pv]) - Scal([pq;pv])));
        
        if delS < tol, break; end
        k = k + 1;
    end
end