function [A, B_vec] = construct_aug_system(G, Btr, rhs, pq, pv, ref, W0)
    % 强制转换为列向量
    pq = pq(:); pv = pv(:);
    n_pq = length(pq); n_pv = length(pv);
    dim = 2*n_pq + 2*n_pv;
    
    A = zeros(dim, dim);
    B_vec = zeros(dim, 1);
    
    % 索引映射: x = [Vre_pq, Vim_pq, Vim_pv, Q_pv]
    re_pq_cols = 1:n_pq;
    im_pq_cols = n_pq+1 : 2*n_pq;
    im_pv_cols = 2*n_pq+1 : 2*n_pq+n_pv;
    q_pv_cols  = 2*n_pq+n_pv+1 : dim;

    % --- 1. PQ 节点方程行 ---
    % 实部方程行 (1:n_pq): G*Vre - B*Vim = Re(rhs)
    A(1:n_pq, re_pq_cols) = G(pq, pq);
    A(1:n_pq, im_pq_cols) = -Btr(pq, pq);
    A(1:n_pq, im_pv_cols) = -Btr(pq, pv);
    
    % 虚部方程行 (n_pq+1:2*n_pq): B*Vre + G*Vim = Im(rhs)
    A(n_pq+1:2*n_pq, re_pq_cols) = Btr(pq, pq);
    A(n_pq+1:2*n_pq, im_pq_cols) = G(pq, pq);
    A(n_pq+1:2*n_pq, im_pv_cols) = G(pq, pv);

    % --- 2. PV 节点方程行 ---
    % 实部方程行 (2*n_pq+1 : 2*n_pq+n_pv)
    row_pv_re = 2*n_pq+1 : 2*n_pq+n_pv;
    A(row_pv_re, re_pq_cols) = G(pv, pq);
    A(row_pv_re, im_pq_cols) = -Btr(pv, pq);
    A(row_pv_re, im_pv_cols) = -Btr(pv, pv);
    
    % 虚部方程行 (2*n_pq+n_pv+1 : dim)
    row_pv_im = 2*n_pq+n_pv+1 : dim;
    A(row_pv_im, re_pq_cols) = Btr(pv, pq);
    A(row_pv_im, im_pq_cols) = G(pv, pq);
    A(row_pv_im, im_pv_cols) = G(pv, pv);
    % 核心修正：Q 的系数是实部 W0 的对角阵
    A(row_pv_im, q_pv_cols)  = diag(real(W0(pv)));

    % 填充 RHS
    B_vec(1:n_pq) = real(rhs(pq));
    B_vec(n_pq+1:2*n_pq) = imag(rhs(pq));
    B_vec(row_pv_re) = real(rhs(pv));
    B_vec(row_pv_im) = imag(rhs(pv));
end