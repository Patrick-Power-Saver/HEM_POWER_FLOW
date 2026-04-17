function val = pade_eval_hem(c_arr, M)
%% pade_eval_hem  ―  [M-1/M] Padé approximant at alpha = 1
%
%  Implements the Direct Method (DM) described in Li (2021) §3.1.3-3.1.4.
%
%  INPUT
%    c_arr  : row/column vector of Maclaurin coefficients in 1-based MATLAB
%             ordering, i.e.  c_arr(k) = c[k-1]  (the (k-1)-th order coeff)
%    M      : denominator polynomial degree; numerator degree = M-1
%
%  OUTPUT
%    val    : value of the [M-1/M] PA evaluated at alpha = 1
%
%  THEORY (from eq. 3.28-3.29 in dissertation, with L = M-1)
%
%  Denominator coefficients b[1],...,b[M]  (b[0] = 1) are solved from
%
%      C * [b1;...;bM] = d ,   where
%
%      C(i,j)  = c[ M-1-(i-1)+(j-1) ] = c[ M-i+j-1 ]   →  c_arr(M-i+j)
%      d(i)    = -c[ (M-1)+M+1-i ]     = -c[ 2M-i   ]   → -c_arr(2M-i+1)
%
%  Numerator coefficients a[0],...,a[M-1]:
%
%      a[k] = Σ_{j=0}^{min(k,M)}  c[k-j] * b[j]   (b[0]=1)
%
%  PA at alpha=1:   val = Σ a[k] / Σ b[j]

n_c = length(c_arr);        % number of available coefficients

% ── Adjust M if not enough coefficients ───────────────────────────────
while M > 1 && 2*M > n_c - 1
    M = M - 1;
end
if M <= 0
    val = sum(c_arr);
    return
end

L = M - 1;           % numerator polynomial degree
c = c_arr(:);        % ensure column vector, 1-indexed (c(k) = c[k-1])

% ── Build Hankel coefficient matrix (M x M) ───────────────────────────
Ch = zeros(M, M);
for i = 1:M
    for j = 1:M
        idx = M - i + j;          % 1-based index into c
        if idx >= 1 && idx <= n_c
            Ch(i, j) = c(idx);
        end
    end
end

% ── Build RHS vector d (M x 1) ────────────────────────────────────────
d = zeros(M, 1);
for i = 1:M
    idx = 2*M - i + 1;           % 1-based index: c[2M-i]
    if idx >= 1 && idx <= n_c
        d(i) = -c(idx);
    end
end

% ── Solve for denominator coefficients ────────────────────────────────
if rcond(Ch) < 1e-14
    % Near-singular: fall back to partial-sum evaluation
    val = sum(c_arr(1:min(n_c, 2*M + 1)));
    return
end

b_coeff = Ch \ d;               % b[1],...,b[M]
b_full  = [1; b_coeff];         % b[0]=1, b[1],...,b[M]

% ── Compute numerator coefficients a[0],...,a[M-1] (eq. 3.29) ─────────
a_coeff = zeros(L + 1, 1);
for k = 0:L
    for j = 0:min(k, M)
        ic = k - j + 1;          % c_arr index for c[k-j]
        ib = j + 1;              % b_full index for b[j]
        if ic >= 1 && ic <= n_c
            a_coeff(k + 1) = a_coeff(k + 1) + c(ic) * b_full(ib);
        end
    end
end

% ── Evaluate PA at alpha = 1 ──────────────────────────────────────────
num_val = sum(a_coeff);
den_val = sum(b_full);

if abs(den_val) < 1e-12
    val = sum(c_arr(1:min(n_c, 2*M + 1)));   % fallback
else
    val = num_val / den_val;
end

end