function res = HEM_DCOPF(mpc, max_order)
%% HEM_DCOPF  ―  HEM Interior-Point Method for DC Optimal Power Flow
%
%  Implements the HEIPM-DCOPF algorithm of Songyan Li (ASU Dissertation,
%  2021), Chapters 5-6.
%
%  KEY STEPS
%  1. Holomorphic embedding of the DCOPF problem with parameter alpha
%  2. Symmetric embedding of inequality constraints
%  3. Recursive computation of Maclaurin series coefficients (eq. 5.126-5.141)
%  4. Padé approximants [M-1/M] evaluated at alpha = 1
%  5. Iterative correction when slack variables turn negative (Ch. 6)
%
%  INPUT
%    mpc        : MATPOWER case struct
%    max_order  : (optional) maximum series order, default = 40
%
%  OUTPUT
%    res.Pg       : generator real power (MW)
%    res.Va_deg   : bus voltage angle (degrees)
%    res.Pf       : branch power flows (MW)
%    res.lmp      : nodal LMPs ($/MWh)
%    res.cost     : total generation cost ($/h)
%    res.all_c    : raw Maclaurin coefficient matrix (diagnostic)
%
%  REFERENCE
%    Li, S. (2021). HEIPM applied to DCOPF. §5.7-5.8, ASU Dissertation.

if nargin < 2, max_order = 40; end

define_constants;

%% ─── 1. Load & convert case ─────────────────────────────────────────────
[baseMVA, bus, gen, branch, gencost] = loadcase(mpc);
[~, bus, gen, branch] = ext2int(bus, gen, branch);
[ref, ~, ~]           = bustypes(bus, gen);

nb = size(bus,    1);
nl = size(branch, 1);
ng = size(gen,    1);

%% ─── 2. Extract system parameters ───────────────────────────────────────
Pd    = bus(:, PD)   / baseMVA;           % nodal load (pu)
PGmin = gen(:, PMIN) / baseMVA;           % gen lower limits (pu)
PGmax = gen(:, PMAX) / baseMVA;           % gen upper limits (pu)

% Line flow limits (pu); 0 in MATPOWER means unlimited → replace with big M
Fmax        = branch(:, RATE_A) / baseMVA;
Fmax(Fmax == 0) = 1e4;

% Quadratic & linear cost coefficients (convert to pu-consistent units)
%   Original gencost: a*(MW)^2 + b*(MW) + c0  [$/h]
%   In pu:            a*baseMVA^2*(pu)^2 + b*baseMVA*(pu) + c0
a_cost = gencost(:, COST)     * baseMVA^2;  % $/h/pu^2
b_cost = gencost(:, COST + 1) * baseMVA;    % $/h/pu
c_cost = gencost(:, COST + 2);              % $/h

% Add small regularisation to prevent near-zero quadratic coefficients
% (linear-only costs make aD singular; see dissertation §5.9)
a_cost = max(a_cost, 1e-4 * baseMVA^2);

%% ─── 3. Build system matrices ────────────────────────────────────────────
gbus = gen(:, GEN_BUS);

% IgG  (NB×NG): 1 if generator j is at bus i
Cg  = full(sparse(gbus, 1:ng, 1, nb, ng));

% IbusG (NG×NB): 1 if gen i at bus j, excluding the slack bus (eq. 5.113-5.114)
IbusG = zeros(ng, nb);
for i = 1:ng
    if gbus(i) ~= ref
        IbusG(i, gbus(i)) = 1;
    end
end

% Branch susceptances and line matrix L (NL×NB)
f_bus    = branch(:, F_BUS);
t_bus    = branch(:, T_BUS);
b_branch = 1 ./ branch(:, BR_X);
b_branch(~isfinite(b_branch)) = 0;

Cft = sparse([(1:nl)';(1:nl)'], [f_bus; t_bus], ...
             [ones(nl,1); -ones(nl,1)], nl, nb);
L   = full(sparse([(1:nl)';(1:nl)'], [f_bus; t_bus], ...
           [b_branch; -b_branch], nl, nb));

% DC susceptance matrix Bbus (NB×NB)
Bbus = full(Cft' * sparse(1:nl,1:nl,b_branch,nl,nl) * Cft);

% Modified B for slack-bus handling (eq. 5.7.1):
%   Row ref → [0…1…0] (enforces θ_ref = 0 as an equality)
B_mod        = Bbus;
B_mod(ref,:) = 0;
B_mod(ref, ref) = 1;
B_mod_T      = B_mod';          % NOTE: B_mod is NOT symmetric

% Zero out slack-bus row of Cg in power-balance equation
%   (power at slack is captured by total-balance constraint, eq. 5.111)
Cg_eq1          = Cg;
Cg_eq1(ref, :)  = 0;
Pd_eq1          = Pd;
Pd_eq1(ref)     = 0;           % no load contribution at slack row

PL_total = sum(Pd);            % total load for eq. 5.141

%% ─── 4. Reference state at alpha = 0 (eq. 5.161) ────────────────────────
%  All voltage angles, power generation, and Lagrange multipliers start at 0.
%  Slack variables are initialised symmetrically about the constraint midpoint.

s1_0     = (PGmax - PGmin) / 2;
s2_0     = s1_0;
s3_0     = Fmax;
s4_0     = Fmax;

sigma1_0 = 1 ./ s1_0;
sigma2_0 = 1 ./ s2_0;
sigma3_0 = 1 ./ s3_0;
sigma4_0 = 1 ./ s4_0;

gamma2_0 = sigma2_0;           % from KKT at alpha=0 (eq. 5.153)
gamma4_0 = sigma4_0;

%% ─── 5. Build constant coefficient matrix A (Appendix A / eq. 5.142) ────
%
%  Variable ordering (total dimension = 2*NB + 7*NG + 6*NL + 1):
%    theta(NB), PG(NG), lambda(NB), gamma1(NG), gamma2(NG),
%    gamma3(NL), gamma4(NL), s1(NG), s2(NG), s3(NL), s4(NL),
%    sigma1(NG), sigma2(NG), sigma3(NL), sigma4(NL), lambda_slack(1)

dim = 2*nb + 7*ng + 6*nl + 1;

% Index ranges (1-based MATLAB columns in A)
idx.theta   = 1:nb;
idx.PG      = nb                + (1:ng);
idx.lambda  = nb + ng           + (1:nb);
idx.gamma1  = 2*nb + ng         + (1:ng);
idx.gamma2  = 2*nb + 2*ng       + (1:ng);
idx.gamma3  = 2*nb + 3*ng       + (1:nl);
idx.gamma4  = 2*nb + 3*ng + nl  + (1:nl);
idx.s1      = 2*nb + 3*ng + 2*nl + (1:ng);
idx.s2      = 2*nb + 4*ng + 2*nl + (1:ng);
idx.s3      = 2*nb + 5*ng + 2*nl + (1:nl);
idx.s4      = 2*nb + 5*ng + 3*nl + (1:nl);
idx.sigma1  = 2*nb + 5*ng + 4*nl + (1:ng);
idx.sigma2  = 2*nb + 6*ng + 4*nl + (1:ng);
idx.sigma3  = 2*nb + 7*ng + 4*nl + (1:nl);
idx.sigma4  = 2*nb + 7*ng + 5*nl + (1:nl);
idx.lslack  = 2*nb + 7*ng + 6*nl + 1;

% Row ranges mirror variable indices exactly (square system)
row.eq1  = idx.theta;          % Bθ − Cg·PG = rhs         (5.126)
row.eq2  = idx.PG;             % PG + s1 = rhs             (5.127)
row.eq3  = idx.gamma1;         % s1 + s2 = 0               (5.128)
row.eq4  = idx.gamma3;         % Lθ + s3 = 0               (5.129)
row.eq5  = idx.gamma4;         % s3 + s4 = 0               (5.130)
row.eq6  = idx.gamma2;         % ∇PG KKT                   (5.131)
row.eq7  = idx.lambda;         % B^T λ + L^T γ3 = 0        (5.132)
row.eq8  = idx.s1;             % γ1+γ2−σ1 = rhs            (5.133)
row.eq9  = idx.s2;             % γ2−σ2 = rhs               (5.134)
row.eq10 = idx.s3;             % γ3+γ4−σ3 = rhs            (5.135)
row.eq11 = idx.s4;             % γ4−σ4 = rhs               (5.136)
row.eq12 = idx.sigma1;         % Dσ1·s1 + Ds1·σ1 = rhs    (5.137)
row.eq13 = idx.sigma2;         % Dσ2·s2 + Ds2·σ2 = rhs    (5.138)
row.eq14 = idx.sigma3;         % Dσ3·s3 + Ds3·σ3 = rhs    (5.139)
row.eq15 = idx.sigma4;         % Dσ4·s4 + Ds4·σ4 = rhs    (5.140)
row.eq16 = idx.lslack;         % e^T PG = rhs              (5.141)

A = zeros(dim, dim);

% Eq1: B_mod·θ[n] − Cg_eq1·PG[n] = rhs
A(row.eq1,  idx.theta)  = B_mod;
A(row.eq1,  idx.PG)     = -Cg_eq1;

% Eq2: PG[n] + s1[n] = rhs
A(row.eq2,  idx.PG)     = eye(ng);
A(row.eq2,  idx.s1)     = eye(ng);

% Eq3: s1[n] + s2[n] = 0
A(row.eq3,  idx.s1)     = eye(ng);
A(row.eq3,  idx.s2)     = eye(ng);

% Eq4: L·θ[n] + s3[n] = 0
A(row.eq4,  idx.theta)  = L;
A(row.eq4,  idx.s3)     = eye(nl);

% Eq5: s3[n] + s4[n] = 0
A(row.eq5,  idx.s3)     = eye(nl);
A(row.eq5,  idx.s4)     = eye(nl);

% Eq6: aD·PG[n] − IbusG·λ[n] + γ1[n] + e·λslack[n] = rhs  (5.131)
A(row.eq6,  idx.PG)     = diag(a_cost);
A(row.eq6,  idx.lambda) = -IbusG;
A(row.eq6,  idx.gamma1) = eye(ng);
A(row.eq6,  idx.lslack) = ones(ng, 1);

% Eq7: B_mod^T·λ[n] + L^T·γ3[n] = 0  (5.132)
A(row.eq7,  idx.lambda) = B_mod_T;
A(row.eq7,  idx.gamma3) = L';

% Eq8: γ1[n] + γ2[n] − σ1[n] = −σ1[n-1]
A(row.eq8,  idx.gamma1) = eye(ng);
A(row.eq8,  idx.gamma2) = eye(ng);
A(row.eq8,  idx.sigma1) = -eye(ng);

% Eq9: γ2[n] − σ2[n] = −σ2[n-1]
A(row.eq9,  idx.gamma2) = eye(ng);
A(row.eq9,  idx.sigma2) = -eye(ng);

% Eq10: γ3[n] + γ4[n] − σ3[n] = −σ3[n-1]
A(row.eq10, idx.gamma3) = eye(nl);
A(row.eq10, idx.gamma4) = eye(nl);
A(row.eq10, idx.sigma3) = -eye(nl);

% Eq11: γ4[n] − σ4[n] = −σ4[n-1]
A(row.eq11, idx.gamma4) = eye(nl);
A(row.eq11, idx.sigma4) = -eye(nl);

% Eq12-15: Dσi[0]·si[n] + Dsi[0]·σi[n] = conv_rhs   (5.137-5.140)
A(row.eq12, idx.s1)     = diag(sigma1_0);
A(row.eq12, idx.sigma1) = diag(s1_0);

A(row.eq13, idx.s2)     = diag(sigma2_0);
A(row.eq13, idx.sigma2) = diag(s2_0);

A(row.eq14, idx.s3)     = diag(sigma3_0);
A(row.eq14, idx.sigma3) = diag(s3_0);

A(row.eq15, idx.s4)     = diag(sigma4_0);
A(row.eq15, idx.sigma4) = diag(s4_0);

% Eq16: e^T·PG[n] = δn1·PL_total
A(row.eq16, idx.PG)     = ones(1, ng);

fprintf('  [HEM-DCOPF] System size: %d vars. Condition of A: %.2e\n', ...
        dim, condest(A));

%% ─── 6. Initialise coefficient storage ──────────────────────────────────
%  all_c(:, n+1) stores the n-th order Maclaurin coefficient (n = 0,1,…)
all_c = zeros(dim, max_order + 1);

% n = 0  (reference state)
all_c(idx.s1,     1) = s1_0;
all_c(idx.s2,     1) = s2_0;
all_c(idx.s3,     1) = s3_0;
all_c(idx.s4,     1) = s4_0;
all_c(idx.sigma1, 1) = sigma1_0;
all_c(idx.sigma2, 1) = sigma2_0;
all_c(idx.sigma3, 1) = sigma3_0;
all_c(idx.sigma4, 1) = sigma4_0;
all_c(idx.gamma2, 1) = gamma2_0;
all_c(idx.gamma4, 1) = gamma4_0;
% All other variables are 0 at n=0

%% ─── 7. Maclaurin recursion (eq. 5.126-5.141) ───────────────────────────
for n = 1:max_order
    rhs = zeros(dim, 1);

    % Eq1 RHS: −δn1·Pd_eq1
    if n == 1, rhs(row.eq1) = -Pd_eq1; end

    % Eq2 RHS: (PGmax+PGmin)/2 · δn1
    if n == 1, rhs(row.eq2) = (PGmax + PGmin) / 2; end

    % Eq3-5 RHS: zero

    % Eq6 RHS: −δn1·b_cost
    if n == 1, rhs(row.eq6) = -b_cost; end

    % Eq7 RHS: zero

    % Eq8-11 RHS: −σi[n-1]   (stored at column n in 1-based indexing)
    rhs(row.eq8)  = -all_c(idx.sigma1, n);
    rhs(row.eq9)  = -all_c(idx.sigma2, n);
    rhs(row.eq10) = -all_c(idx.sigma3, n);
    rhs(row.eq11) = -all_c(idx.sigma4, n);

    % Eq12-15 RHS: −Σ_{k=1}^{n-1} si[k] .* σi[n-k]
    %   Element-wise convolution of product si(α)·σi(α) = 1
    conv12 = zeros(ng, 1);   conv13 = zeros(ng, 1);
    conv14 = zeros(nl, 1);   conv15 = zeros(nl, 1);
    for k = 1:n-1
        kp1 = k + 1;   nk1 = n - k + 1;   % 1-based coefficient indices
        conv12 = conv12 + all_c(idx.s1, kp1) .* all_c(idx.sigma1, nk1);
        conv13 = conv13 + all_c(idx.s2, kp1) .* all_c(idx.sigma2, nk1);
        conv14 = conv14 + all_c(idx.s3, kp1) .* all_c(idx.sigma3, nk1);
        conv15 = conv15 + all_c(idx.s4, kp1) .* all_c(idx.sigma4, nk1);
    end
    rhs(row.eq12) = -conv12;
    rhs(row.eq13) = -conv13;
    rhs(row.eq14) = -conv14;
    rhs(row.eq15) = -conv15;

    % Eq16 RHS: δn1·PL_total
    if n == 1, rhs(row.eq16) = PL_total; end

    % Solve for n-th coefficients
    all_c(:, n + 1) = A \ rhs;
end

%% ─── 8. Padé approximants → evaluate at alpha = 1 ───────────────────────
M = floor(max_order / 2);      % [M-1/M] near-diagonal PA

PG_pade    = zeros(ng, 1);
theta_pade = zeros(nb, 1);
lam_pade   = zeros(nb, 1);

for i = 1:ng
    PG_pade(i)    = pade_eval_hem(all_c(idx.PG(i),     :), M);
end
for i = 1:nb
    theta_pade(i) = pade_eval_hem(all_c(idx.theta(i),  :), M);
    lam_pade(i)   = pade_eval_hem(all_c(idx.lambda(i), :), M);
end

%% ─── 9. Iterative correction (Chapter 6) ─────────────────────────────────
%  If any PG falls outside [PGmin, PGmax], fix it at the violated bound,
%  remove the corresponding generator from the free set, and re-solve.
%  This corrects the branch-cut crossing described in §6.2.

max_iter = 5;
fixed_g  = false(ng, 1);
PG_fixed = zeros(ng, 1);

for iter = 1:max_iter
    violated = false;
    for i = 1:ng
        if ~fixed_g(i)
            if PG_pade(i) < PGmin(i) - 1e-6
                PG_fixed(i) = PGmin(i);
                fixed_g(i)  = true;
                violated     = true;
            elseif PG_pade(i) > PGmax(i) + 1e-6
                PG_fixed(i) = PGmax(i);
                fixed_g(i)  = true;
                violated     = true;
            end
        end
    end

    if ~violated, break; end

    % Re-solve with fixed generators: reduce effective load for free gens
    if any(~fixed_g)
        free_idx  = find(~fixed_g);
        fixed_idx = find(fixed_g);

        % Adjust total load for free generators
        fixed_inj = sum(PG_fixed(fixed_idx));
        PL_adj    = Pd;
        for fi = fixed_idx(:)'
            bus_fi = gbus(fi);
            PL_adj(bus_fi) = PL_adj(bus_fi) - PG_fixed(fi);
        end
        PL_adj(ref) = 0;

        % Build reduced system for free generators only
        ng_f = length(free_idx);
        Cg_f = Cg_eq1(:, free_idx);
        IbusG_f = IbusG(free_idx, :);
        a_f = a_cost(free_idx);
        b_f = b_cost(free_idx);
        PGmin_f = PGmin(free_idx);
        PGmax_f = PGmax(free_idx);
        PL_total_f = sum(PL_adj);

        res_f = hem_dcopf_core(B_mod, B_mod_T, L, Cg_f, IbusG_f, ...
                               a_f, b_f, PL_adj, PL_total_f, ...
                               PGmin_f, PGmax_f, Fmax, ...
                               nb, ng_f, nl, ref, max_order);

        PG_pade(free_idx) = res_f.PG_pu;
        theta_pade        = res_f.theta_pu;
    end
end

% Final clip to limits (safety)
PG_pade = max(PGmin, min(PGmax, PG_pade));

%% ─── 10. Assemble output ─────────────────────────────────────────────────
res.Pg      = PG_pade  * baseMVA;            % MW
res.PG_pu   = PG_pade;
res.Va      = theta_pade;                     % radians
res.Va_deg  = theta_pade * 180 / pi;         % degrees
res.Pf      = L * theta_pade * baseMVA;      % MW (from-bus flows)
res.lmp     = lam_pade * baseMVA;            % $/MWh
res.cost    = sum(gencost(:,COST)   .* res.Pg.^2 + ...
                  gencost(:,COST+1) .* res.Pg    + ...
                  gencost(:,COST+2));              % $/h
res.all_c   = all_c;
res.M_pade  = M;

fprintf('  [HEM-DCOPF] Cost = %.4f $/h  |  PA order [%d/%d]\n', ...
        res.cost, M-1, M);
end


%% ═══════════════════════════════════════════════════════════════════════════
function res = hem_dcopf_core(B_mod, B_mod_T, L, Cg, IbusG, ...
                              a_cost, b_cost, Pd_eq1, PL_total, ...
                              PGmin, PGmax, Fmax, nb, ng, nl, ref, max_order)
%% hem_dcopf_core  ─  inner HEIPM-DCOPF solver for reduced generator set
%  Used by the iterative correction loop in HEM_DCOPF.

s1_0 = (PGmax - PGmin)/2;  s2_0 = s1_0;
s3_0 = Fmax;               s4_0 = Fmax;
sigma1_0 = 1./s1_0;  sigma2_0 = 1./s2_0;
sigma3_0 = 1./s3_0;  sigma4_0 = 1./s4_0;

dim = 2*nb + 7*ng + 6*nl + 1;

ix.theta  = 1:nb;
ix.PG     = nb               + (1:ng);
ix.lambda = nb+ng            + (1:nb);
ix.g1     = 2*nb+ng          + (1:ng);
ix.g2     = 2*nb+2*ng        + (1:ng);
ix.g3     = 2*nb+3*ng        + (1:nl);
ix.g4     = 2*nb+3*ng+nl     + (1:nl);
ix.s1     = 2*nb+3*ng+2*nl   + (1:ng);
ix.s2     = 2*nb+4*ng+2*nl   + (1:ng);
ix.s3     = 2*nb+5*ng+2*nl   + (1:nl);
ix.s4     = 2*nb+5*ng+3*nl   + (1:nl);
ix.si1    = 2*nb+5*ng+4*nl   + (1:ng);
ix.si2    = 2*nb+6*ng+4*nl   + (1:ng);
ix.si3    = 2*nb+7*ng+4*nl   + (1:nl);
ix.si4    = 2*nb+7*ng+5*nl   + (1:nl);
ix.ls     = 2*nb+7*ng+6*nl+1;

rw = ix;  % rows same as columns for square system

A = zeros(dim, dim);
A(rw.theta, ix.theta)  = B_mod;
A(rw.theta, ix.PG)     = -Cg;
A(rw.PG,    ix.PG)     = eye(ng);  A(rw.PG,  ix.s1)  = eye(ng);
A(rw.g1,    ix.s1)     = eye(ng);  A(rw.g1,  ix.s2)  = eye(ng);
A(rw.g3,    ix.theta)  = L;        A(rw.g3,  ix.s3)  = eye(nl);
A(rw.g4,    ix.s3)     = eye(nl);  A(rw.g4,  ix.s4)  = eye(nl);
A(rw.g2,    ix.PG)     = diag(a_cost);
A(rw.g2,    ix.lambda) = -IbusG;
A(rw.g2,    ix.g1)     = eye(ng);
A(rw.g2,    ix.ls)     = ones(ng,1);
A(rw.lambda,ix.lambda) = B_mod_T;
A(rw.lambda,ix.g3)     = L';
A(rw.s1, ix.g1) = eye(ng); A(rw.s1, ix.g2) = eye(ng); A(rw.s1, ix.si1) = -eye(ng);
A(rw.s2, ix.g2) = eye(ng); A(rw.s2, ix.si2) = -eye(ng);
A(rw.s3, ix.g3) = eye(nl); A(rw.s3, ix.g4) = eye(nl); A(rw.s3, ix.si3) = -eye(nl);
A(rw.s4, ix.g4) = eye(nl); A(rw.s4, ix.si4) = -eye(nl);
A(rw.si1, ix.s1) = diag(sigma1_0); A(rw.si1, ix.si1) = diag(s1_0);
A(rw.si2, ix.s2) = diag(sigma2_0); A(rw.si2, ix.si2) = diag(s2_0);
A(rw.si3, ix.s3) = diag(sigma3_0); A(rw.si3, ix.si3) = diag(s3_0);
A(rw.si4, ix.s4) = diag(sigma4_0); A(rw.si4, ix.si4) = diag(s4_0);
A(rw.ls,  ix.PG) = ones(1, ng);

all_c = zeros(dim, max_order+1);
all_c(ix.s1,1)=s1_0; all_c(ix.s2,1)=s2_0; all_c(ix.s3,1)=s3_0; all_c(ix.s4,1)=s4_0;
all_c(ix.si1,1)=sigma1_0; all_c(ix.si2,1)=sigma2_0;
all_c(ix.si3,1)=sigma3_0; all_c(ix.si4,1)=sigma4_0;
all_c(ix.g2,1)=sigma2_0;  all_c(ix.g4,1)=sigma4_0;

for n = 1:max_order
    rhs = zeros(dim,1);
    if n==1
        rhs(rw.theta) = -Pd_eq1;
        rhs(rw.PG)    = (PGmax+PGmin)/2;
        rhs(rw.g2)    = -b_cost;
        rhs(rw.ls)    = PL_total;
    end
    rhs(rw.s1)  = -all_c(ix.si1, n);
    rhs(rw.s2)  = -all_c(ix.si2, n);
    rhs(rw.s3)  = -all_c(ix.si3, n);
    rhs(rw.s4)  = -all_c(ix.si4, n);

    c12=zeros(ng,1); c13=zeros(ng,1); c14=zeros(nl,1); c15=zeros(nl,1);
    for k = 1:n-1
        kp1=k+1; nk1=n-k+1;
        c12=c12+all_c(ix.s1,kp1).*all_c(ix.si1,nk1);
        c13=c13+all_c(ix.s2,kp1).*all_c(ix.si2,nk1);
        c14=c14+all_c(ix.s3,kp1).*all_c(ix.si3,nk1);
        c15=c15+all_c(ix.s4,kp1).*all_c(ix.si4,nk1);
    end
    rhs(rw.si1)=-c12; rhs(rw.si2)=-c13; rhs(rw.si3)=-c14; rhs(rw.si4)=-c15;

    all_c(:, n+1) = A \ rhs;
end

M = floor(max_order/2);
res.PG_pu    = arrayfun(@(i) pade_eval_hem(all_c(ix.PG(i),:),    M), 1:ng)';
res.theta_pu = arrayfun(@(i) pade_eval_hem(all_c(ix.theta(i),:), M), 1:nb)';

end