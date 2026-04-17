function res = Lossless_dcopf(mpc, varargin)
% DC OPF

%% 导入系统数据
define_constants;
[baseMVA, bus, gen, branch, gencost] = loadcase(mpc);
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
[ref, pv, pq] = bustypes(bus, gen); % reference bus index
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);

%% 定义常量
nb = size(bus, 1);          % number of buses
nl = size(branch, 1);       % number of lines
ng = size(gen, 1);          % number of gens

Pd = bus(:, PD);            % bus load (MW)
Fmax = branch(:, RATE_A);   % branch flow limit (MW)
Pgmin = gen(:, PMIN);       % gen max. power output (MW)
Pgmax = gen(:, PMAX);       % gen max. power output (MW)
C2 = gencost(:, COST);                          % gen injection cost ($/MWh^2)
C1 = gencost(:, COST + 1);                      % gen injection cost ($/MWh)
C0 = gencost(:, COST + 2);                      % gen injection cost ($/h)

% 单位换算
% -------- unit conversion --------- %
Pd = Pd / baseMVA;
Fmax = Fmax / baseMVA;          % p.u.
Pgmin = Pgmin / baseMVA;        % p.u.
Pgmax = Pgmax / baseMVA;        % p.u.
C2 = C2 * (baseMVA^2);          % $/h
C1 = C1 * baseMVA;              % $/h

%% 建立矩阵
gbus = gen(:, GEN_BUS);                 % connection matrix for gen
Cg = sparse(gbus, (1:ng)', 1, nb, ng);  % element i, j is 1 if gen(j) at bus i

f = branch(:, F_BUS);               % connection matrix for line and from - to buses
t = branch(:, T_BUS);               % element k, i is 1 if branch k connects "from" bus i
i = [(1:nl)'; (1:nl)'];             % element k, j is -1 if branch k connects "to" bus j
Cf = sparse((1:nl)', f, 1, nl, nb);
Ct = sparse((1:nl)', t, 1, nl, nb);
Cft = full(sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb));
Ctf = sparse(i, [t; f], [ones(nl, 1); -ones(nl, 1)], nl, nb);
b = 1 ./ branch(:, BR_X);%branch(:, BR_X); imag(Ybus)
Bf = spdiags(b, 0, nl, nl) * Cft;   % Bf * Va is the vector of real branch power (p.u.)
Bt = spdiags(b, 0, nl, nl) * Ctf;
Bbus = Cft' * Bf;                   % Bbus * Va is the vector of nodal real power injection (p.u.)
%Bbus = Ctf' * Bt;
% Gbus = real(Ybus);
%PTDF = makePTDF(mpc);

%% 定义决策变量
Pg = sdpvar(ng, 1, 'full');     % generator real power injection (MW)
Va = sdpvar(nb, 1, 'full');     % bus voltage angle (rad)
Pbus = sdpvar(nb, 1, 'full');
Pf =  sdpvar(nl, 1, 'full');
Pt =  sdpvar(nl, 1, 'full');

% -------- assign initial value -------- %
% mpopt = mpoption( 'out.all', 0);
% results = rundcopf(mpc, mpopt);
% Vm_sol = results.bus(:, VM);                    % voltage magnitude (p.u.)
% Va_sol = results.bus(:, VA) * (pi/180);         % voltage angle (rad)
% [Vd_sol, Vq_sol] = pol2cart(Va_sol, Vm_sol);    % complex bus voltage in cartesian coord. (p.u.)
% Pg_sol = results.gen(:, PG) / baseMVA;          % generator active power injection (p.u.)
% Qg_sol = results.gen(:, QG) / baseMVA;          % generator reactive power injection (p.u.)
% assign(Va, Vd_sol + 1i * Vq_sol);
% assign(Pg, Pg_sol + 1i * Qg_sol);

%% 目标函数
%Obj = C1' * Pg;            % linear cost function ($/h)
Obj = C2' * (Pg.^2) + C1' * Pg + ones(ng, 1)' * C0;       
Const = [];

%% 节点功率平衡等式约束
Const = [Const, Pbus + Pd - Cg * Pg == 0];   % nodal real power balance (MW)
Const = [Const, Pbus - Bbus * Va == 0];
%Const = [Const, Pf == PTDF * Pbus];

%% 支路功率与节点电压的等式约束
Const = [Const, Bf * Va == Pf];
Const = [Const, Bt * Va == Pt];
Const = [Const,  Va(ref) == 0];                         % reference bus (rad)
Const = [Const, Pgmin <= Pg];                      % gen real power output limit (MW)
Const = [Const, Pg <= Pgmax]; 

%% 线路潮流上下限不等式约束 （配电网不需要此约束）
if (nb < nl)
    Const = [Const, Fmax >= Pf];        % real power branch flow limit (MW)
    Const = [Const, Pf >= -Fmax];
    Const = [Const, Fmax >= Pt ];        % real power branch flow limit (MW)
    Const = [Const, Pt >= -Fmax];
else
end

Opts = sdpsettings('solver', 'gurobi', 'verbose', 0);
optimize(Const, Obj, Opts);
check(Const);
   
res.Pg = value(Pg)* baseMVA;         % gen real power output (MW)
res.Va = value(Va)/pi*180;         % bus voltage angle (rad)
res.Pbus = value(Pbus);
res.Pf = value(Pf) * baseMVA;
res.Pt = value(Pt) * baseMVA;
res.lmp = dual(Const(1));   % lmp($/MWh)
res.cost = value(Obj);      % gen real power output cost ($/h)
disp('done');
end