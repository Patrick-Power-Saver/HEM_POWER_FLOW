function [MVAbase, bus, gen, gencost, branch, f, success, et] = Lossless_Linear_opf(mpc, varargin)
% Lossless Linear ACOPF
% 11-05-2022

%% 导入系统数据
define_constants;
mpopt = mpoption;
[baseMVA, bus, gen, branch, gencost] = loadcase(mpc);
[i2e, bus, gen, branch] = ext2int(bus, gen, branch);
[ref, pv, pq] = bustypes(bus, gen); % reference bus index
t0 = clock;

%% 参数设置
branch(:, BR_STATUS) = 1;
br = find(branch(:, BR_STATUS));            %% in-service branches
no =bus(:,1);
nb = size(bus, 1);          % number of buses
nl = size(branch, 1);       % number of lines
ng = size(gen, 1);          % number of gens
wbuses = find((bus(:,2) == PV |bus(:,2) == PQ ) & (bus(:,PD) > 0));%bus(:,2) == REF |

Pd = bus(:, PD);              % bus load (MVA)
Qd = bus(:, QD);
Fmax = branch(:, RATE_A);                       % branch flow limit (MVA)

Vmin = bus(:, VMIN);                            % minimum voltage magnitude (p.u.)
Vmax = bus(:, VMAX);                            % maximum voltage magnitude (p.u.)
Sgmin = gen(:, PMIN) + 1i * gen(:, QMIN);       % gen max. power output (MVA)
Sgmax = gen(:, PMAX) + 1i * gen(:, QMAX);       % gen max. power output (MVA)
C2 = gencost(:, COST);                          % gen injection cost ($/MWh^2)
C1 = gencost(:, COST + 1);                      % gen injection cost ($/MWh)
C0 = gencost(:, COST + 2);                      % gen injection cost ($/h)

% 节点-支路关联矩阵和节点导纳矩阵
gbus = gen(:, GEN_BUS);
Cg = full(sparse(gbus, (1:ng)', 1, nb, ng));
Cw = sparse(wbuses, (1:length(wbuses))', 1, nb, length(wbuses));
f = branch(:, F_BUS);
t = branch(:, T_BUS);
i = [(1:nl)'; (1:nl)'];             % element k, j is -1 if branch k connects "to" bus j
Cf = full(sparse((1:nl)', f, 1, nl, nb));
Ct = full(sparse((1:nl)', t, 1, nl, nb));
Cft = full(sparse(i, [f; t], [ones(nl, 1); -ones(nl, 1)], nl, nb));
Ctf = sparse(i, [t; f], [ones(nl, 1); -ones(nl, 1)], nl, nb);

% H = makePTDF(mpc);
% Hg = H(find(branch(:,RATE_A) > 0), gbus);
% Hd = H(find(branch(:,RATE_A) > 0),find(bus(:,PD) > 0));
% Pd_hat = Pd(find(bus(:,PD) > 0));Qd_hat = Qd(find(bus(:,PD) > 0));

[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
Ysh = (bus(:, GS) + 1j * bus(:, BS)) / baseMVA;
Ysc = 1 ./ (branch(:, BR_R) - 1j * branch(:, BR_X));
Gbus = real(Ybus);
Bbus = imag(Ybus);
Gsh = real(Ysh);
Bsh = imag(Ysh);

GP = Gbus; %导纳矩阵的实部
BD = diag(sum(Bbus));
B_mod= Bbus-BD;% B'
BP = -B_mod; %BP=-B'
BQ = -Bbus; %BQ=-B
GQ = -Gbus;                       %GQ approxiately equals -Gbus
Xp_dlpf = [BP GP];
Xq_dlpf = [GQ BQ];
CJ = [Xp_dlpf;Xq_dlpf];

%单位换算
%-------- unit conversion --------- %
Pd = Pd / baseMVA;
Qd = Qd / baseMVA;
Fmax = Fmax / baseMVA;          % p.u.
% FLmax = sparse([f;t],[t,f],[Fmax;-Fmax],nb,nb);
Sgmin = Sgmin / baseMVA;        % p.u.
Sgmax = Sgmax / baseMVA*0.8;        % p.u.
C2 = C2 * (baseMVA^2);          % $/h
C1 = C1 * baseMVA;              % $/h
b = 0;

phi = 0.95;%风电功率因数

%% 定义变量
Vm = sdpvar(nb, 1, 'full');
Vm_2 = sdpvar(nb, 1, 'full');
Va = sdpvar(nb,1, 'full');
Pbus = sdpvar(nb,1, 'full');
Qbus = sdpvar(nb,1, 'full');
Pf = sdpvar(nl,1, 'full');
Qf = sdpvar(nl,1, 'full');
Pt = sdpvar(nl,1, 'full');
Qt = sdpvar(nl,1, 'full');
Pg = sdpvar(ng, 1, 'full');
Qg = sdpvar(ng, 1, 'full');

% 为节点电压赋初值
% -------- assign initial value -------- %
% Matpower solution
% mpopt = mpoption( 'out.all', 0);
% result = runpf(mpc, mpopt);
% assign(Vm_2, result.bus(:,VM).^2);
% assign(Va, result.bus(:,VA)*pi/180);
% Sbus = makeSbus(result.baseMVA, result.bus, result.gen);
% assign(Pbus, real(Sbus));
% assign(Qbus, imag(Sbus));
% assign(Pg, result.gen(:,PG));
% assign(Qg, result.gen(:,QG));

%% 目标函数
obj = C2' * (Pg.^2) + C1' * Pg + ones(ng, 1)' * C0;         % linear cost function ($/h)

%% 约束条件
Constraints = [];

%% A 节点功率平衡等式约束
Constraints = [Constraints; Pbus + Pd - Cg * Pg == 0];
Constraints = [Constraints; Qbus + Qd - Cg * Qg == 0];

%% B 节点注入与节点电压的线性等式约束（B,C二选一）
Constraints=[Constraints;Pbus == Xp_dlpf * [Va;Vm]];
Constraints=[Constraints;Qbus == Xq_dlpf * [Va;Vm]];

%% C 节点注入与支路功率的线性等式约束（B,C二选一）
% Constraints=[Constraints,Pbus == Cft' * Pf + sum(Gbus,2).*Vm ];
% Constraints=[Constraints,Qbus == Cft' * Qf - sum(Bbus,2).*Vm ];
% Constraints=[Constraints,Pbus == Ctf' * Pt + sum(Gbus,2).*Vm ];
% Constraints=[Constraints,Qbus == Ctf' * Qt - sum(Bbus,2).*Vm ];

%% D 支路功率与节点电压的等式约束（D1,D2二选一,注意Vm或Vm_2要与B、C一致）
% D1 (Vi - Vj)
Constraints=[Constraints; Pf == (Vm(f)-Vm(t)).* real(Ysc) - (Va(f)-Va(t)).*imag(Ysc) ];
Constraints=[Constraints; Qf == -(Vm(f)-Vm(t)).* imag(Ysc) - (Va(f)-Va(t)).*real(Ysc)];
Constraints=[Constraints; Pt == (Vm(t)-Vm(f)).* real(Ysc) - (Va(t)-Va(f)).*imag(Ysc)];
Constraints=[Constraints; Qt == -(Vm(t)-Vm(f)).* imag(Ysc) - (Va(t)-Va(f)).*real(Ysc)];

% D2 (Vi.^2 - Vj.^2)/2
% Constraints=[Constraints; -Pf  == (Vm_2(f)-Vm_2(t)).* real(Ysc) / 2 - (Va(f)-Va(t)).*imag(Ysc)];
% Constraints=[Constraints; -Qf  == -(Vm_2(f)-Vm_2(t)).* imag(Ysc)/ 2 - (Va(f)-Va(t)).*real(Ysc)];
% Constraints=[Constraints; -Pt  == (Vm_2(t)-Vm_2(f)).* real(Ysc)/ 2 - (Va(t)-Va(f)).*imag(Ysc)];
% Constraints=[Constraints; -Qt  == -(Vm_2(t)-Vm_2(f)).* imag(Ysc)/ 2 - (Va(t)-Va(f)).*real(Ysc)];

%% 电压上下限不等式约束（Vm或Vm_2与BCD保持一致）
Constraints = [Constraints; Vmin <= Vm <= Vmax];         % voltage magnitude (p.u.)
%Constraints=[Constraints, Vmin.^2 <= Vm_2 <= Vmax.^2];

%% 机组出力上下限不等式约束
Constraints = [Constraints; real(Sgmin) <= Pg <= real(Sgmax)];                                  % complex power output of generator (p.u.)
Constraints = [Constraints; imag(Sgmin) <= Qg <= imag(Sgmax)];

%% 平衡节点电压相角为0
Constraints = [Constraints; Va(ref) == 0];                                  % reference bus angle (rad)

%% 线路容量不等式约束（配电网不需要此约束）
t_poly=10;
if (nb < nl)
    Constraints = [Constraints,  Pf.^2 + Qf.^2 <= Fmax.^2];                  % branch flow limit at "from" end (p.u.)
    Constraints = [Constraints,  Pt.^2 + Qt.^2 <= Fmax.^2];
%     for i=1:t_poly
%         theta_poly = pi * (i-1) / (2*t_poly);
%         Constraints = [Constraints,(-Fmax<= cos(theta_poly) * Pf + sin(theta_poly) * Qf <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * Pf + cos(theta_poly) * Qf <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= cos(theta_poly) * Pt + sin(theta_poly) * Qt <= Fmax)];
%         Constraints = [Constraints,(-Fmax<= -sin(theta_poly) * Pt + cos(theta_poly) * Qt <= Fmax)];
%     end
else
end

%% 求解问题
opt = sdpsettings('solver','gurobi','verbose',0);% ipopt gurobi cplex
varargout=optimize(Constraints,obj,opt);
varargout.info

%% 检查约束违反
check(Constraints);

%% 输出结果
% -------- post-processing -------- %
%res.Ploss = sum(value(Ploss)) / baseMVA;
res.Vm = sqrt(value(Vm_2));
res.Va = value(Va)/pi*180;
res.Pbus = value(Pbus);
res.Qbus = value(Qbus);
res.Pf = value(Pf)* baseMVA;
res.Qf = value(Qf)* baseMVA;
res.Pt = value(Pt)* baseMVA;
res.Qt = value(Qt)* baseMVA;
res.Pg = value(Pg)* baseMVA;         % gen real power output (MW)
res.Qg = value(Qg)* baseMVA;
res.cost = value(obj);      % gen real power output cost ($/h)
res.lmp_energy = dual(Constraints(1));   % lmp($/MWh)
% res.lmp_con = dual(Const('conlow'));   % lmp($/MWh)
% res.muh_con = double(dual(Const('conhigh'))) .* (abs(double(dual(Const('conhigh'))))>1e-6);   % lmp($/MWh)
% res.mul_con = double(dual(Const('conlow'))) .* (abs(double(dual(Const('conlow'))))>1e-6);   % lmp($/MWh)

mpc.bus(:,VM) = res.Vm;
mpc.bus(:,VA) = res.Va;
mpc.gen(:,PG) = res.Pg;
mpc.gen(:,QG) = res.Qg;
mpc.branch(:,PF) = res.Pf;
mpc.branch(:,PT) = res.Pt;
mpc.branch(:,QF) = res.Qf;
mpc.branch(:,QT) = res.Qt;
V = mpc.bus(:,VM) .* exp(sqrt(-1) * pi/180 * mpc.bus(:,VA));
[mpc.bus, mpc.gen, mpc.branch] = pfsoln(mpc.baseMVA, mpc.bus, mpc.gen, mpc.branch, Ybus, Yf, Yt, V, ref, pv, pq, mpopt);
[i2eac, mpc.bus, mpc.gen, mpc.branch] = ext2int(mpc.bus, mpc.gen, mpc.branch);
success = 1;
results.success = success;
results.et = etime(clock, t0);
[mpc.bus, mpc.gen, mpc.branch, mpc.f] = deal(mpc.bus, mpc.gen, mpc.branch, res.cost);
results = int2ext(mpc);

%% zero out result fields of out-of-service gens & branches
if ~isempty(results.order.gen.status.off)
    results.gen(results.order.gen.status.off, [PG QG]) = 0;
end
if ~isempty(results.order.branch.status.off)
    results.branch(results.order.branch.status.off, [PF QF PT QT]) = 0;
end

if nargout == 1 || nargout == 2
    MVAbase = results;
    bus = success;
elseif nargout > 2
    [MVAbase, bus, gen, branch, f, et] = ...
        deal(results.baseMVA, results.bus, results.gen, results.branch, results.f, results.et);
    % else  %% don't define MVAbase, so it doesn't print anything
end
disp('done');

% %% 绘制功率圆
% % 参数设置
% R = 1; % 圆的半径
% num_squares = t_poly; % 旋转正方形的数量
% 
% % 圆的边界点
% theta_circle = linspace(0, 2 * pi, 100);
% x_circle = R * cos(theta_circle);
% y_circle = R * sin(theta_circle);
% 
% % 创建图形窗口
% figure;
% plot(x_circle, y_circle, 'b-', 'LineWidth', 1.5); % 绘制圆的边界
% hold on;
% 
% % 绘制多个旋转正方形
% for i = 1:num_squares
%     % 计算旋转角度
%     theta_ro = pi * (i-1) / (2*num_squares)
% 
%     % 定义正方形的四个顶点，未旋转时位于 x 和 y 轴方向
%     square_x = R * [-1, 1, 1, -1, -1];
%     square_y = R * [-1, -1, 1, 1, -1];
% 
%     % 旋转正方形
%     x_rotated = cos(theta_ro) * square_x - sin(theta_ro) * square_y;
%     y_rotated = sin(theta_ro) * square_x + cos(theta_ro) * square_y;
% 
%     % 绘制旋转后的正方形
%     plot(x_rotated, y_rotated, 'r-', 'LineWidth', 1.2);
% end
% 
% % 设置图形属性
% axis equal;
% xlim([-1.5, 1.5]);
% ylim([-1.5, 1.5]);
% xlabel('Rx');
% ylabel('Ry');
% title('Polygonal Approximation of Circle with Rotated Squares');
% legend('x^2 + y^2 = R^2', 'Rotated squares');
% grid on;
% hold off;
end