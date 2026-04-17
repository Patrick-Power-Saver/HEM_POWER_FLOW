% This program solves the powerflow problem by Holomorphic embedding method.
%It is compared with powerflow  in rectangular coordinates. It accepts the data in matpower format.
clear;
clc;
% matpower data input change the file name as necessary
%mpc = trans_opf_case14_ieee;%
mpc = trans_opf_case118_ieee;
%mpc = trans_opf_case793_goc;
%mpc = trans_opf_case_rts96;

%mpc = dist_opf_case12_ieee;
%mpc = dist_opf_case69_ieee;
%mpc = dist_opf_case141_ieee;

% solves the power floe for tolerance of 1e-8
tol=1e-8;

%% 导入系统数据
define_constants;
mpopt = mpoption;
[baseMVA, bus, gen, branch, gencost] = loadcase(mpc);
[Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch);
[ref, pv, pq] = bustypes(bus, gen);
%[Sbus,Ybus,Yf,Yt,ref,pv,pq]=datapf(baseMVA, bus, gen, branch);
% rk=find(gen(:,1)==ref);
% gen(rk,:)=[];
tic;
%Sbus=Sbus.';
% power flow solved in rectangular coordinates.
%[V_base,cS,Sbus,LFs,LFr,pv1,pq1]=rectpowfl (baseMVA,Sbus,Ybus,Yf,Yt,bus,branch,gen,ref,pv,pq);
opt = mpoption('verbose', 0, 'out.all', 0);
tic
result = runpf(mpc, opt);
Vm = result.bus(:, 8);  % 第8列是 VM (Voltage Magnitude)
Va = result.bus(:, 9);  % 第9列是 VA (Voltage Angle)
Va_rad = Va * (pi / 180);
V = Vm .* exp(1i * Va_rad);
toc
Sbus=V.*conj(Ybus*V);
Sbus(pv) = real(Sbus(pv));
Sbus(ref)=0;
V0=bus(:,8);

% VV1-Bus votages determined from Holomorphic mbedding method.HEM方法的电压
%Vh-Polynomial coefficients os 's' determined. 多项式系数
%k -iteration 迭代次数
%delS-accuracy of the solution.误差
tic
[VV1,Vh,Qh,k,delS] = HEPF(Ybus, Sbus, V0, ref, pv, pq, tol);
toc

% ---  Padé逼近方法 ---
fprintf('正在计算 Padé 加速解...\n');
VV_pade = zeros(size(VV1));

% 核心修正点：确保索引不超出 Vh 的实际列数
[n_bus, max_cols] = size(Vh);
actual_k = min(k, max_cols); 

for i = 1:n_bus
    % 提取第 i 个节点的有效级数系数
    coeffs = Vh(i, 1:actual_k);
    
    % 只有当系数不全为零时才进行 Padé 逼近
    if any(coeffs)
        VV_pade(i) = pade_approx(coeffs);
    else
        VV_pade(i) = VV1(i); % 若无系数则保留原值
    end
end

% --- 结果对比 ---
fprintf('V_Matpower    V_HEM(直接求和)    V_HEM(Padé加速)\n')
disp([V(1:5), VV1(1:5), VV_pade(1:5)]); % 仅展示前5个节点

err_normal = max(abs(V - VV1));
err_pade = max(abs(V - VV_pade));
fprintf('普通 HEM 最大误差: %.2e\n', err_normal);
fprintf('Padé 加速后最大误差: %.2e\n', err_pade);