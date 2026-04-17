%% Test ACOPF Models
%  12-05-2022

clc
clear all

%% 定义选项
mpopt = mpoption( 'out.all', 0);%, 'opf.ac.solver', 'ipopt'
% mpopt = mpoption;

%% 导入系统数据
case_name = {'T-1';'T-2';'T-3';'T-4';'T-5';'T-6';'T-7';'T-8';'T-9'};

mpc = {loadcase('trans_opf_case9_ieee'); loadcase('trans_opf_case14_ieee'); loadcase('trans_opf_case_rts79'); ...
    loadcase('trans_opf_case30_ieee'); loadcase('trans_opf_case_rts96');loadcase('trans_opf_case57_ieee');...
    loadcase('trans_opf_case39_ieee');loadcase('trans_opf_case118_ieee');loadcase('trans_opf_case300_ieee')};

%% 不同模型求解
for c = 2
    % ACOPF
    Res_BaseAC{c} = runopf(mpc{c}, mpopt);
    
    % DCOPF
    %Res_BaseDC{c} = rundcopf(mpc{c}, mpopt);
    
    % DCOPF
    Res_LDC{c} = Lossless_dcopf(mpc{c}, mpopt);
    
    % Lossless Linear ACOPF Model
    mpc{c} = ext2int(mpc{c});
    tuni_less1 = clock;
    unit_less1 = cputime;
    Res_LL{c} = Lossless_Linear_opf(mpc{c}, mpopt);
    tuni_less2 = clock;
    unit_less2 = cputime;
    tuni_less = etime(tuni_less2,tuni_less1);
    unit_less = unit_less2 - unit_less1;
    
end