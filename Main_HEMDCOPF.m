%% Test HEM-DCOPF Models

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
    Res_BaseDC{c} = rundcopf(mpc{c}, mpopt);
    %Res_LDC{c} = Lossless_dcopf(mpc{c}, mpopt);
    
    %HEM-DCOPF
    mpc{c} = ext2int(mpc{c});
    Res_HEMDC{c} = HEM_DCOPF(mpc{c},10);
    
end