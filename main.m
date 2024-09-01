clc;clear;
fun = @object;
nonlcon = @nonlcon;
TenBarTruss = @TenBarTruss;

% 初始點&線性限制不等式
r0  =  [0.1,0.1];
A = [];
b = [];
% 線性限制等式
Aeq = [];
beq = [];
% 目標值上下限 & 運算方式
lb = [0.001,0.001];
ub = [0.5,0.5];
options = optimoptions('fmincon','Display','final','Algorithm','sqp');

% 最佳化運算 & 結果顯示
[r,fval,exitflag] = fmincon(fun,r0,A,b,Aeq,beq,lb,ub,nonlcon,options);
[Q, stress,R,K] = TenBarTruss(r)