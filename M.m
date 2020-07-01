global sigma nu betta b alppha deltta pi S tau k_n 
addpath C:\dynare\4.6.1\matlab
!psexec

nu      =0.37;
betta   =0.996;
b       =0.7;
alppha  =0.3;
sigma   =2;
deltta  =0.025;
S       =1.00645;
tau     =1;
k_n     =2;
pi      =1.01;

options.Display = 'iter';

Mss0=ones(13,1);
[Mss,exitflag]=fsolve(@M_ss,Mss0,options);
exitflag;

dynare ME.mod noclearall parallel 

