function K = fun_truss1d_K(E,alpha,L,A,p,theta)
%FUN_TRUSS1D_K
%    K = FUN_TRUSS1D_K(E,ALPHA,L,A,P,THETA)

%    This function was generated by the Symbolic Math Toolbox version 8.7.
%    14-Jul-2022 09:51:21

t2 = 1.0./L;
t3 = A.*E.*t2;
t4 = -t3;
K = reshape([t3,t4,t4,t3],[2,2]);
