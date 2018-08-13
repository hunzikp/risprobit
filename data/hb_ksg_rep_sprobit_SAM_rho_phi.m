%HB_KSG Replication PSRM 
%Franzese, Hays, Cook 2014
clear all;
tic; 

r=100;

W=csvread('C:\Users\scott\Dropbox\Spatial Probit PSRM\PSRM Conditional Accept\Replication Folder\Application - Civil War\Data\hb_ksg_contig_W.csv');%non-rs
WNT=normw(W);


D=csvread('C:\Users\scott\Dropbox\Spatial Probit PSRM\PSRM Conditional Accept\Replication Folder\Application - Civil War\Data\hb_ksg_rep_Africa.csv');

Y=D(:,4);%Civil War Incidence 

nt = length(Y);  

XNT=D(:,[8:15]);%neighpol neighpolsq neighlgdpl lnpop polity2l polity2sq postcoldw lgdp96l

k=size(XNT,2)+1;

pyrs = D(:,5); 
pyrs2 = pyrs.^2; 
pyrs3 = pyrs.^3; 

py = [pyrs pyrs2 pyrs3];

TL=csvread('C:\Users\scott\Dropbox\Spatial Probit PSRM\PSRM Conditional Accept\Replication Folder\Application - Civil War\Data\hb_ksg_TL.csv');

i_nt = eye(nt); 

rand_mat1 = rand(nt,(r/2));
rand_mat2 = (1 - rand_mat1);
rand_mat = [rand_mat1, rand_mat2];


z = 1-2*Y;
Z = diag(z);


ytl = TL*Y; 

niave_ysl = WNT*Y; 

%Niave Probit (RS)
cons = ones(nt,1); 
[b,dev,stats] = glmfit([cons(31:nt,:) XNT(31:nt,:) niave_ysl(31:nt,:) ytl(31:nt,:)],Y(31:nt,:),'binomial','link','probit','constant','off');
beta_rs = stats.beta;
se_rs = stats.se; 
temp_rs = [beta_rs, se_rs, beta_rs./se_rs] ;
printmat(temp_rs,'Naive Results - RS','constant neighpol neighpolsq neighlgdpl lnpop polity2l polity2sq postcoldw lgdp96l spat_lag timelag', 'beta se t-stat')


%Niave Probit (DTH)
cons = ones(nt,1); 
[b,dev,stats] = glmfit([cons(31:nt,:) XNT(31:nt,:) niave_ysl(31:nt,:) py(31:nt,:)],Y(31:nt,:),'binomial','link','probit','constant','off');
beta_dth = stats.beta;
se_dth = stats.se; 
temp_dth = [beta_dth, se_dth, beta_dth./se_dth]; 
printmat(temp_dth,'Naive Results - DTH','constant neighpol neighpolsq neighlgdpl lnpop polity2l polity2sq postcoldw lgdp96l spat_lag py1 py2 py3', 'beta se t-stat')
 
p0 = zeros(1,k+2); 

lb = [-Inf; -Inf; -Inf; -Inf; -Inf; -Inf; -Inf; -Inf; -Inf; -0.85; -0.85];
ub = [Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; Inf; 0.85; 0.85]; 
A = [0 0 0 0 0 0 0 0 0  1 1;
     0 0 0 0 0 0 0 0 0 -1 -1];
b = [0.9; 0.9];

options = optimset;
% Modify options setting
options = optimset(options,'Display' ,'iter');
options = optimset(options,'UseParallel','always'); 
%options = optimset(options,'FunValCheck' ,'on');
options = optimset(options,'PlotFcns' ,{  @optimplotx @optimplotfval @optimplotfirstorderopt });
options = optimset(options,'Diagnostics' ,'on');
options = optimset(options,'HessUpdate' ,'bfgs'); %dfp
%options = optimset(options,'InitialHessType' ,'scaled-identity');
options = optimset(options,'Algorithm' ,'interior-point');
options = optimset(options,'LargeScale' ,'on');
options = optimset(options,'MaxFunEvals' ,10e4);
options = optimset(options,'MaxIter' ,10e2);
options = optimset(options,'TolX', 10e-8);
[p,fval,exitflag,output,lambda,grad,hessian] = ...
fmincon(@(p)sprobit_llf_hbksg_RR_rho_phi(p,Z,rand_mat,nt,WNT,XNT,TL,i_nt,r),p0,A,b,[],[],lb,ub,[],options);
se = sqrt(diag(inv(hessian)));
pt = p';
vcov_hat = inv(hessian); 

temp_stpmsl = [pt, se, pt./se] ;
printmat(temp_stpmsl,'STP-MSL','constant neighpol neighpolsq neighlgdpl lnpop polity2l polity2sq postcoldw lgdp96l timelag spat_lag', 'beta se t-stat')


%save('hb_ksg_RR_rho_phi_results');
