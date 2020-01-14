%% Heat calculation 
% shane 2019/11/27

%% …………………………华丽的分隔线…………………………%%
% This method can be consist of three models:
%        NO.1:calculate the heat generation
%        NO.2:find the heat resistence
%        NO.3:sloving the equations

clc
clear 
close all
global R_heat Q T0 T8
%% bearing heat calculation
%% Improving Palmgren's Way
% global way of calculating the heat generation
n=20000;     %rotational sped:RPM
v=90.7;     %Dynamic viscosity:Pa*s
f_0=0.9;    %ecofficient relate to lubricant and bearing type
dm=53.90;    %pitch diameter
Zn=18;       %number of balls 
temp=v*n;
if temp>=2000
    M_v=1e-7*f_0*temp^(2/3)*dm^3;  %  N*mm
else
    M_v=160e-7*f_0*dm^3;
end
% silding fraction 
% F=[5.68911385499705;10.0636041623618;14.4308602079465;18.7341404608281;22.9705437451725;27.1444771391662;31.2611549675728;35.3253844428405;39.3413730961347;41.6929570868337;47.2427202865732;51.1339978757883;54.9890205447014;58.8099319964538;62.5986413787958;66.3568594316891;70.0861277701651;73.7878427713478;77.4632752099748;81.1135865257848;84.7398424038546;88.3430241948852;91.9240385856011;95.4837258407382;99.0228668704729;102.542189325249;106.042372879891;109.524053837694;112.987828793780;116.434260013360];
% mu=[0.0238481119503919,0.00938684354080753,0.00596488604487738,0.00441155840498280,0.00349691577832388,0.00289624914760768,0.00248122571663885,0.00217760137478595,0.00194655329792208,0.00183470884311368,0.00162118814121175,0.00150430510177064,0.00140911584017210,0.00133158137623968,0.00126889042944305,0.00121907755859035,0.00118078239391033,0.00115309323705877,0.00113544306687890,0.00112753920828962,0.00112931529095026,0.00114089838300126,0.00116258673649697,0.00119483514614172,0.00123824591619295,0.00129356407214605,0.00136167811317730,0.00144361285595259,0.00154054404825442,0.00165379743786914];
% for j=1:6
% DT=[];
% for i=1:length(mu)
% P=23;
% mu_mix=0.0103;
P_mu=[0.885431624;
0.813211808;
0.839611691;
0.863582547;
0.905419543;
0.923767893;
0.94066572;
0.751754447;
0.769257536;
0.814566881;
0.94992498;
1.008961057;
1.063260335];
DT=[];
for i=1:length(P_mu)
M_s=0.5*dm*P_mu(i)*Zn;
H_tot=1.01e-5*(M_v*1+M_s)*n;     % Heat power: W

%% Bearing heat generation
% shane's way (we only consider the heat generate by
%      sild fraction and air-oil drag force)



% Heat distribution
Q1=0.5*H_tot;
Q2=0.25*H_tot;
Q3=Q2;

%% Heat resistance calculation
% Contains 3 convective resistances,4 conduction...
%    resistances and 2 unknowable resistances.

% convective resistances
alpha=[177.27 201.82 253.61]; % convection coefficient
D_b=7.3e-3;
DO_o=68e-3;
DO_i=58.8e-3;
DI_o=50.95e-3;
DI_i=40e-3;
B=15e-3;
R_03=1/(pi*B*alpha(2)*DO_i);
R_04=1/(pi*alpha(3)*D_b^2);
R_05=1/(pi*B*alpha(1)*DI_o);

%conduction resistances
K_b=60.5;       % conduction coefficient
R_34=1/(pi*K_b*D_b);
R_45=R_34;
R_23=log(DO_o/DO_i)/(2*pi*B*K_b);
R_56=log(DI_o/DI_i)/(2*pi*B*K_b);

% unknowable resistances
R_28=20;
R_68=20;


%% Parameters transition
R_heat=[R_34 R_04 R_45 R_03 R_23 R_05 R_56 R_28 ...
        R_68];
Q=[Q1 Q2 Q3];
T0=35;
T8=20;


%% Slove equations 
% Not special
T1=[10 12 13 14 25];%初值向量：5*1
options =optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[Y,fval,exitflag]=fsolve(@fun_heat,T1,options);
if exitflag==1
    disp('结果收敛');
else
    disp('所设初值无效');
end
DT(i,1)=Y(1)-T0;
DT(i,2)=Y(3)-T0;
end
%% …………………………华丽的分隔线…………………………%%
% disp('Suppose conditions:')
% sprintf('Environment Temperature:%0.2d C \nLubricant Temperature:%0.2d C \n',T8,T0)
% sprintf('Heat generation:%0.2d W \n',Q)
% sprintf('Heat resistance:%0.2d \n',R_heat)
% sprintf('Ball Temperature raise:%0.3f C \nOuter ring Temperature:%0.3f C \nInner ring Temperature raise:%0.3f C \n',Y(1)-T0,Y(2),Y(3)-T0)


