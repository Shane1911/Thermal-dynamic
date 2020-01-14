%%    The dynamic model with consider the heat...
%           expansion of bearing elements 
%   shane 2019/11/26

%% …………………………华丽的分隔线…………………………%%
clc
clear
close all
global alfa A Zn F K_oj K_ij f D omega_i m I dm Q_1j...
    Q_2j alfa_i alfa_o Ic omega_Rj omega_cj dF
%% 参数输入
% f=[0.523 0.523];%沟道曲率
dm=53.9;        %节圆直径 单位：mm
D=7.3;          %钢球直径 单位：mm
alfa=15;        %初始接触角 /°度
Zn=18;          %钢球个数
den=7800;       %钢球密度 kg/mm^3
m=pi*den*(D*1e-3)^3/6;      %钢球质量 kg
I=m*(D*1e-3)^2/10;          %钢球的转动惯量 kg*m^2
Ic=I+0.25*m*(dm*1e-3)^2;    %钢球绕轴承轴线的转动惯量
K_ij =7.5806e+05;           %单位 N*mm^1.5
K_oj =7.7985e+05;           %单位 N*mm^1.5
% A=(f(1)+f(2)-1)*D;          %曲率中心距离 mm

%% 外部输入
rate=[];
omega_i=15000*pi/30;%内圈转速 rad/s
T=2*pi/omega_i;%内圈转动周期 s
omega_cj=0.5*omega_i*(1-D*cosd(alfa)/dm);
alfa_j=atan(sind(alfa)/(cosd(alfa)+D/dm));
omega_Rj=0.5*omega_i*dm/D*(1-(D*cosd(alfa)/dm)^2)*sin(alfa_j);
W1=3*[450 1600 -250 227.9*0.9 2*pi/Zn 82*0.9 0];%角速度初值向量：7*1
W=W1;
output=[];le=0;
fin_t=0.2/0.2/T;
W_x=[];W_y=[];W_z=[];W_c=[];alfa_oi=[];
Q12=[];Deta_a=[];HC=[];MU=[];
for dF=1:2;
F=100;
% fi=[0.528 0.533 0.538];
f=[0.523 0.523 ];%沟道曲率
A=(f(1)+f(2)-1)*D;   
for i=1:fin_t
% 初值设置
deta_a=0.0103;
X=[0.0447 0.1622 0.0024 0.0018 deta_a];%初值向量：5*1
options =optimoptions('fsolve','Algorithm','Levenberg-Marquardt');
[Y,fval,exitflag]=fsolve(@fun,X,options);
if exitflag==1
    disp('结果收敛');
else
    disp('所设初值无效');
end
Q_1j=K_oj*Y(3)^1.5;%单位/N
Q_2j=K_ij*Y(4)^1.5;%单位/N
A_2j=A*cosd(alfa);
cos_1j=Y(2)/((f(1)-0.5)*D+Y(3));
cos_2j=(A_2j-Y(2))/((f(2)-0.5)*D+Y(4));
alfa_i=acos(cos_2j);
alfa_o=acos(cos_1j);
tspan = 0:1e-5:0.2*T;
[t,y] = ode15s(@fun1,tspan,W,options);
num=length(t);
W=y(num,:);%角速度初值向量：7*1
omega_cj=y(num,4);
omega_Rj=y(num,1);
output(1+le:le+num,:)=y(1:num,:);
le=le+num;
end
temp=length(output(:,1));
W_x(dF)=output(temp,1);
W_y(dF)=output(temp,2);
W_z(dF)=output(temp,3);
W_c(dF)=output(temp,4);
alfa_oi(dF,1)=alfa_o;
alfa_oi(dF,2)=alfa_i;
Y_fin=output(temp,:);
suboutput=fun4(Y_fin);
HC(dF,:)=suboutput(1:2);
MU(dF,1)=suboutput(3);
MU(dF,2)=suboutput(4);
Q12(dF,1)=Q_1j;
Q12(dF,2)=Q_2j;
Deta_a(dF)=Y(5);
dF=dF+1;
end
%%
num2=length(W_x);
uox=[];uix=[];
for k=1:num2
  alfao=alfa_oi(k,1);
  alfai=alfa_oi(k,2);
  W=[W_x(k) W_y(k) W_z(k) W_c(k)];
  uox(k)=0.5*W(4)*1e-3*(dm)+0.5*D*1e-3*(W(3)*cos(alfao)-W(2)*sin(alfao)+W(4)*cos(alfao));
  uix(k)=0.5*1e-3*(-W(4)+omega_i)*(dm)+0.5*D*1e-3*(W(3)*cos(alfai)-W(2)*sin(alfai)-(-W(4)+omega_i)*cos(alfai));
end
F_y=20:20:600;
plot(F_y,uox,'-*');
hold on
plot(F_y,uix,'-S');
grid on
% num1=size(W_x,1);
% num2=size(W_x,2);
% wso=[];wsi=[];
% for k=1:num2
%   alfa_o=alfa_oi(1,k);
%   alfa_i=alfa_oi(2,k);
%   W=[W_x(num1,k) W_y(num1,k) W_z(num1,k) W_c(num1,k)];
%   wso(k)=W(3)*sin(alfa_o)+W(2)*cos(alfa_o)+W(4)*sin(alfa_o);
%   wsi(k)=W(3)*sin(alfa_i)+W(2)*cos(alfa_i)+(W(4)-omega_i)*sin(alfa_i);
% end
% F_y=20:20:600;
% plot(F_y,-wso,'b-*');
% hold on
% plot(F_y,-wsi,'r-*')
% grid on
% num2=length(W_x);
% uox=[];uix=[];
% for k=1:num2
%   alfao=alfa_oi(k,1);
%   alfai=alfa_oi(k,2);
%   W=[W_x(k) W_y(k) W_z(k) W_c(k)];
%   uox(k)=W(3)*sin(alfa_o)+W(2)*cos(alfa_o)+W(4)*sin(alfa_o);
%   uix(k)=W(3)*sin(alfa_i)+W(2)*cos(alfa_i)+(W(4)-omega_i)*sin(alfa_i);
% end
% F_y=60:20:200;
% plot(F_y,uox,'-*');
% hold on
% plot(F_y,uix,'-S');
% grid on
