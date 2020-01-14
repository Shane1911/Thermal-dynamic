function [ Y ] = fun4(W)
global  Q_1j Q_2j alfa_i alfa_o dm D omega_i alfa hc 
cos_1j=cos(alfa_o);sin_1j=sin(alfa_o);
cos_2j=cos(alfa_i);sin_2j=sin(alfa_i);
%求hc
Q_12=[Q_1j Q_2j];
u=[];
u(1)=0.5*W(4)*1e-3*(dm+D*cos(alfa_o))+0.5*D*1e-3*(W(2)*sin(alfa_o)-W(3)*cos(alfa_o));
u(2)=0.5*1e-3*(-W(4)+omega_i)*(dm-D*cos(alfa_i))+0.5*D*1e-3*(W(2)*sin(alfa_i)-W(3)*cos(alfa_i));
Rx=[];
Rx(1)=0.5*D*1e-3*(1+D*cosd(alfa)/dm);
Rx(2)=0.5*D*1e-3*(1-D*cosd(alfa)/dm);
E1=2.075e11;
E=E1/(1-0.3^2);
yita_0=0.04667;
alfa_0=1.2*1e-8;
U=yita_0.*u./(2*E.*Rx);
G=E*alfa_0;
W_l=Q_12./(E.*Rx.^2);
k=[6.9731 8.2435];
hc=2.69*Rx.*U.^(0.67).*G.^(0.53).*W_l.^(-0.067).*(1-0.61*exp(-0.73.*k));%单位：m
% V=0.026;deta=6.3e-5;
% Hc=3.672*W_l.^(-0.045*k.^0.18).*U.^(0.663*k.^0.025).*G.^(0.502*k.^0.064)...
%     .*(1-0.572*exp(-0.74.*k)).*(1+0.025*deta^1.248.*V^0.119.*W_l.^(-0.133).*U.^(-0.884)...
%     .*G.^(-0.977).*k.^0.081);
% hc=Hc.*Rx;
%计算滚动阻力
deta1=[1.0269 1.0174];
% kr=0.5968./(deta1-1.0003);
% F_R=2.86*E.*Rx.^2.*kr.^(0.348).*G.^(0.022).*U.^(0.66).*W_l.^(0.47);
%求摩擦力和摩擦力矩
%求椭圆参数
R1=1e-3*[3.9320 3.0557];;%m
a=(6*k.^2.*deta1.*Q_12.*R1/(pi*E)).^(1/3);%椭圆长轴 m
b=a./k;%椭圆短半轴 m
W_j=W(1:4);
ymax1=@(x)a(1)*sqrt(1-(x/b(1)).^2);
ymax2=@(x)a(2)*sqrt(1-(x/b(2)).^2);
fox=4/hc(1)*integral2(@(x,y)fun2(x,y,W_j,a,b,1),0,b(1),0,ymax1);
c1=fox/Q_1j;
foy=4/hc(1)*integral2(@(x,y)fun2(x,y,W_j,a,b,2),0,b(1),0,ymax1);
fix=4/hc(2)*integral2(@(x,y)fun2(x,y,W_j,a,b,3),0,b(2),0,ymax2);
c2=fix/Q_2j;
fiy=4/hc(2)*integral2(@(x,y)fun2(x,y,W_j,a,b,4),0,b(2),0,ymax2);
Y=[hc c1 c2];
% Y=[hc fix fiy];
end

