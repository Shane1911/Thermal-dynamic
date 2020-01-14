function [ Y ] = fun1(t, W )
global  Q_1j Q_2j alfa_i alfa_o dm D omega_i alfa hc I Zn Ic 
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
R1=1e-3*[3.9320 3.0557];%m
a=(6*k.^2.*deta1.*Q_12.*R1/(pi*E)).^(1/3);%椭圆长轴 m
b=a./k;%椭圆短半轴 m
W_j=W(1:4);
ymax1=@(x)a(1)*sqrt(1-(x/b(1)).^2);
ymax2=@(x)a(2)*sqrt(1-(x/b(2)).^2);
fox=4/hc(1)*integral2(@(x,y)fun2(x,y,W_j,a,b,1),0,b(1),0,ymax1);%F_R(1);
foy=4/hc(1)*integral2(@(x,y)fun2(x,y,W_j,a,b,2),0,b(1),0,ymax1);
fix=4/hc(2)*integral2(@(x,y)fun2(x,y,W_j,a,b,3),0,b(2),0,ymax2);%-F_R(2);
fiy=4/hc(2)*integral2(@(x,y)fun2(x,y,W_j,a,b,4),0,b(2),0,ymax2);
mo=4/hc(1)*integral2(@(x,y)fun2(x,y,W_j,a,b,5),0,b(1),0,ymax1);
mi=4/hc(2)*integral2(@(x,y)fun2(x,y,W_j,a,b,6),0,b(2),0,ymax2);
Mx=0.5*1e-3*D*(foy-fiy);
My=0.5*1e-3*D*(fix*sin(alfa_i)+fox*sin(alfa_o))-(mo*cos(alfa_o)+mi*cos(alfa_i));
Mz=0.5*1e-3*D*(-fix*cos(alfa_i)-fox*cos(alfa_o))-(mo*sin(alfa_o)+mi*sin(alfa_i));
f_cage=0; 
WC=0.5*omega_i*(1-D*cosd(alfa)/dm);
f_v=(0.05*pi*860*(D*1e-3)^2*(dm*1e-3*WC)^(1.95))/32/9.8;
V1=0.5*52*1e-3*(omega_i-W(4));
FM_cage=2*pi*yita_0*1e-3*(0.5*52)^2*V1*8/0.25/Zn/(0.5*dm);
Tao=diag([I I I Ic 1 Ic 1]);
MF=[];
MF(1)=Mx+I*W(4)*W(2);
MF(2)=My-I*W(4)*W(1);
MF(3)=Mz-0.02*FM_cage*0.5*D*1e-3;
MF(4)=0.5*1e-3*(-46*fox+31*fix)-0.5*1e-3*dm*(f_cage+f_v+FM_cage);
MF(5)=W(4);
MF(6)=-0.5*1e-3*(-70.1*fox+43.1*fix)+0.5*1e-3*dm*(f_v);
MF(7)=W(6);
Y=Tao\MF';
end


