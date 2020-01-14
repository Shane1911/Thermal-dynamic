function [ Y ] = fun( X );
global alfa A Zn F K_oj K_ij f D m I dm...
omega_cj omega_Rj dF
    Expansion=[0 0;
0.0003504	0.009265625;];      
    eb=Expansion(dF,1);
    er=Expansion(dF,2);
    Pd=0;
%     eb=0;
%     er=0;
    A_1j=A*sind(alfa)+X(5);
    A_2j=A*cosd(alfa)+er-Pd;
    cos_1j=X(2)/((f(1)-0.5)*D+X(3)-eb);
    sin_1j=X(1)/((f(1)-0.5)*D+X(3)-eb);
    cos_2j=(A_2j-X(2))/((f(2)-0.5)*D+X(4)-eb);
    sin_2j=(A_1j-X(1))/((f(2)-0.5)*D+X(4)-eb);
    Q_1j=K_oj*X(3)^1.5;%单位/N
    Q_2j=K_ij*X(4)^1.5;%单位/N
    F_cj=0.5*dm*1e-3*m*omega_cj^2;
    M_gj=I*omega_Rj*omega_cj;
    F_1j=M_gj*1e3/D;
    F_2j=M_gj*1e3/D;
    Y(1,1)=-Q_1j*sin_1j+Q_2j*sin_2j-F_1j*cos_1j+F_2j*cos_2j;
    Y(2,1)=-Q_1j*cos_1j+Q_2j*cos_2j+F_1j*sin_1j-F_2j*sin_2j+F_cj;
    Y(3,1)=(A_1j-X(1))^2+(A_2j-X(2))^2-((f(2)-0.5)*D+X(4)-eb)^2;
    Y(4,1)=X(1)^2+X(2)^2-((f(1)-0.5)*D+X(3)-eb)^2;
    Y(5,1)=F-Zn*(Q_2j*sin_2j+F_2j*cos_2j);
end
