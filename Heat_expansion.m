clc
close all
clear 
% …………………………华丽的分隔线………………………… %
Tep=[20 15;
    30 25;
    60 55];
Tep_b=Tep(:,1);
Tep_r=Tep(:,2);
r_b=3.65;
R_i=29.65;
alpha_b=3.2e-6;
alpha_r=12.5e-6;
expansion_b=alpha_b*r_b.*Tep_b;
sprintf('滚子热膨胀量为: %0.2d mm',expansion_b)
expansion_r=alpha_r*R_i.*Tep_r;
sprintf('内圈热膨胀量为: %0.2d mm',expansion_r)
y=[];
for i=1:length(Tep_b)
    y(i,:)=[expansion_b(i) expansion_r(i)];
end