function [t,air_infot]=check_err(air_info,target)
% if err_hor and err_ver all out of bound, then t = 1
% else t = 0
% air_info1: air infomation at target node
% air_info: air infomation at origin node
% target: target nodeNo
global gf parmbnd delta theta
% cal distance
d=norm(gf(air_info(1),2:4)-gf(target,2:4),2);
% cal err
err=d*delta;

air_infot=[target,air_info(2:3)+err,air_info(4)+d,air_info(5)];
t=(air_infot(2)>theta)|(air_infot(3)>theta);
% correction
if gf(target,5)==0 && air_infot(2)<=parmbnd(4) && air_infot(3)<=parmbnd(3) % horizon
    air_infot(2)=0;
    air_infot(5)=air_infot(5)+1;
elseif gf(target,5)==1 && air_infot(2)<=parmbnd(2) && air_infot(3)<=parmbnd(1) % vertical
    air_infot(3)=0;
    air_infot(5)=air_infot(5)+1;
end
end