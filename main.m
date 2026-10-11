clc;
clear all;
format long;
warning on;
%%
gvar;
load('data2imuavp.mat');
%% 地球重力模型
C = zeros(11); S = C;
C(1:2:11,1) = [1;
    -4.84116e-4;   % WGS-84 系数
    7.90303e-7;
    -1.687249e-9;
    3.460524e-12;
    -2.650e-15];
earth.GM=GM; earth.Re=Re; earth.wie=wie; earth.ff=ff; earth.C=C; earth.S=S;
earth.e2=e2;
%% 统一格式 列状态
avp = avp';
imu = imu';
%%
n_sample = 2; %子样数 —— 原始对比为单子样+前一周期
n_tmp=1;
imu_buff = zeros(7,2);
Step = size(imu,2);

avp_buff = zeros(10,fix(Step/n_sample)+1); %俯仰、横滚、方位角，东-北-天速度，纬经高，时间t
avp_buff(:,1) = [avp0;0];
q = a2q(avp0(1:3,1));
q_before=q;
%% 外推缓存
wnie_before=zeros(3,1); wnen_before=zeros(3,1);
gn_before=zeros(3,1); Mpv_before = zeros(3,3);
%%
for tmp1 = 1:Step
    if mod(tmp1,n_sample) == 0
        t=tmp1/n_sample;
        %继承上一时刻姿态
        q_before=q;
        %传感器离散时间
        if t==1
            dT=imu(end,tmp1);
        else
            dT=imu(end,tmp1)-imu(end,tmp1-2);
        end
        %数据提出
        imu_buff(:,n_tmp) = imu(:,tmp1); %x\y\z陀螺角增量，x\y\z加计速度增量，t采样时间
        n_tmp=1;
        %地球数据解算
        earth = eath_param_cal(avp_buff(7:9,t),avp_buff(4:6,t),earth);
        %圆锥/桨划、旋转补偿后的运动增量
        [phi_b,dv_scul,dv_rot] = correct_rot_mo(imu_buff,n_sample);
        %姿态更新
        q=q_update(q_before,phi_b,earth.wnin,dT);
        %姿态四元数到欧拉角
        avp_buff(1:3,t+1)=q2a(q);
        % Mpv
        Mpv = [0 1/earth.RMh 0; sec(avp_buff(7,t))/earth.RNh 0 0;0 0 1];
        %外推
        if t==1
            v_mid=avp_buff(4:6,t); 
            wnie_mid=earth.wnie; wnen_mid=earth.wnen; 
            gn_mid=earth.gn;
            Mpv_mid = Mpv;
        else
            v_mid=linear_push(avp_buff(4:6,t-1),avp_buff(4:6,t)); 
            wnie_mid=linear_push(wnie_before,earth.wnie); wnen_mid=linear_push(wnen_before,earth.wnen); 
            gn_mid=linear_push(gn_before,earth.gn);
            Mpv_mid = linear_push(Mpv_before,Mpv);
        end
        wnie_before=earth.wnie; wnen_before=earth.wnen;
        gn_before=earth.gn;
        Mpv_before=Mpv;
        %速度更新
        avp_buff(4:6,t+1)=v_update(avp_buff(4:6,t),imu_buff(4:6,:),dv_scul,dv_rot,q_before,dT,v_mid,wnie_mid,wnen_mid,gn_mid);
        %位置更新
        avp_buff(7:9,t+1) = p_update(avp_buff(7:9,t),avp_buff(4:6,t),avp_buff(4:6,t+1),dT,Mpv_mid);
        %时间戳记录
        avp_buff(10,t+1)=imu(end,tmp1);
    else
        imu_buff(:,n_tmp) = imu(:,tmp1);
        n_tmp=n_tmp+1;
    end
end

avp_buff_results = avp_buff(:,2:end);
%%
plot_here((avp(:,n_sample:n_sample:end))', avp_buff_results);
%% 位置更新
function p_m = p_update(p_before,v_before,v_m,dT,Mpv_mid)
p_m=p_before+Mpv_mid*(v_before+v_m)*dT/2;
end
%% 速度更新
function v_m=v_update(v_before,dv_imu,dv_scul,dv_rot,q_before,dT,v_mid,wnie_mid,wnen_mid,gn_mid)
dv1=dv_imu(:,1); dv2=dv_imu(:,2); dv_m=sum(dv_imu,2);
tmp = q_rot(q_before,dv1+5*dv2);
wnin_mid=wnie_mid+wnen_mid;
dv_sf = q_rot(q_before,dv_m)-dT/6*cross(wnin_mid,tmp)+q_rot(q_before,dv_scul+dv_rot);
dv_cor=dT*(gn_mid-cross(2*wnie_mid+wnen_mid,v_mid));
v_m = v_before+dv_sf+dv_cor;
end
%% 姿态更新
function q=q_update(q_before,phi_b,wnin,dT)
phi_n=dT*wnin;
q_n=phi2q(phi_n); q_b=phi2q(phi_b);
q=q_multi(q_multi(q_inv(q_n),q_before),q_b);
end
%% 地球信息解算
function earth = eath_param_cal(pos,vn,earth) %pos: 纬、经(弧度)、高 vn: 东北天ENU速度 e2: 偏心率平方 ff: 扁率
slat = sin(pos(1)); clat = cos(pos(1)); tlat = tan(pos(1));
tmp = 1-earth.e2*slat^2;
RN = earth.Re/sqrt(tmp); RM = RN*(1-earth.e2)/tmp;
earth.RMh=RM+pos(3); earth.RNh=RN+pos(3);
earth.wnie=[0 earth.wie*clat earth.wie*slat]';
earth.wnen=[-vn(2)/earth.RMh vn(1)/earth.RNh vn(1)*tlat/earth.RNh]';
earth.wnin=earth.wnen+earth.wnie;
earth.gn=egm(earth.GM, earth.Re, earth.wie, earth.ff, earth.C, earth.S, pos(1), pos(2), pos(3));
end
%% 圆锥、桨划补偿
function [phi_b,dv_scul,dv_rot] = correct_rot_mo(imu_buff,n_sample)
dtheta = imu_buff(1:3,:); dv = imu_buff(4:6,:);
dtheta_m = sum(imu_buff(1:3,:),2); dv_m = sum(imu_buff(4:6,:),2);

dphi=zeros(3,1); dv_scul=zeros(3,1);
c_table = zeros(4,4); %补偿系数表格
c_table(2,1)=2/3; c_table(3,1:2)=[27/20 9/20]; c_table(4,1:3)=[214/105 92/105 54/105];
for tmp1=1:n_sample-1
    dphi=dphi+c_table(n_sample,n_sample-tmp1)*cross(dtheta(:,tmp1),dtheta(:,n_sample));
    dv_scul=dv_scul+c_table(n_sample,n_sample-tmp1)*cross(dtheta(:,tmp1),dv(:,n_sample))+...
        c_table(n_sample,n_sample-tmp1)*cross(dv(:,tmp1),dtheta(:,n_sample));
end
dv_rot = cross(dtheta_m,dv_m)/2;
phi_b = dtheta_m+dphi;
end

%% 欧拉角到四元数 
function q=a2q(a) %a: pitch roll yaw （弧度）
sa = sin(a/2); ca = cos(a/2);
q = [ca(3)*ca(1)*ca(2)-sa(3)*sa(1)*sa(2);
     ca(3)*sa(1)*ca(2)-sa(3)*ca(1)*sa(2);
     sa(3)*sa(1)*ca(2)-ca(3)*ca(1)*sa(2);
     sa(3)*ca(1)*ca(2)-ca(3)*sa(1)*sa(2);];
end
%% 四元数到欧拉角
function a=q2a(q) %a: pitch roll yaw （弧度）
tmp = 2*(q(3)*q(4)+q(1)*q(2));
a=zeros(3,1);
a(1,1) = asin(tmp);
if abs(tmp) <= 1-1e-5
    a(2,1) = -atan2(2*(q(2)*q(4)-q(1)*q(3)),q(1)^2-q(2)^2-q(3)^2+q(4)^2);
    a(3,1) = -atan2(2*(q(2)*q(3)-q(1)*q(4)),q(1)^2-q(2)^2+q(3)^2-q(4)^2);
else
    a(2,1) = atan2(2*(q(2)*q(4)+q(1)*q(3)),q(1)^2+q(2)^2-q(3)^2-q(4)^2);
    a(3,1) = 0;
end
end
%% 等效旋转矢量到四元数
function q=phi2q(phi)
phi_norm=norm(phi);
q=[cos(phi_norm/2);phi/phi_norm*sin(phi_norm/2)];
end
%% 四元数旋转三维状态
function r_out=q_rot(q,r)
q_i=q_inv(q);
r=[0;r];
r_out=q_multi(q_multi(q,r),q_i);
r_out = r_out(2:4,1);
end
%% 四元数乘法
function q_out=q_multi(q1,q2)
q1_0 = q1(1); q1_v = q1(2:4,1);
q2_0 = q2(1); q2_v = q2(2:4,1);
q_out = [q1_0*q2_0-q1_v'*q2_v;
         q1_0*q2_v+q2_0*q1_v+cross(q1_v,q2_v)];
end
%% 四元数取对偶
function q_out=q_inv(q)
q_out=[q(1);-q(2:4,1)];
end
%% 线性外推计算
function x_half=linear_push(x1,x2)
x_half = (3*x2-x1)/2; %x1为前2时刻，x2为前1时刻
end