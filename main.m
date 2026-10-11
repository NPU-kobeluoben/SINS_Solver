%% 单子样+前一周期
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
%% 统一格式 列状态
avp = avp';
imu = imu';
avp0=[avp0;0]; %补齐初始时间戳
%% 对齐数据
pos_align=[];
pos_align=align_init(pos_align,60);
%% 对初始姿态添加小的偏转扰动
bias = [1e-5;1e-6;1e-5;zeros(7,1)];
avp0_with_bias=avp0+bias;
%% 仿真缓存
Step = size(imu,2);
% avp_buff = zeros(10,Step+1); %俯仰、横滚、方位角，东-北-天速度，纬经高，时间t
% avp_buff(:,1) = [avp0;0];
t_align_stop = 0; %终止对准的后一步数
avp_cnt = 1;
%%
sins=[];
sins_without_align=[];
sins_without_align=sins_init(sins_without_align,avp0_with_bias,GM,Re,wie,e2,ff,S,C);

for t_setp = 1:Step
    t = imu(end,t_setp);
    sins_without_align=sins_update(sins_without_align,imu(:,t_setp));
    avp_without_align(:,t_setp)=sins_without_align.avp;
    
    if t<=pos_align.max_end_T
        pos_align.imu_sum = pos_align.imu_sum+imu(1:6,t_setp);
    else
        if ~pos_align.flag
            pos_align.flag = true;
            pos_align.t_end_T = imu(end,t_setp-1);
            pos_align=align_cal(pos_align);
            t_align_stop = t_setp;
            
            avp0_align = [c2a(pos_align.C);0;0;0;avp0(7:9,1);pos_align.t_end_T];
            sins=sins_init(sins,avp0_align,GM,Re,wie,e2,ff,S,C);
        end
        sins=sins_update(sins,imu(:,t_setp));
        avp_buff(:,avp_cnt)=sins.avp;
        avp_cnt = avp_cnt+1;
    end
end
%% 参考avp、无对齐avp提出
avp_ref=avp(:,t_align_stop:end);
avp_none_a=avp_without_align(:,t_align_stop:end);
%% imu 对齐时数据可视化 —— 确保静态对准
plot_imu_data(imu(1:6,1:t_align_stop-1),imu(end,1:t_align_stop-1));
%% 对比可视化
plot_here(avp_ref', avp_buff);
%% 是否对齐可视化对比
plot_err((avp_none_a(1:end-1,:)-avp_ref(1:end-1,:))', avp_buff(1:end-1,:)-avp_ref(1:end-1,:),avp_buff(end,:));

%% 对齐数据初始化
function pos_align=align_init(pos_align,end_T)
pos_align.imu_sum = zeros(6,1); %角增量和、线速度增量和
pos_align.max_end_T = end_T; %最大终止对齐时间
pos_align.t_end_T = end_T; %实际终止对齐时间
pos_align.imu_mean = zeros(6,1);
pos_align.C=zeros(3,3); %姿态阵 C^n_b
pos_align.flag=false;
end
%% 双矢量定姿粗对准
function pos_align=align_cal(pos_align)
pos_align.imu_mean=pos_align.imu_sum./pos_align.t_end_T;
w = pos_align.imu_mean(1:3,:);f = pos_align.imu_mean(4:6,:);
tmp1 = cross(f,w); tmp2 = cross(tmp1,f);
pos_align.C=[-tmp1'./norm(tmp1);
             tmp2'./norm(tmp2);
             f'./norm(f)];
end
%% SINS 初始化
function sins=sins_init(sins,avp0,GM,Re,wie,e2,ff,S,C)
sins.earth.GM=GM; sins.earth.Re=Re; sins.earth.wie=wie; sins.earth.ff=ff; sins.earth.C=C; sins.earth.S=S;
sins.earth.e2=e2;

sins.avp=[avp0(1:end-1,1);0];
sins.start_t = avp0(end,1); %结束对准时间戳
sins.imu_buff=zeros(7,2); %六轴+时间戳
sins.q = a2q(avp0(1:3,1));
sins.q_before=sins.q;
sins.init=false;
% 外推缓存
sins.wnie_before=zeros(3,1); sins.wnen_before=zeros(3,1);
sins.gn_before=zeros(3,1); sins.Mpv_before = zeros(3,3);
sins.v_before=zeros(3,1);
end
%% 单子样+前一时刻 SINS
function sins=sins_update(sins,imu)
if ~sins.init
    sins.imu_buff(:,1)=imu(:,1); %初始时刻前一字样保持相同——实为单子样
    dT=imu(end,1)-sins.start_t; %注意时间戳
else
    sins.imu_buff(:,1)=sins.imu_buff(:,2);
    dT=imu(end,1)-sins.imu_buff(end,1);
end
%数据提出
sins.imu_buff(:,2) = imu(:,1); %x\y\z陀螺角增量，x\y\z加计速度增量，t采样时间
%继承上一时刻姿态(四元数)
sins.q_before=sins.q;
%地球数据解算
sins.earth = eath_param_cal(sins.avp(7:9,1),sins.avp(4:6,1),sins.earth);
%圆锥/桨划、旋转补偿后的运动增量 —— 单子样+前一周期
[phi_b,dv_scul,dv_rot] = correct_rot_mo(sins.imu_buff);
%姿态更新(四元数)
sins.q=q_update(sins.q_before,phi_b,sins.earth.wnin,dT);
%姿态四元数到欧拉角
sins.avp(1:3,1)=q2a(sins.q);
% Mpv
sins.Mpv = [0 1/sins.earth.RMh 0; sec(sins.avp(7,1))/sins.earth.RNh 0 0;0 0 1];
%外推
if ~sins.init
    v_mid=sins.avp(4:6,1);
    wnie_mid=sins.earth.wnie; wnen_mid=sins.earth.wnen;
    gn_mid=sins.earth.gn;
    Mpv_mid = sins.Mpv;
    
    sins.init = true;
else
    v_mid=linear_push(sins.v_before,sins.avp(4:6,1));
    wnie_mid=linear_push(sins.wnie_before,sins.earth.wnie); wnen_mid=linear_push(sins.wnen_before,sins.earth.wnen);
    gn_mid=linear_push(sins.gn_before,sins.earth.gn);
    Mpv_mid = linear_push(sins.Mpv_before,sins.Mpv);
end
sins.wnie_before=sins.earth.wnie; sins.wnen_before=sins.earth.wnen;
sins.gn_before=sins.earth.gn;
sins.Mpv_before=sins.Mpv; sins.v_before=sins.avp(4:6,1);
%速度更新
sins.avp(4:6,1)=v_update(sins.avp(4:6,1),sins.imu_buff(4:6,2),dv_scul,dv_rot,sins.q_before,dT,v_mid,wnie_mid,wnen_mid,gn_mid);
%位置更新
sins.avp(7:9,1) = p_update(sins.avp(7:9,1),sins.v_before,sins.avp(4:6,1),dT,Mpv_mid);
%时间戳记录
sins.avp(10,1)=imu(end,1);
end

%% 位置更新
function p_m = p_update(p_before,v_before,v_m,dT,Mpv_mid)
p_m=p_before+Mpv_mid*(v_before+v_m)*dT/2;
end
%% 速度更新
function v_m=v_update(v_before,dv_imu,dv_scul,dv_rot,q_before,dT,v_mid,wnie_mid,wnen_mid,gn_mid)
wnin_mid=wnie_mid+wnen_mid;
I=eye(3);
dv_sf = (I-dT/2*M_cross(wnin_mid))*q_rot(q_before,dv_imu)+q_rot(q_before,dv_scul+dv_rot);
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
%% 圆锥、桨划补偿 —— 单字样+前一周期
function [phi_b,dv_scul,dv_rot] = correct_rot_mo(imu_buff)
dtheta = imu_buff(1:3,:); dv = imu_buff(4:6,:);
dtheta_m = imu_buff(1:3,2); dv_m = imu_buff(4:6,2);

dphi=cross(dtheta(:,1),dtheta(:,2))/12;
dv_scul=(cross(dtheta(:,1),dv(:,2))+cross(dv(:,1),dtheta(:,2)))/12;

dv_rot = cross(dtheta_m,dv_m)/2;
phi_b = dtheta_m+dphi;
end
%% 姿态阵到四元数
function q=c2q(C)
q=zeros(4,1);
if C(1,1)>=C(2,2)+C(3,3)
    q(2) = 0.5*sqrt(1+C(1,1)-C(2,2)+C(3,3));
    q(1) = (C(3,2)-C(2,3))/(4*q(2)); q(3)=(C(1,2)+C(2,1))/(4*q(2)); q(4)=(C(1,3)+C(3,1))/(4*q(2));
elseif C(2,2)>=C(1,1)+C(3,3)
    q(3) = 0.5*sqrt(1-C(1,1)+C(2,2)-C(3,3));
    q(1) = (C(1,3)-C(3,1))/(4*q(3)); q(2)=(C(1,2)+C(2,1))/(4*q(3)); q(4)=(C(2,3)+C(3,2))/(4*q(3));
elseif C(3,3)>=C(1,1)+C(2,2)
    q(4) = 0.5*sqrt(1-C(1,1)-C(2,2)+C(3,3));
    q(1) = (C(2,1)-C(1,2))/(4*q(4)); q(2)=(C(1,3)+C(3,1))/(4*q(4)); q(3)=(C(2,3)+C(3,2))/(4*q(4));
else
    q(1) = 0.5*sqrt(1+C(1,1)+C(2,2)+C(3,3));
    q(2) = (C(3,2)-C(2,3))/(4*q(1)); q(3)=(C(1,3)+C(3,1))/(4*q(4)); q(4)=(C(2,1)+C(1,2))/(4*q(4));
end
end
%% 姿态阵到欧拉角
function eu = c2a(C)
eu=zeros(3,1);
eu(1,1) = asin(C(3,2));
if abs(C(3,2))<=1-1e-5
    eu(2,1)=-atan2(C(3,1),C(3,3));
    eu(3,1)=-atan2(C(1,2),C(2,2));
else
    eu(2,1)=atan2(C(1,3),C(1,1));
    eu(3,1)=0;
end
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
%% 反对称矩阵生成
function m_out = M_cross(m)
m_out = [0 -m(3) m(2);
         m(3) 0  -m(1);
         -m(2) m(1) 0];
end
%% 线性外推计算
function x_half=linear_push(x1,x2)
x_half = (3*x2-x1)/2; %x1为前2时刻，x2为前1时刻
end