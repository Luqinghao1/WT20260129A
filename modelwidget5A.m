function Dp = Composite_shale_oil_reservoir_fitfun(x,t)
% Ver:2025-12-1
% 2026-1-29:改fs1、fs2
% -------------------
x = [1e-2, 10, 1000, 20, 500, 0.4, 0.08,  1e-3, 1e-4, 20000];
t = logspace(-3,4,100);
clear global
% -------------------
global hDp hDp_der
% 拟合未知参数
kf = x(1);       % 内区的渗透率
M12 = x(2);      % 流度比 
L = x(3);        % 水平井长
Lf = x(4);       % 裂缝半长
rm = x(5);       % 改造区半径
omga1 = x(6);    % 内区储容比
omga2 = x(7);    % 外区储容比
remda1 = x(8);   % 内区窜流系数
remda2 = x(9);   % 外区窜流系数
re = x(10);      % 外区半径


eta12 = 0.2;      % 导压系数比
C     = 1e-7;  % 有因次井筒储存系数 (m3/MPa), 可作为输入或拟合参数
S     = 0.01;
 
% 可知参数
nf = 4;   % 裂缝条数
xwD = linspace(-0.9,0.9,nf);  % 裂缝位置
phi = 0.05;   % 基质孔隙度   
h = 20;     % 有效厚度
mu = 0.5;   % 粘度
B = 1.2;   % 体积系数
Ct = 5e-4;  % 综合压缩系数
q = 5;      % 定产生产(m3/s)
rw = 0.1;   % 井筒

% 无因次定义
CD = 0.159*C/ ( phi  * h * Ct * rw^2); % 无因次 CD
LfD = Lf/L;      % 无因次裂缝长度
rmD = rm/L;      % 无因次改造区半径
reD = re/L;      % 无因次外区半径

%关于时间tD的Laplace逆变换   %时间t以天为单位
tD =  14.4*kf*t/(phi*mu*Ct*L^2);  % t = phi_m*mu*Ct*L^2*tD/(14.4*kf);
PD=zeros(size(tD));

N=4;
% n = length(tD);
for i=1:length(tD)
     PD(i)=0; 
     y = 0;
     
     for j=1:N
        z=(log(2)/tD(i))*j;   
        v=0;
        for k=floor((j+1)/2):min(N/2,j)
            v=v+power(k,N/2+1)*factorial(2*k)/(factorial(N/2-k)*(factorial(k))^2*factorial(j-k)*factorial(2*k-j));
        end
        if(mod(N/2+j,2)==1), v=-v; else,   v=v; end

        fs1 = (omga1*(1-omga1)*z + remda1)/((1-omga1)*z + remda1);
        fs2 = eta12*(omga2*(1-omga2)*eta12*z + remda2)/((1-omga2)*eta12*z + remda2);
       
        y = PWD_inf(z,fs1,fs2,M12,LfD,rmD,nf,xwD,reD,S,CD);
        pf = y(nf+1);
        % ------------
        PD(i) = PD(i)+v*pf*log(2)/tD(i);
     end
     % -------------------
     % 摄动法考虑压敏
     gamaD = 0.02;
     if gamaD == 0
         PD(i) = PD(i);  else
         PD(i) = -1/gamaD*log(1-gamaD*PD(i));
     end
end

% t = phi_m*mu*Ct*L^2*tD/(14.4*kf);
Dp = 1.842e-3*q*mu*B*PD/(kf*h);    
Dp(isnan(Dp)) = 0;   
Dp_der = t(2:end).*diff(Dp)./diff(t);
if isempty(hDp)
    hDp = loglog(t,Dp,'-');  hold on
    hDp_der = loglog(t(2:end),Dp_der,'-');  
else
    set(hDp,'ydata',Dp);
    set(hDp_der,'ydata',Dp_der);
end
drawnow
% grid on
% xlabel('t,h');
% ylabel('\Deltap&d(\Deltapp)/d(ln(\Deltapp)),MPa');
% title('有因次压力响应')
% legend('压力','压力导数')

end 

% ===================================================================
function pf = PWD_inf(z,fs1,fs2,M12,LfD,rmD,nf,xwD,reD,S,CD)

% xwD = linspace(-0.9,0.9,nf);
ywD = zeros(size(xwD));

gama1 = sqrt(z*fs1);
gama2 = sqrt(z*fs2);
% --------------------------------------
% mAB=0;   %% 无穷大外边界
% mAB=besselk(1,gama2*reD)/besseli(1,gama2*reD);  %% 封闭边界
mAB=-besselk(0,gama2*reD)/besseli(0,gama2*reD); %% 定压边界
% --------------------------------------
Acup = M12*gama1*besselk(1,gama1*rmD)*(mAB*besseli(0,gama2*rmD)+besselk(0,gama2*rmD))+gama2*besselk(0,gama1*rmD)*(mAB*besseli(1,gama2*rmD)-besselk(1,gama2*rmD));
Acdown = M12*gama1*besseli(1,gama1*rmD)*(mAB*besseli(0,gama2*rmD)+besselk(0,gama2*rmD))-gama2*besseli(0,gama1*rmD)*(mAB*besseli(1,gama2*rmD)-besselk(1,gama2*rmD));
Ac = Acup/Acdown;

% 压力表达式
A = zeros(nf+1,nf+1);
for i=1:nf
    for j=1:nf
        y11=@(a)besselk(0,gama1*sqrt(((xwD(i)-xwD(j)-a).^2+(ywD(i)-ywD(j)).^2)))+Ac*besseli(0,gama1*sqrt(((xwD(i)-xwD(j)-a).^2+(ywD(i)-ywD(j)).^2)));
        pfD(i,j)=quadgk(y11,-LfD,LfD)/(M12*z*2*LfD);
        A(i,j)=z*pfD(i,j);
    end
end
A(:,nf+1) = -1;
A(nf+1,:) = z;
A(nf+1,nf+1) = 0;

b = zeros(nf+1,0);
b(nf+1) = 1;
pf = (A+eps(nf))\b';

% -----------------------
% 考虑井筒储存CD和表皮S
pf(nf+1)=(z*pf(nf+1)+S)/(z+CD*z^2*(z*pf(nf+1)+S));
% -------------------------
end  % function pf=PWD(z,fs1,fs2,M12)