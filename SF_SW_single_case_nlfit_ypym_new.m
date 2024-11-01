tic

% dual cascade n=357238
% n=366500;
% n=357238;
% n=357238;
n=363000;

filename = ['PHI_p_m_0' num2str(n) '.txt'];
[t, yp, ym, wp, wm, a, b, c, d]=textread(filename,'%n%n%n%n%n%n%n%n%n');

r =t(1:600);
Nr=length(r);

% filter

% how many small scale points we want to keep, pick N1=11

% from the 11 to the end point, how many points we choose -- N2=47
% this number is calculated using log(t(11))-log(t(10)) as the step size in
% log t coordinates
N1 = 11;
N2 = 47;
N = N1+N2;

R = zeros(N,1);

R(1:N1) = r(1:N1);

rini = r(N1);
rend = r(end);
d_lnr = (log(rend)-log(rini))/N2;
R2 = exp([log(rini)+d_lnr/2:d_lnr:log(rend)-d_lnr/2])';
r_cut = exp([log(rini):d_lnr:log(rend)]);
nr_cut = zeros(N2+1,1);
nr_cut(1) = 11;
nr_cut(end) = length(r);
nr_now = 2;

for j=12:Nr-1
    if r(j)>r_cut(nr_now)
        nr_cut(nr_now)= j-1;
        nr_now=nr_now+1;
    end
end


R(N1+1:end) = R2;




% N =20; % how many data points after filtered
% rini = r(1);
% rend = r(end);
% d_lnr = (log(rend)-log(rini))/N;
% R = exp([log(rini)+d_lnr/2:d_lnr:log(rend)-d_lnr/2])';
% r_cut = exp([log(rini):d_lnr:log(rend)]);
% nr_cut = zeros(N+1,1);
% nr_cut(1) = 1;
% nr_cut(end) = length(r);
% nr_now = 2;
% 
% for j=2:Nr-1
%     if r(j)>r_cut(nr_now)
%         nr_cut(nr_now)= j-1;
%         nr_now=nr_now+1;
%     end
% end

v = -(yp)-ym;
V = zeros(N,1);
V(1:N1) = v(1:N1);
for j=1:N2
    V(N1+j) = mean( v( nr_cut(j):nr_cut(j+1) ) );
end

% for j=1:7
% %     R(j)=r(j);
%     V(j)=v(j);
% end

V0= V;
% ns=1;
% ne=Nr-Nr/8;
% R=r(ns:ne);
% V=yp(ns:ne)+ym(ns:ne);

NR=length(R);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  nonlinear fitting code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r = R/R(end); % normalize the scale
% the range of kf
k = flip(1./r);
% normalize the data
V_data = V./R; % note the division is R
Nfit = 10; % # of fitting scales

% now in 3D, e(2) = eps/kf^3

fun = @(e,r) 4/3*e(1) ...
    - 4*e(2)/e(3)^3*( sin(e(3)*r) - e(3)*r.*cos(e(3)*r) )./r.^3 ...
    - 4*e(4)/e(5)^3*( sin(e(5)*r) - e(5)*r.*cos(e(5)*r) )./r.^3 ...
    - 4*e(6)/e(7)^3*( sin(e(7)*r) - e(7)*r.*cos(e(7)*r) )./r.^3 ...
    - 4*e(8)/e(9)^3*( sin(e(9)*r) - e(9)*r.*cos(e(9)*r) )./r.^3 ...
    - 4*e(10)/e(11)^3*( sin(e(11)*r) - e(11)*r.*cos(e(11)*r) )./r.^3 ...
    - 4*e(12)/e(13)^3*( sin(e(13)*r) - e(13)*r.*cos(e(13)*r) )./r.^3 ...
    - 4*e(14)/e(15)^3*( sin(e(15)*r) - e(15)*r.*cos(e(15)*r) )./r.^3 ...
    - 4*e(16)/e(17)^3*( sin(e(17)*r) - e(17)*r.*cos(e(17)*r) )./r.^3 ...
    - 4*e(18)/e(19)^3*( sin(e(19)*r) - e(19)*r.*cos(e(19)*r) )./r.^3 ...
    - 4*e(20)/e(21)^3*( sin(e(21)*r) - e(21)*r.*cos(e(21)*r) )./r.^3 ;
%     - 4*e(22)/e(23)^3*( sin(e(23)*r) - e(23)*r.*cos(e(23)*r) )./r.^3 ;

% fun = @(e,r) 4/3*e(1) ...
%     - 4*e(2)*( sin(e(3)*r) - e(3)*r.*cos(e(3)*r) )./r.^3 ...
%     - 4*e(4)*( sin(e(5)*r) - e(5)*r.*cos(e(5)*r) )./r.^3 ...
%     - 4*e(6)*( sin(e(7)*r) - e(7)*r.*cos(e(7)*r) )./r.^3 ...
%     - 4*e(8)*( sin(e(9)*r) - e(9)*r.*cos(e(9)*r) )./r.^3 ...
%     - 4*e(10)*( sin(e(11)*r) - e(11)*r.*cos(e(11)*r) )./r.^3 ...
%     - 4*e(12)*( sin(e(13)*r) - e(13)*r.*cos(e(13)*r) )./r.^3 ;


% initial condition
x0 = ones(2*Nfit+1,1);

% set the bounds for kf and epsilon

lb =ones(1,2*Nfit+1);
ub =ones(1,2*Nfit+1);

% extend the scale range 
ext = 0.3;

lb(2) = -10*lb(2);
ub(2) =  10*ub(2);
lb(3) = min(k)*ext;
ub(3) = min(k)*(1+ext)/2;
% ub(3) = min(k);

lb(4) = -10*lb(4);
ub(4) =  10*ub(4);
lb(5) = min(k)*(1+ext)/2;
ub(5) = min(k);

lb(2*(Nfit-1)) = -10*lb(2*(Nfit-1));
ub(2*(Nfit-1)) =  10*ub(2*(Nfit-1));
lb(2*(Nfit-1)+1) = max(k);
ub(2*(Nfit-1)+1) = max(k)*(1+1/ext)/2;

lb(2*Nfit) = -10*lb(2*Nfit);
ub(2*Nfit) =  10*ub(2*Nfit);
lb(2*Nfit+1) = max(k)*(1+1/ext)/2;
% lb(2*Nfit+1) = max(k);
ub(2*Nfit+1) = max(k)/ext;



for j=3:Nfit-2
%     dk_ln = (log(max(k)) - log(min(k)))/(Nfit+1);
%     lb(2*j+1) = exp( log(min(k))+ (j-1)*dk_ln );
%     ub(2*j+1) = exp( log(min(k))+ (j)*dk_ln );
    
% tune the range of kf
    
    dk_ln = (log(max(k)) - log(min(k)))/(Nfit-4);
%     dk_ln = (log(max(k)) - log(min(k)*ext))/(Nfit+1);
    lb(2*j+1) = exp( log(min(k))+ (j-3)*dk_ln );
    ub(2*j+1) = exp( log(min(k))+ (j-2)*dk_ln );
    
    
    lb(2*j) = -10*lb(2*j);
    ub(2*j) = 10*ub(2*j);
end
lb(1) = -10;
ub(1) = 10;

% nonlinear fitting
e = lsqcurvefit(fun,x0,r,V_data,lb,ub);

% V_fit = 4/3*e(1)*r ...
%     - 4*e(2)/e(3)^3*( sin(e(3)*r) - e(3)*r.*cos(e(3)*r) )./r.^2 ...
%     - 4*e(4)/e(5)^3*( sin(e(5)*r) - e(5)*r.*cos(e(5)*r) )./r.^2 ...
%     - 4*e(6)/e(7)^3*( sin(e(7)*r) - e(7)*r.*cos(e(7)*r) )./r.^2 ...
%     - 4*e(8)/e(9)^3*( sin(e(9)*r) - e(9)*r.*cos(e(9)*r) )./r.^2 ...
%     - 4*e(10)/e(11)^3*( sin(e(11)*r) - e(11)*r.*cos(e(11)*r) )./r.^2 ...
%     - 4*e(12)/e(13)^3*( sin(e(13)*r) - e(13)*r.*cos(e(13)*r) )./r.^2 ...
%     - 4*e(14)/e(15)^3*( sin(e(15)*r) - e(15)*r.*cos(e(15)*r) )./r.^2 ...
%     - 4*e(16)/e(17)^3*( sin(e(17)*r) - e(17)*r.*cos(e(17)*r) )./r.^2 ...
%     - 4*e(18)/e(19)^3*( sin(e(19)*r) - e(19)*r.*cos(e(19)*r) )./r.^2 ...
%     - 4*e(20)/e(21)^3*( sin(e(21)*r) - e(21)*r.*cos(e(21)*r) )./r.^3 ...
%     - 4*e(22)/e(23)^3*( sin(e(23)*r) - e(23)*r.*cos(e(23)*r) )./r.^3 ;

% V_fit = 4/3*e(1)*r ...
%     - 4*e(2)*( sin(e(3)*r) - e(3)*r.*cos(e(3)*r) )./r.^2 ...
%     - 4*e(4)*( sin(e(5)*r) - e(5)*r.*cos(e(5)*r) )./r.^2 ...
%     - 4*e(6)*( sin(e(7)*r) - e(7)*r.*cos(e(7)*r) )./r.^2 ...
%     - 4*e(8)*( sin(e(9)*r) - e(9)*r.*cos(e(9)*r) )./r.^2 ...
%     - 4*e(10)*( sin(e(11)*r) - e(11)*r.*cos(e(11)*r) )./r.^2 ...
%     - 4*e(12)*( sin(e(13)*r) - e(13)*r.*cos(e(13)*r) )./r.^2 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% figure
% % loglog(t,yp,'b+','linewidth',1.3,'markersize',7)
% % hold on
% % loglog(t,-yp,'bo','linewidth',1.3,'markersize',7)
% loglog(r,V0./R,'g+','linewidth',1.3,'markersize',7)
% hold on
% loglog(r,-V0./R,'go','linewidth',1.3,'markersize',7)
% loglog(r,V_fit./r,'r','linewidth',1.3)
% loglog(r,-V_fit./r,'r--','linewidth',1.3)
% % loglog(r,(-2*ed*r+10e-7*r.^3),'k')
% % loglog(r,-(-2*ed*r+10e-7*r.^3),'k--')
% leg=legend('data +','data -','filtered data +','filtered data -','$$V_{fit}+$$','$$V_{fit}-$$');
% set(leg,'interpreter','latex')
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$Y^+$$','interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')

R_ext_ini = R(1)*ext;
R_ext_end = R(end)/ext;
% R_ext_end = R(end);
R_ext = exp([log(R_ext_ini):(log(R_ext_end)-log(R_ext_ini))/1e2:log(R_ext_end)])';

r_ext_ini = r(1)*ext;
r_ext_end = r(end)/ext;
% r_ext_end = r(end);
r_ext = exp([log(r_ext_ini):(log(r_ext_end)-log(r_ext_ini))/1e2:log(r_ext_end)])';

V_fit = 4/3*e(1)*r_ext ...
    - 4*e(2)/e(3)^3*( sin(e(3)*r_ext) - e(3)*r_ext.*cos(e(3)*r_ext) )./r_ext.^2 ...
    - 4*e(4)/e(5)^3*( sin(e(5)*r_ext) - e(5)*r_ext.*cos(e(5)*r_ext) )./r_ext.^2 ...
    - 4*e(6)/e(7)^3*( sin(e(7)*r_ext) - e(7)*r_ext.*cos(e(7)*r_ext) )./r_ext.^2 ...
    - 4*e(8)/e(9)^3*( sin(e(9)*r_ext) - e(9)*r_ext.*cos(e(9)*r_ext) )./r_ext.^2 ...
    - 4*e(10)/e(11)^3*( sin(e(11)*r_ext) - e(11)*r_ext.*cos(e(11)*r_ext) )./r_ext.^2 ...
    - 4*e(12)/e(13)^3*( sin(e(13)*r_ext) - e(13)*r_ext.*cos(e(13)*r_ext) )./r_ext.^2 ...
    - 4*e(14)/e(15)^3*( sin(e(15)*r_ext) - e(15)*r_ext.*cos(e(15)*r_ext) )./r_ext.^2 ...
    - 4*e(16)/e(17)^3*( sin(e(17)*r_ext) - e(17)*r_ext.*cos(e(17)*r_ext) )./r_ext.^2 ...
    - 4*e(18)/e(19)^3*( sin(e(19)*r_ext) - e(19)*r_ext.*cos(e(19)*r_ext) )./r_ext.^2 ...
    - 4*e(20)/e(21)^3*( sin(e(21)*r_ext) - e(21)*r_ext.*cos(e(21)*r_ext) )./r_ext.^2 ;
%     - 4*e(22)/e(23)^3*( sin(e(23)*r_ext) - e(23)*r_ext.*cos(e(23)*r_ext) )./r_ext.^2 ;
% figure
% loglog(t,v,'b+','linewidth',1.3,'markersize',7)
% hold on
% loglog(t,-v,'bo','linewidth',1.3,'markersize',7)
% loglog(R,V0,'g+','linewidth',1.3,'markersize',7)
% hold on
% loglog(R,-V0,'go','linewidth',1.3,'markersize',7)
% % loglog(R,V_fit*R(end),'r','linewidth',1.3)
% % loglog(R,-V_fit*R(end),'r--','linewidth',1.3)
% loglog(R_ext,V_fit*R(end),'r','linewidth',1.3)
% loglog(R_ext,-V_fit*R(end),'r--','linewidth',1.3)
% % loglog(r,(-2*ed*r+10e-7*r.^3),'k')
% % loglog(r,-(-2*ed*r+10e-7*r.^3),'k--')
% leg=legend('data +','data -','filtered data +','filtered data -','$$V_{fit}+$$','$$V_{fit}-$$');
% set(leg,'interpreter','latex')
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$Y^++Y^-$$','interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')


figure
% semilogx(750*t,v./t/7e5,'bO','linewidth',1.3,'markersize',7)
% hold on
semilogx(750*R,V0./R/7e5,'kO','linewidth',1.3,'markersize',7)
hold on
% semilogx(R,V_fit*R(end)./R,'r','linewidth',1.3)
semilogx(750*R_ext,V_fit./r_ext/7e5,'k','linewidth',1.3)
semilogx(750*R_ext,0*V_fit./r_ext/7e5,'--k','linewidth',1.3)
% loglog(r,(-2*ed*r+10e-7*r.^3),'k')
% loglog(r,-(-2*ed*r+10e-7*r.^3),'k--')
% leg=legend('data','filtered data','$$V_{fit}$$');
leg=legend('data','$$V_{fit}$$');
set(leg,'interpreter','latex')
xlabel('$$r\,\, (km)$$','interpreter','latex')
% ylabel('$$(Y^++Y^-)/r$$','interpreter','latex')
ylabel('$$(Y^++Y^-)/r\,\,(m^2s^{-3})$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')
xlim([min(750*t) max(750*t(600))])


% figure
% loglog(750*t,v,'b+','linewidth',1.3,'markersize',7)
% hold on
% loglog(750*t,-v,'bO','linewidth',1.3,'markersize',7)
% loglog(750*R,V0,'g+','linewidth',1.3,'markersize',7)
% loglog(750*R,-V0,'gO','linewidth',1.3,'markersize',7)
% loglog(750*R_ext,V_fit*R(end),'r-','linewidth',1.3)
% loglog(750*R_ext,-V_fit*R(end),'r--','linewidth',1.3)
% % loglog(r,(-2*ed*r+10e-7*r.^3),'k')
% % loglog(r,-(-2*ed*r+10e-7*r.^3),'k--')
% leg=legend('data+','data-','filtered data +','filtered data -','$$V_{fit}+$$','$$V_{fit}-$$');
% set(leg,'interpreter','latex')
% xlabel('$$r$$','interpreter','latex')
% % ylabel('$$(Y^++Y^-)/r$$','interpreter','latex')
% ylabel('$$Y^++Y^-$$','interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')
% 
% figure
% % semilogx(t,v./t,'bO','linewidth',1.3,'markersize',7)
% % hold on
% semilogx(R,V0./700/R,'bO','linewidth',1.3,'markersize',7)
% hold on
% % semilogx(R,V_fit*R(end)./R,'r','linewidth',1.3)
% semilogx(R_ext,V_fit./700/r_ext,'r','linewidth',1.3)
% % loglog(r,(-2*ed*r+10e-7*r.^3),'k')
% % loglog(r,-(-2*ed*r+10e-7*r.^3),'k--')
% % leg=legend('data','filtered data','$$V_{fit}$$');
% leg=legend('data','$$V_{fit}$$');
% set(leg,'interpreter','latex')
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$(Y^++Y^-)/r$$','interpreter','latex')
% % ylabel('$$(Y^+)/r$$','interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')

% figure
% % semilogx(k,e,'bO','LineWidth',1.3)
% % hold on
% for j=1:Nfit
%     semilogx(1./e(2*j+1)*R(end),e(2*j)*e(2*j+1)^3,'rd','LineWidth',1.3)
%     hold on
% end
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$\epsilon_i$$','interpreter','latex')
% % leg=legend('prescribed','fitted');
% % set(leg,'interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')

% figure
% % semilogx(k,e,'bO','LineWidth',1.3)
% % hold on
% % for j=1:Nfit
% %     semilogx(1./e(2*j+1)*R(end),e(2*j),'rd','LineWidth',1.3)
% %     hold on
% % end
% semilogx(1./e(3:2:end)*R(end),e(2:2:end),'rd','LineWidth',1.3)
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$\epsilon_i$$','interpreter','latex')
% % leg=legend('prescribed','fitted');
% % set(leg,'interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')


%%%%%%% construct the flux

F = zeros(1,2*Nfit+2);
F (1) = -e(1);
for j=1:Nfit
        F(2*j) = F(2*j-1);
        F(2*j+1) = F(2*j-1) + e(2*j);
end
F(end) = F(end-1);

kf = zeros(1,2*Nfit+2);
kf(1) = e(3)/2;
kf(end) = e(end)*2;

for j=2:2*Nfit+1
    if mod(j,2)==1
        kf(j) = e(j);
    end
end

for j=2:2*Nfit+1
    if mod(j,2)==0
        kf(j) = kf(j+1)-0.01*(kf(j+1)-kf(j-1));
    end
end
% 
% figure 
% semilogx(kf*pi/R(end),F,'b','LineWidth',1.3)
% hold on
% semilogx(kf*pi/R(end),0*F,'k--','LineWidth',1.3)
% xlabel('$$K$$','interpreter','latex')
% ylabel('$$F$$','interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')
% 
% 
% %%%%%%% construct the flux in R space
% 
% figure 
% semilogx(flip(1./kf(3:end-2)),flip(F(3:end-2)),'b','LineWidth',1.3)
% hold on
% semilogx(flip(1./kf(3:end-2)),0*flip(F(3:end-2)),'k--','LineWidth',1.3)
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$F$$','interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')
% 
% figure 
% semilogx(flip(R(end)./kf),flip(F),'b','LineWidth',1.3)
% hold on
% semilogx(flip(R(end)./kf),0*flip(F),'k--','LineWidth',1.3)
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$F$$','interpreter','latex')
% set(gca,'FontSize',14,'FontName','Times')
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F1=F;
for j=1:4
    F1(j)=F(5);
    F(end+1-j) = F(end-4);
end

figure 
semilogx(flip(750*R(end)./kf),1e9*flip(F1/7e5),'k','LineWidth',1.3)
hold on
semilogx(flip(750*R(end)./kf),0*flip(F1/7e5),'k--','LineWidth',1.3)
xlabel('$$r\,\, (km)$$','interpreter','latex')
ylabel('$$F\,\,(Jkg^{-1}s^{-1})$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')
xlim([min(750*t) max(750*t(1:600))])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FLUX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 
% 
% figure
% loglog(t,yp,'b-','LineWidth',1.3)
% hold on
% loglog(t,-yp,'b--','LineWidth',1.3)
% loglog(t,ym,'r-','LineWidth',1.3)
% loglog(t,-ym,'r-','LineWidth',1.3)
% loglog(t,t,'k','LineWidth',1.3)
% leg=legend('$$Y^+$$','$$-Y^+$$','$$Y^-$$','$$-Y^-$$');
% set(leg,'interpreter','latex')
% xlabel('$$t$$','interpreter','latex')
% set(gca,'fontsize',14,'fontname','times')
% 
% 
% figure
% loglog(t,yp+ym,'g-','LineWidth',1.3)
% hold on
% loglog(t,-(yp+ym),'g--','LineWidth',1.3)
% loglog(t,t,'k','LineWidth',1.3)
% leg=legend('$$Y^++Y^-$$','$$-Y^+-Y^-$$');
% set(leg,'interpreter','latex')
% xlabel('$$t$$','interpreter','latex')
% set(gca,'fontsize',14,'fontname','times')