Nf = 6;

load(['energy_sw_tot_nlfit_ytot_Nf_' num2str(Nf) '_t600.mat'])


N=40000;


% l=480*(1:1:1024);

Nl=256;
lmin = 480;
lmax = 480*600;
log_l_min = log(lmin);
log_l_max = log(lmax);
d_ln_l = (log_l_max - log_l_min)/(Nl-1);

l = exp([log_l_min:d_ln_l:log_l_max]);

eps_6 = zeros(N,Nl);


for n=1:N
    
    
%     e0 = 0;
    
    s_now=1;
    for p = Nl:-1:1
        
%         p_now = 0;
        
        if l(p) > 288000/k_inj_tot_ytot(n,s_now)
%             eps(n,p)=e0;
            p_now = p;
        else 
            eps_6(n,p_now) = e_inj_tot_ytot(n,s_now);
            eps_6(n,p_now-1) = e_inj_tot_ytot(n,s_now);
            s_now=s_now+1;
        end
        
        
        if s_now > Nf
%             eps(n,1:p) = zeros(1,p)*e0;
            break
        end
    end
end


% figure
% pcolor(tt(:,10:10:end),log(ll(:,10:10:end)),eps(10:10:end,:)')
% title('injection rate','interpreter','latex')
% xlabel('$$t$$','interpreter','latex')
% ylabel('$$ln(dt)$$','interpreter','latex')
% shading flat
% colorbar
% colormap bluewhitered
% 
% Flux_av = movmean(eps,128,1);
% 
% figure
% pcolor(tt(:,10:10:end),log(ll(:,10:10:end)),Flux_av(10:10:end,:)')
% title('Averaged injection rate over 128 datasets','interpreter','latex')
% xlabel('$$t$$','interpreter','latex')
% ylabel('$$ln(dt)$$','interpreter','latex')
% shading flat
% colorbar
% colormap bluewhitered


eps_av_6 = movmean(eps_6,1024,1);

% figure
% pcolor(tt(:,10:10:end),ll(:,10:10:end),eps_av(10:10:end,:)')
% title('Averaged injection rate over 1024 datasets','interpreter','latex')
% xlabel('$$t$$','interpreter','latex')
% ylabel('$$ln(dt)$$','interpreter','latex')
% set(gca,'yscal','log')
% shading flat
% colorbar
% colormap bluewhitered

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nf = 7;

load(['energy_sw_tot_nlfit_ytot_Nf_' num2str(Nf) '_t600.mat'])


N=40000;


change = zeros(N,Nf); 


isup = zeros(N,1);
isdown = zeros(N,1);


% l=480*(1:1:1024);

Nl=256;
lmin = 480;
lmax = 480*600;
log_l_min = log(lmin);
log_l_max = log(lmax);
d_ln_l = (log_l_max - log_l_min)/(Nl-1);

l = exp([log_l_min:d_ln_l:log_l_max]);

eps_7 = zeros(N,Nl);


for n=1:N
    
    
%     e0 = 0;
    
    s_now=1;
    for p = Nl:-1:1
        
%         p_now = 0;
        
        if l(p) > 288000/k_inj_tot_ytot(n,s_now)
%             eps(n,p)=e0;
            p_now = p;
        else 
            eps_7(n,p_now) = e_inj_tot_ytot(n,s_now);
            eps_7(n,p_now-1) = e_inj_tot_ytot(n,s_now);
            s_now=s_now+1;
        end
        
        
        if s_now > Nf
%             eps(n,1:p) = zeros(1,p)*e0;
            break
        end
    end
end



[tt, ll] = meshgrid(time_tot_ytot,700*l');


% figure
% pcolor(tt(:,10:10:end),log(ll(:,10:10:end)),eps(10:10:end,:)')
% title('injection rate','interpreter','latex')
% xlabel('$$t$$','interpreter','latex')
% ylabel('$$ln(dt)$$','interpreter','latex')
% shading flat
% colorbar
% colormap bluewhitered
% 
% Flux_av = movmean(eps,128,1);
% 
% figure
% pcolor(tt(:,10:10:end),log(ll(:,10:10:end)),Flux_av(10:10:end,:)')
% title('Averaged injection rate over 128 datasets','interpreter','latex')
% xlabel('$$t$$','interpreter','latex')
% ylabel('$$ln(dt)$$','interpreter','latex')
% shading flat
% colorbar
% colormap bluewhitered


eps_av_7 = movmean(eps_7,1024,1);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eps_av_0 = (eps_6+eps_7)/2/7e5;
% eps_av_0 = eps_6/7e5;

max_eps_av_0 = max(max(abs(eps_av_0)));

eps_av = (eps_av_6+eps_av_7)/2/7e5;

max_eps_av = max(max(abs(eps_av)));

eps_av = eps_av*max_eps_av_0/max_eps_av;

tt = (tt-336570)*0.05/9;

figure
pcolor(tt(:,10:10:end),ll(:,10:10:end),eps_av(10:10:end,:)')
title('Energy injection rate','interpreter','latex')
xlabel('Day of Year 1996','interpreter','latex')
ylabel('$$r$$','interpreter','latex')
set(gca,'yscal','log')
shading flat
colorbar
colormap bluewhitered
set(gca,'FontSize',22,'FontName','Times')

eps = (eps_6+eps_7)/2/7e5;
eps_mean = movmean(eps,floor(256/20),2);

max_eps_mean = max(max(abs(eps_mean)));

eps_mean=eps_mean*max_eps_av_0/max_eps_mean;


figure
pcolor(tt(:,10:10:end),ll(:,10:10:end),1e9*eps_mean(10:10:end,:)')
title('Energy injection rate','interpreter','latex')
xlabel('Day of Year 1996','interpreter','latex')
ylabel('$$r$$','interpreter','latex')
set(gca,'yscal','log')
shading flat
colorbar
colormap bluewhitered
set(gca,'FontSize',22,'FontName','Times')
% 
% figure
% pcolor(tt(:,10:10:end),ll(:,10:10:end),-1e9*eps_mean(10:10:end,:)')
% title('-Energy injection rate','interpreter','latex')
% xlabel('Day of Year 1996','interpreter','latex')
% ylabel('$$r$$','interpreter','latex')
% set(gca,'yscal','log')
% shading flat
% colorbar
% colormap bluewhitered
% set(gca,'FontSize',22,'FontName','Times')


e_av_mean = movmean(eps_av,floor(256/20),2);



figure
pcolor(tt(:,10:10:end),ll(:,10:10:end),e_av_mean(10:10:end,:)')
title('Averaged injection rate over 1024 datasets','interpreter','latex')
xlabel('$$t$$','interpreter','latex')
ylabel('$$r\,(km)$$','interpreter','latex')
set(gca,'yscal','log')
shading flat
colorbar
colormap bluewhitered
set(gca,'FontSize',14,'FontName','Times')
ylim([lmin*750 lmax*750])
%     
%     if F0<0
%         isup(n)=1;
%     end
%     
%     for s = 1:Nf
%         
%         F= F0 + e_inj_tot_ytot(n,s);
%         
%         if F*F0<0
%             if F0<0
%                 change(n,s)=1;
%             else
%                 change(n,s)=-1;
%             end
%             
%         end
%         
%         F0=F;
%     end
%     
%     if F0>0
%         isdown(n)=1;
%     end
% end


NL=75;
Lmin = 480;
Lmax = 480*600;
log_L_min = log(Lmin);
log_L_max = log(Lmax);
d_ln_L = (log_L_max - log_L_min)/(NL-1);

L = exp([log_L_min:d_ln_L:log_L_max]);

eps=1e9*eps_av_0;

max_eps = max(max(eps));

Ne=75;
e = (1:1:Ne)/Ne*max_eps;

pdf_eps=zeros(NL,Ne);

for j=1:Nl
    for p=1:N
        if eps(p,j)>1/3e2*max_eps
        q=ceil(eps(p,j)/max_eps*Ne);
        s=ceil((log(l(j))-log_L_min)/(log_L_max - log_L_min)*NL);
        pdf_eps(s,q)=pdf_eps(s,q)+1;
        end
    end
end

[ee,LL]=meshgrid(e,L*750);

figure
pcolor(ee,LL,pdf_eps/sum(sum(pdf_eps)))
title('injection pdf','interpreter','latex')
xlabel('$$\epsilon$$','interpreter','latex')
ylabel('$$r$$','interpreter','latex')
set(gca,'yscal','log')
shading flat
colorbar
% colormap bluewhitered
set(gca,'FontSize',14,'FontName','Times')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% log epsilon axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_eps = max(max(eps));
min_eps = 1/3e2*max_eps;

pos_eps=zeros(N,Nl);
for j=1:Nl
    for p=1:N
        if eps(p,j)>min_eps
        pos_eps(p,j)=eps(p,j);
        end
    end
end

ln_pos_eps = log(pos_eps);

% Ne=100;
dln_e = (log(max_eps)-log(min_eps))/(Ne-1);
ln_e = log(min_eps):dln_e:log(max_eps);
e=exp(ln_e);

pdf_eps=zeros(NL,Ne);

for j=1:Nl
    for p=1:N
        if eps(p,j)>1/3e2*max_eps
        q=ceil(ln_pos_eps(p,j)/log(max_eps)*Ne);
        s=ceil((log(l(j))-log_L_min)/(log_L_max - log_L_min)*NL);
        pdf_eps(s,q)=pdf_eps(s,q)+1;
        end
    end
end

[ee,LL]=meshgrid(e,L*750);

figure
subplot(1,2,1)
pcolor(ee,LL,pdf_eps/sum(sum(pdf_eps)))
title('injection pdf','interpreter','latex')
xlabel('$$\epsilon(Jkg^{-1}s^{-1})$$','interpreter','latex')
ylabel('$$r(km)$$','interpreter','latex')
set(gca,'yscal','log')
set(gca,'xscal','log')
shading flat
colorbar
xlim([1e2 max_eps])
caxis([0 3.5e-3])
% colormap bluewhitered
set(gca,'FontSize',14,'FontName','Times')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% log epsilon axis; disspation pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

max_eps = max(max(-eps));
min_eps = 1/3e2*max_eps;

neg_eps=zeros(N,NL);
for j=1:Nl
    for p=1:N
        if -eps(p,j)>min_eps
        neg_eps(p,j)=-eps(p,j);
        end
    end
end

ln_pos_eps = log(neg_eps);

% Ne=200;
dln_e = (log(max_eps)-log(min_eps))/(Ne-1);
ln_e = log(min_eps):dln_e:log(max_eps);
e=exp(ln_e);

pdf_eps=zeros(NL,Ne);

for j=1:Nl
    for p=1:N
        if -eps(p,j)>1/3e2*max_eps
        q=ceil(ln_pos_eps(p,j)/log(max_eps)*Ne);
        s=ceil((log(l(j))-log_L_min)/(log_L_max - log_L_min)*NL);
        pdf_eps(s,q)=pdf_eps(s,q)+1;
        end
    end
end

[ee,LL]=meshgrid(e,L*750);

subplot(1,2,2)
pcolor(ee,LL,pdf_eps/sum(sum(pdf_eps)))
title('dissipation pdf','interpreter','latex')
xlabel('$$-\epsilon(Jkg^{-1}s^{-1})$$','interpreter','latex')
ylabel('$$r(km)$$','interpreter','latex')
set(gca,'yscal','log')
set(gca,'xscal','log')
shading flat
colorbar
xlim([1e2 max_eps])
caxis([0 3.5e-3])
% colormap bluewhitered
set(gca,'FontSize',14,'FontName','Times')
