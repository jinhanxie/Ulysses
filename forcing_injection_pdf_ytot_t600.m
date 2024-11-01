
Nf = 6;

load(['energy_sw_tot_nlfit_ytot_Nf_' num2str(Nf) '_t600.mat'])

% l_min = 1/975.7850*480*1024;
% l_max = 1*480*1024;




Nl=128;
lmin = 480;
lmax = 480*600;
log_l_min = log(lmin);
log_l_max = log(lmax);
d_ln_l = (log_l_max - log_l_min)/(Nl-1);

l = exp([log_l_min:d_ln_l:log_l_max]);


Nedge = 20;
d_ln_edge = (log(lmax)-log(lmin))/Nedge;
edges = exp(log(lmin):d_ln_edge:log(lmax));


N = 40000;
eps = zeros(N,Nl);


for n=1:N
    
    
%     e0 = 0;
    
    s_now=1;
    for p = Nl:-1:1
        
%         p_now = 0;
        
        if l(p) > 480*600/k_inj_tot_ytot(n,s_now)
%             eps(n,p)=e0;
            p_now = p;
        else 
            eps(n,p_now) = e_inj_tot_ytot(n,s_now);
%             eps(n,p_now-1) = e_inj_tot_ytot(n,s_now);
            s_now=s_now+1;
        end
        
        
        if s_now > Nf
%             eps(n,1:p) = zeros(1,p)*e0;
            break
        end
    end
end


H_e_pos = zeros(N,Nl);
H_e_neg = zeros(N,Nl);

Nerror = 4000;

for p = 1:N
    e_max = max(eps(p,:));
    e_min = min(eps(p,:));
    for q = 1:Nl
        if eps(p,q)>0 && eps(p,q)>e_max/Nerror
            H_e_pos(p,q) = 1;
        end
        if eps(p,q)<0 && eps(p,q)<e_min/Nerror
            H_e_neg(p,q) = 1;
        end
    end
end

Neav = floor(Nl/(Nedge-1));

e_sum_1 = sum(eps,1);
emean_1 = movmean(e_sum_1,Neav);

e_sum_pos_1 = sum(eps.*H_e_pos,1);
emean_pos_1 = movmean(e_sum_pos_1,Neav);
e_sum_neg_1 = sum(eps.*H_e_neg,1);
emean_neg_1 = movmean(e_sum_neg_1,Neav);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Nf = 7;

load(['energy_sw_tot_nlfit_ytot_Nf_' num2str(Nf) '_t600.mat'])

% l_min = 1/975.7850*480*1024;
% l_max = 1*480*1024;
% 
% 
% Nedge = 20;
% d_ln_edge = (log(l_max)-log(l_min))/Nedge;
% edges = exp(log(l_min):d_ln_edge:log(l_max));

% Nl=256;
% lmin = 480;
% lmax = 480*1024;
% log_l_min = log(lmin);
% log_l_max = log(lmax);
% d_ln_l = (log_l_max - log_l_min)/(Nl-1);
% 
% l = exp([log_l_min:d_ln_l:log_l_max]);
% 
% 
% N = 40000;
eps = zeros(N,Nl);


for n=1:N
    
    
%     e0 = 0;
    
    s_now=1;
    for p = Nl:-1:1
        
%         p_now = 0;
        
        if l(p) > 480*600/k_inj_tot_ytot(n,s_now)
%             eps(n,p)=e0;
            p_now = p;
        else 
            eps(n,p_now) = e_inj_tot_ytot(n,s_now);
%             eps(n,p_now-1) = e_inj_tot_ytot(n,s_now);
            s_now=s_now+1;
        end
        
        
        if s_now > Nf
%             eps(n,1:p) = zeros(1,p)*e0;
            break
        end
    end
end


H_e_pos = zeros(N,Nl);
H_e_neg = zeros(N,Nl);

Nerror = 4000;

for p = 1:N
    e_max = max(eps(p,:));
    e_min = min(eps(p,:));
    for q = 1:Nl
        if eps(p,q)>0 && eps(p,q)>e_max/Nerror
            H_e_pos(p,q) = 1;
        end
        if eps(p,q)<0 && eps(p,q)<e_min/Nerror
            H_e_neg(p,q) = 1;
        end
    end
end

Neav = floor(Nl/(Nedge-1));

e_sum_2 = sum(eps,1);
emean_2 = movmean(e_sum_2,Neav);

e_sum_pos_2 = sum(eps.*H_e_pos,1);
emean_pos_2 = movmean(e_sum_pos_2,Neav);
e_sum_neg_2 = sum(eps.*H_e_neg,1);
emean_neg_2 = movmean(e_sum_neg_2,Neav);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

emean = (emean_1+emean_2)/2;
emean_pos = (emean_pos_1+emean_pos_2)/2;
emean_neg = (emean_neg_1+emean_neg_2)/2;

emean=movmean(emean,Neav);
emean_pos=movmean(emean_pos,Neav);
emean_neg=movmean(emean_neg,Neav);

Hpos = zeros(Nf,40000);
Hneg = zeros(Nf,40000);

% when e<e/Nsmall will be considered as noise

Nsmall = 4;

for j=1:40000
%     e_max = max(e_inj_tot_ytot(j,:));
%     e_min = min(e_inj_tot_ytot(j,:));
    e_abs = max(e_inj_tot_ytot(j,:)) - min(e_inj_tot_ytot(j,:));
    for p = 1:Nf
%         if e_inj_tot_ytot(j,p)>e_max/Nsmall
        if e_inj_tot_ytot(j,p)>e_abs/Nsmall
            Hpos(p,j)=1;
        end
%         if e_inj_tot_ytot(j,p)<e_min/Nsmall
        if e_inj_tot_ytot(j,p)<e_abs/Nsmall
            Hneg(p,j)=1;
        end
    end

end


figure

h = histogram(Hpos.*1./k_inj_tot_ytot'*480*1024,edges);
dBin=h.BinEdges(2)-h.BinEdges(1);
h_L=h.BinEdges(2:end)-dBin/2;
l_pos_pdf = h.Values/sum(sum(Hpos));
semilogx(h_L,h.Values/sum(sum(Hpos)),'linewidth',1.3)
hold on
semilogx(h_L,movmean(l_pos_pdf,3),'linewidth',1.3)
title('pdf of forcing scale (positive)','interpreter','latex')
xlabel('$$r$$','interpreter','latex')
ylabel('$$pdf$$','interpreter','latex')

% figure
% h = histogram(Hpos.*1./k_inj_tot_ytot'*480*1024,edges,'Normalization','probability');
% title('pdf of forcing scale (positive)','interpreter','latex')
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$pdf$$','interpreter','latex')
% set(gca,'xscal','log')
% 
% figure
% 
% h = histogram(Hneg.*1./k_inj_tot_ytot'*480*1024,edges);
% dBin=h.BinEdges(2)-h.BinEdges(1);
% h_L=h.BinEdges(2:end)-dBin/2;
% l_neg_pdf = h.Values/sum(sum(Hneg));
% semilogx(h_L,h.Values/sum(sum(Hneg)),'linewidth',1.3)
% title('pdf of forcing scale (negative)','interpreter','latex')
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$pdf$$','interpreter','latex')
% 
% figure
% h = histogram(Hneg.*1./k_inj_tot_ytot'*480*1024,edges,'Normalization','probability');
% title('pdf of forcing scale (negative)','interpreter','latex')
% xlabel('$$r$$','interpreter','latex')
% ylabel('$$pdf$$','interpreter','latex')
% set(gca,'xscal','log')

% energy injection rate




figure
% semilogx(h_L,movmean(l_pos_pdf,1),'LineWidth',2)
% hold on
% semilogx(h_L,movmean(l_neg_pdf,1),'LineWidth',2)
% semilogx(l,emean/sum(emean)*length(l)/length(h_L),'LineWidth',2)
% hold on
% semilogx(l,emean_pos/sum(emean_pos)*length(l)/length(h_L),'LineWidth',2)
% semilogx(l,emean_neg/sum(emean_neg)*length(l)/length(h_L),'LineWidth',2)
% semilogx(l*750,emean*length(l)/length(h_L)/N/7e5,'LineWidth',2)
% hold on
semilogx(l*750,1e9*emean_pos*length(l)/length(h_L)/N/7e5,'r','LineWidth',2)
hold on
semilogx(l*750,1e9*emean_neg*length(l)/length(h_L)/N/7e5,'b','LineWidth',2)
semilogx(l*750,0*emean_pos*length(l)/length(h_L)/N/7e5,'k--','LineWidth',1)
xlim([lmin*750 lmax*750])
% leg=legend('total injection distribution','injection distribution ($$\epsilon>0$$)','dissipation distribution ($$\epsilon<0$$)');
leg=legend('injection distribution','dissipation distribution');
set(leg,'interpreter','latex')
xlabel('$$r\,(km)$$','interpreter','latex')
ylabel('$$\epsilon\,(m^2s^{-3})$$','interpreter','latex')
set(gca,'fontname','times','fontsize',14)


figure
semilogy(1e9*emean_pos*length(l)/length(h_L)/N/7e5,l*750,'r','LineWidth',2)
hold on
semilogy(1e9*emean_neg*length(l)/length(h_L)/N/7e5,l*750,'b','LineWidth',2)
semilogy(0*emean_pos*length(l)/length(h_L)/N/7e5,l*750,'k--','LineWidth',1)
ylim([lmin*750 lmax*750])
% leg=legend('total injection distribution','injection distribution ($$\epsilon>0$$)','dissipation distribution ($$\epsilon<0$$)');
leg=legend('injection distribution','dissipation distribution');
set(leg,'interpreter','latex')
% ylabel('$$r\,(km)$$','interpreter','latex')
xlabel('$$\epsilon\,(m^2s^{-3})$$','interpreter','latex')
set(gca,'fontname','times','fontsize',14)
