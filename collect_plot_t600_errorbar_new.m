load('energy_sw_tot_nlfit_ytot_Nf_6_t600.mat')
tt = time_tot_ytot;
einj6 = sum(max(e_inj_tot_ytot,0),2);
edis6 = sum(min(e_inj_tot_ytot,0),2);
Fup6 = e_up_tot_ytot;
Fdown6 = einj6+edis6-Fup6;

load('energy_sw_tot_nlfit_ytot_Nf_7_t600.mat')

einj7 = sum(max(e_inj_tot_ytot,0),2);
edis7 = sum(min(e_inj_tot_ytot,0),2);
Fup7 = e_up_tot_ytot;
Fdown7 = einj7+edis7-Fup7;

einj = (einj6+einj7)/2/7e5;
edis = (edis6+edis7)/2/7e5;
Fup = (Fup6+Fup7)/2/7e5;
Fdown = (Fdown6+Fdown7)/2/7e5;

% 
% einj = [einj6; einj7]/7e5;
% edis = [edis6; edis7]/7e5;
% Fup = [Fup6; Fup7]/7e5;
% Fdown = [Fdown6; Fdown7]/7e5;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% four cases
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = length(einj); 
einj_dual = NaN(1,N);
edis_dual = NaN(1,N);
Fup_dual = NaN(1,N);
Fdown_dual = NaN(1,N);
H_dual = zeros(1,N);
for j=1:N
    if Fup(j)>0 && Fdown(j)>0
        einj_dual(j) = einj(j);
        edis_dual(j) = edis(j);
        Fup_dual(j) = Fup(j);
        Fdown_dual(j) = Fdown(j);
        H_dual(j)=1;
    end
end

einj_dual = einj_dual(~isnan(einj_dual));
edis_dual = edis_dual(~isnan(edis_dual));
Fup_dual = Fup_dual(~isnan(Fup_dual));
Fdown_dual = Fdown_dual(~isnan(Fdown_dual));


N = length(einj); 
einj_ndual = NaN(1,N);
edis_ndual = NaN(1,N);
Fup_ndual = NaN(1,N);
Fdown_ndual = NaN(1,N);

for j=1:N
    if Fup(j)<0 && Fdown(j)<0
        einj_ndual(j) = einj(j);
        edis_ndual(j) = edis(j);
        Fup_ndual(j) = Fup(j);
        Fdown_ndual(j) = Fdown(j);
    end
end

einj_ndual = einj_ndual(~isnan(einj_ndual));
edis_ndual = edis_ndual(~isnan(edis_ndual));
Fup_ndual = Fup_ndual(~isnan(Fup_ndual));
Fdown_ndual = Fdown_ndual(~isnan(Fdown_ndual));


N = length(einj); 
einj_forw = NaN(1,N);
edis_forw = NaN(1,N);
Fup_forw = NaN(1,N);
Fdown_forw = NaN(1,N);

for j=1:N
    if Fup(j)<0 && Fdown(j)>0
        einj_forw(j) = einj(j);
        edis_forw(j) = edis(j);
        Fup_forw(j) = Fup(j);
        Fdown_forw(j) = Fdown(j);
    end
end

einj_forw = einj_forw(~isnan(einj_forw));
edis_forw = edis_forw(~isnan(edis_forw));
Fup_forw = Fup_forw(~isnan(Fup_forw));
Fdown_forw = Fdown_forw(~isnan(Fdown_forw));

N = length(einj); 
einj_inv = NaN(1,N);
edis_inv = NaN(1,N);
Fup_inv = NaN(1,N);
Fdown_inv = NaN(1,N);

for j=1:N
    if Fup(j)>0 && Fdown(j)<0
        einj_inv(j) = einj(j);
        edis_inv(j) = edis(j);
        Fup_inv(j) = Fup(j);
        Fdown_inv(j) = Fdown(j);
    end
end

einj_inv = einj_inv(~isnan(einj_inv));
edis_inv = edis_inv(~isnan(edis_inv));
Fup_inv = Fup_inv(~isnan(Fup_inv));
Fdown_inv = Fdown_inv(~isnan(Fdown_inv));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% line with error bar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n_count = 1000;



einj_edis_dual = [einj_dual; edis_dual];
einj_edis_dual_sort = sortrows(einj_edis_dual')';
N = length(einj_dual);
n = floor(N/n_count);
einj_dual_n=zeros(1,n); 
edis_dual_n=zeros(1,n); 
var_einj_edis_dual_n=zeros(1,n); 
for j=1:n-1
    einj_dual_n(j) = mean(einj_edis_dual_sort(1,(j-1)*n_count+1:j*n_count));
    edis_dual_n(j)= mean(einj_edis_dual_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_edis_dual_n(j) = mean((einj_edis_dual_sort(2,(j-1)*n_count+1:j*n_count)-edis_dual_n(j)).^2)^0.5;
end
einj_dual_n(n) = mean(einj_edis_dual_sort(1,(n-1)*n_count+1:end));
edis_dual_n(n)= mean(einj_edis_dual_sort(2,(n-1)*n_count+1:end));
var_einj_edis_dual_n(n) = mean((einj_edis_dual_sort(2,(n-1)*n_count+1:end)-edis_dual_n(n)).^2)^0.5;


einj_edis_forw = [einj_forw; edis_forw];
einj_edis_forw_sort = sortrows(einj_edis_forw')';
N = length(einj_forw);
n = floor(N/n_count);
einj_forw_n=zeros(1,n); 
edis_forw_n=zeros(1,n); 
var_einj_edis_forw_n=zeros(1,n); 
for j=1:n-1
    einj_forw_n(j) = mean(einj_edis_forw_sort(1,(j-1)*n_count+1:j*n_count));
    edis_forw_n(j)= mean(einj_edis_forw_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_edis_forw_n(j) = mean((einj_edis_forw_sort(2,(j-1)*n_count+1:j*n_count)-edis_forw_n(j)).^2)^0.5;
end
einj_forw_n(n) = mean(einj_edis_forw_sort(1,(n-1)*n_count+1:end));
edis_forw_n(n)= mean(einj_edis_forw_sort(2,(n-1)*n_count+1:end));
var_einj_edis_forw_n(n) = mean((einj_edis_forw_sort(2,(n-1)*n_count+1:end)-edis_forw_n(n)).^2)^0.5;

einj_edis_inv = [einj_inv; edis_inv];
einj_edis_inv_sort = sortrows(einj_edis_inv')';
N = length(einj_inv);
n = floor(N/n_count);
einj_inv_n=zeros(1,n); 
edis_inv_n=zeros(1,n); 
var_einj_edis_inv_n=zeros(1,n); 
for j=1:n-1
    einj_inv_n(j) = mean(einj_edis_inv_sort(1,(j-1)*n_count+1:j*n_count));
    edis_inv_n(j)= mean(einj_edis_inv_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_edis_inv_n(j) = mean((einj_edis_inv_sort(2,(j-1)*n_count+1:j*n_count)-edis_inv_n(j)).^2)^0.5;
end
einj_inv_n(n) = mean(einj_edis_inv_sort(1,(n-1)*n_count+1:end));
edis_inv_n(n)= mean(einj_edis_inv_sort(2,(n-1)*n_count+1:end));
var_einj_edis_inv_n(n) = mean((einj_edis_inv_sort(2,(n-1)*n_count+1:end)-edis_inv_n(n)).^2)^0.5;

einj_edis_ndual = [einj_ndual; edis_ndual];
einj_edis_ndual_sort = sortrows(einj_edis_ndual')';
N = length(einj_ndual);
n = floor(N/n_count);
einj_ndual_n=zeros(1,n); 
edis_ndual_n=zeros(1,n); 
var_einj_edis_ndual_n=zeros(1,n); 
for j=1:n-1
    einj_ndual_n(j) = mean(einj_edis_ndual_sort(1,(j-1)*n_count+1:j*n_count));
    edis_ndual_n(j)= mean(einj_edis_ndual_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_edis_ndual_n(j) = mean((einj_edis_ndual_sort(2,(j-1)*n_count+1:j*n_count)-edis_ndual_n(j)).^2)^0.5;
end
einj_ndual_n(n) = mean(einj_edis_ndual_sort(1,(n-1)*n_count+1:end));
edis_ndual_n(n)= mean(einj_edis_ndual_sort(2,(n-1)*n_count+1:end));
var_einj_edis_ndual_n(n) = mean((einj_edis_ndual_sort(2,(n-1)*n_count+1:end)-edis_ndual_n(n)).^2)^0.5;

figure
plot(einj_dual_n,-edis_dual_n,'o','linewidth',1.3)
hold on
plot(einj_forw_n,-edis_forw_n,'o','linewidth',1.3)
plot(einj_inv_n,-edis_inv_n,'o','linewidth',1.3)
plot(einj_ndual_n,-edis_ndual_n,'o','linewidth',1.3)
patch([einj_dual_n flip(einj_dual_n)], [-edis_dual_n-var_einj_edis_dual_n flip(-edis_dual_n+var_einj_edis_dual_n)],[0 0.4470 0.7410], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([einj_forw_n flip(einj_forw_n)], [-edis_forw_n-var_einj_edis_forw_n flip(-edis_forw_n+var_einj_edis_forw_n)],[0.8500 0.3250 0.0980], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([einj_inv_n flip(einj_inv_n)], [-edis_inv_n-var_einj_edis_inv_n flip(-edis_inv_n+var_einj_edis_inv_n)],[0.9290 0.6940 0.1250], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([einj_ndual_n flip(einj_ndual_n)], [-edis_ndual_n-var_einj_edis_ndual_n flip(-edis_ndual_n+var_einj_edis_ndual_n)],[0.4940 0.1840 0.5560], 'FaceAlpha',0.2, 'EdgeColor','none')
leg=legend('dual (divergent)','forward','inverse','dual (convergent)');
set(leg,'interpreter','latex')
xlabel('$$\epsilon_{inj}$$','interpreter','latex')
ylabel('$$\epsilon_{diss}$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

einj_Fup_dual = [einj_dual; Fup_dual];
einj_Fup_dual_sort = sortrows(einj_Fup_dual')';
N = length(einj_dual);
n = floor(N/n_count);
einj_dual_n=zeros(1,n); 
Fup_dual_n=zeros(1,n); 
var_einj_Fup_dual_n=zeros(1,n); 
for j=1:n-1
    einj_dual_n(j) = mean(einj_Fup_dual_sort(1,(j-1)*n_count+1:j*n_count));
    Fup_dual_n(j)= mean(einj_Fup_dual_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_Fup_dual_n(j) = mean((einj_Fup_dual_sort(2,(j-1)*n_count+1:j*n_count)-Fup_dual_n(j)).^2)^0.5;
end
einj_dual_n(n) = mean(einj_Fup_dual_sort(1,(n-1)*n_count+1:end));
Fup_dual_n(n)= mean(einj_Fup_dual_sort(2,(n-1)*n_count+1:end));
var_einj_Fup_dual_n(n) = mean((einj_Fup_dual_sort(2,(n-1)*n_count+1:end)-Fup_dual_n(n)).^2)^0.5;


einj_Fup_forw = [einj_forw; Fup_forw];
einj_Fup_forw_sort = sortrows(einj_Fup_forw')';
N = length(einj_forw);
n = floor(N/n_count);
einj_forw_n=zeros(1,n); 
Fup_forw_n=zeros(1,n); 
var_einj_Fup_forw_n=zeros(1,n); 
for j=1:n-1
    einj_forw_n(j) = mean(einj_Fup_forw_sort(1,(j-1)*n_count+1:j*n_count));
    Fup_forw_n(j)= mean(einj_Fup_forw_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_Fup_forw_n(j) = mean((einj_Fup_forw_sort(2,(j-1)*n_count+1:j*n_count)-Fup_forw_n(j)).^2)^0.5;
end
einj_forw_n(n) = mean(einj_Fup_forw_sort(1,(n-1)*n_count+1:end));
Fup_forw_n(n)= mean(einj_Fup_forw_sort(2,(n-1)*n_count+1:end));
var_einj_Fup_forw_n(n) = mean((einj_Fup_forw_sort(2,(n-1)*n_count+1:end)-Fup_forw_n(n)).^2)^0.5;

einj_Fup_inv = [einj_inv; Fup_inv];
einj_Fup_inv_sort = sortrows(einj_Fup_inv')';
N = length(einj_inv);
n = floor(N/n_count);
einj_inv_n=zeros(1,n); 
Fup_inv_n=zeros(1,n); 
var_einj_Fup_inv_n=zeros(1,n); 
for j=1:n-1
    einj_inv_n(j) = mean(einj_Fup_inv_sort(1,(j-1)*n_count+1:j*n_count));
    Fup_inv_n(j)= mean(einj_Fup_inv_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_Fup_inv_n(j) = mean((einj_Fup_inv_sort(2,(j-1)*n_count+1:j*n_count)-Fup_inv_n(j)).^2)^0.5;
end
einj_inv_n(n) = mean(einj_Fup_inv_sort(1,(n-1)*n_count+1:end));
Fup_inv_n(n)= mean(einj_Fup_inv_sort(2,(n-1)*n_count+1:end));
var_einj_Fup_inv_n(n) = mean((einj_Fup_inv_sort(2,(n-1)*n_count+1:end)-Fup_inv_n(n)).^2)^0.5;

einj_Fup_ndual = [einj_ndual; Fup_ndual];
einj_Fup_ndual_sort = sortrows(einj_Fup_ndual')';
N = length(einj_ndual);
n = floor(N/n_count);
einj_ndual_n=zeros(1,n); 
Fup_ndual_n=zeros(1,n); 
var_einj_Fup_ndual_n=zeros(1,n); 
for j=1:n-1
    einj_ndual_n(j) = mean(einj_Fup_ndual_sort(1,(j-1)*n_count+1:j*n_count));
    Fup_ndual_n(j)= mean(einj_Fup_ndual_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_Fup_ndual_n(j) = mean((einj_Fup_ndual_sort(2,(j-1)*n_count+1:j*n_count)-Fup_ndual_n(j)).^2)^0.5;
end
einj_ndual_n(n) = mean(einj_Fup_ndual_sort(1,(n-1)*n_count+1:end));
Fup_ndual_n(n)= mean(einj_Fup_ndual_sort(2,(n-1)*n_count+1:end));
var_einj_Fup_ndual_n(n) = mean((einj_Fup_ndual_sort(2,(n-1)*n_count+1:end)-Fup_ndual_n(n)).^2)^0.5;

figure
plot(einj_dual_n,Fup_dual_n,'o','linewidth',1.3)
hold on
plot(einj_forw_n,Fup_forw_n,'o','linewidth',1.3)
plot(einj_inv_n,Fup_inv_n,'o','linewidth',1.3)
plot(einj_ndual_n,Fup_ndual_n,'o','linewidth',1.3)
patch([einj_dual_n flip(einj_dual_n)], [Fup_dual_n-var_einj_Fup_dual_n flip(Fup_dual_n+var_einj_Fup_dual_n)],[0 0.4470 0.7410], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([einj_forw_n flip(einj_forw_n)], [Fup_forw_n-var_einj_Fup_forw_n flip(Fup_forw_n+var_einj_Fup_forw_n)],[0.8500 0.3250 0.0980], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([einj_inv_n flip(einj_inv_n)], [Fup_inv_n-var_einj_Fup_inv_n flip(Fup_inv_n+var_einj_Fup_inv_n)],[0.9290 0.6940 0.1250], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([einj_ndual_n flip(einj_ndual_n)], [Fup_ndual_n-var_einj_Fup_ndual_n flip(Fup_ndual_n+var_einj_Fup_ndual_n)],[0.4940 0.1840 0.5560], 'FaceAlpha',0.2, 'EdgeColor','none')
leg=legend('dual (divergent)','forward','inverse','dual (convergent)');
set(leg,'interpreter','latex')
xlabel('$$\epsilon_{inj}$$','interpreter','latex')
ylabel('$$F_{up}$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')

einj_Fdown_dual = [einj_dual; Fdown_dual];
einj_Fdown_dual_sort = sortrows(einj_Fdown_dual')';
N = length(einj_dual);
n = floor(N/n_count);
einj_dual_n=zeros(1,n); 
Fdown_dual_n=zeros(1,n); 
var_einj_Fdown_dual_n=zeros(1,n); 
for j=1:n-1
    einj_dual_n(j) = mean(einj_Fdown_dual_sort(1,(j-1)*n_count+1:j*n_count));
    Fdown_dual_n(j)= mean(einj_Fdown_dual_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_Fdown_dual_n(j) = mean((einj_Fdown_dual_sort(2,(j-1)*n_count+1:j*n_count)-Fdown_dual_n(j)).^2)^0.5;
end
einj_dual_n(n) = mean(einj_Fdown_dual_sort(1,(n-1)*n_count+1:end));
Fdown_dual_n(n)= mean(einj_Fdown_dual_sort(2,(n-1)*n_count+1:end));
var_einj_Fdown_dual_n(n) = mean((einj_Fdown_dual_sort(2,(n-1)*n_count+1:end)-Fdown_dual_n(n)).^2)^0.5;


einj_Fdown_forw = [einj_forw; Fdown_forw];
einj_Fdown_forw_sort = sortrows(einj_Fdown_forw')';
N = length(einj_forw);
n = floor(N/n_count);
einj_forw_n=zeros(1,n); 
Fdown_forw_n=zeros(1,n); 
var_einj_Fdown_forw_n=zeros(1,n); 
for j=1:n-1
    einj_forw_n(j) = mean(einj_Fdown_forw_sort(1,(j-1)*n_count+1:j*n_count));
    Fdown_forw_n(j)= mean(einj_Fdown_forw_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_Fdown_forw_n(j) = mean((einj_Fdown_forw_sort(2,(j-1)*n_count+1:j*n_count)-Fdown_forw_n(j)).^2)^0.5;
end
einj_forw_n(n) = mean(einj_Fdown_forw_sort(1,(n-1)*n_count+1:end));
Fdown_forw_n(n)= mean(einj_Fdown_forw_sort(2,(n-1)*n_count+1:end));
var_einj_Fdown_forw_n(n) = mean((einj_Fdown_forw_sort(2,(n-1)*n_count+1:end)-Fdown_forw_n(n)).^2)^0.5;

einj_Fdown_inv = [einj_inv; Fdown_inv];
einj_Fdown_inv_sort = sortrows(einj_Fdown_inv')';
N = length(einj_inv);
n = floor(N/n_count);
einj_inv_n=zeros(1,n); 
Fdown_inv_n=zeros(1,n); 
var_einj_Fdown_inv_n=zeros(1,n); 
for j=1:n-1
    einj_inv_n(j) = mean(einj_Fdown_inv_sort(1,(j-1)*n_count+1:j*n_count));
    Fdown_inv_n(j)= mean(einj_Fdown_inv_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_Fdown_inv_n(j) = mean((einj_Fdown_inv_sort(2,(j-1)*n_count+1:j*n_count)-Fdown_inv_n(j)).^2)^0.5;
end
einj_inv_n(n) = mean(einj_Fdown_inv_sort(1,(n-1)*n_count+1:end));
Fdown_inv_n(n)= mean(einj_Fdown_inv_sort(2,(n-1)*n_count+1:end));
var_einj_Fdown_inv_n(n) = mean((einj_Fdown_inv_sort(2,(n-1)*n_count+1:end)-Fdown_inv_n(n)).^2)^0.5;

einj_Fdown_ndual = [einj_ndual; Fdown_ndual];
einj_Fdown_ndual_sort = sortrows(einj_Fdown_ndual')';
N = length(einj_ndual);
n = floor(N/n_count);
einj_ndual_n=zeros(1,n); 
Fdown_ndual_n=zeros(1,n); 
var_einj_Fdown_ndual_n=zeros(1,n); 
for j=1:n-1
    einj_ndual_n(j) = mean(einj_Fdown_ndual_sort(1,(j-1)*n_count+1:j*n_count));
    Fdown_ndual_n(j)= mean(einj_Fdown_ndual_sort(2,(j-1)*n_count+1:j*n_count));
    var_einj_Fdown_ndual_n(j) = mean((einj_Fdown_ndual_sort(2,(j-1)*n_count+1:j*n_count)-Fdown_ndual_n(j)).^2)^0.5;
end
einj_ndual_n(n) = mean(einj_Fdown_ndual_sort(1,(n-1)*n_count+1:end));
Fdown_ndual_n(n)= mean(einj_Fdown_ndual_sort(2,(n-1)*n_count+1:end));
var_einj_Fdown_ndual_n(n) = mean((einj_Fdown_ndual_sort(2,(n-1)*n_count+1:end)-Fdown_ndual_n(n)).^2)^0.5;

figure
plot(einj_dual_n,Fdown_dual_n,'o','linewidth',1.3)
hold on
plot(einj_forw_n,Fdown_forw_n,'o','linewidth',1.3)
plot(einj_inv_n,Fdown_inv_n,'o','linewidth',1.3)
plot(einj_ndual_n,Fdown_ndual_n,'o','linewidth',1.3)
patch([einj_dual_n flip(einj_dual_n)], [Fdown_dual_n-var_einj_Fdown_dual_n flip(Fdown_dual_n+var_einj_Fdown_dual_n)],[0 0.4470 0.7410], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([einj_forw_n flip(einj_forw_n)], [Fdown_forw_n-var_einj_Fdown_forw_n flip(Fdown_forw_n+var_einj_Fdown_forw_n)],[0.8500 0.3250 0.0980], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([einj_inv_n flip(einj_inv_n)], [Fdown_inv_n-var_einj_Fdown_inv_n flip(Fdown_inv_n+var_einj_Fdown_inv_n)],[0.9290 0.6940 0.1250], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([einj_ndual_n flip(einj_ndual_n)], [Fdown_ndual_n-var_einj_Fdown_ndual_n flip(Fdown_ndual_n+var_einj_Fdown_ndual_n)],[0.4940 0.1840 0.5560], 'FaceAlpha',0.2, 'EdgeColor','none')
leg=legend('dual (divergent)','forward','inverse','dual (convergent)');
set(leg,'interpreter','latex')
xlabel('$$\epsilon_{inj}$$','interpreter','latex')
ylabel('$$F_{down}$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')

edis_Fup_dual = [edis_dual; Fup_dual];
edis_Fup_dual_sort = sortrows(edis_Fup_dual')';
N = length(edis_dual);
n = floor(N/n_count);
edis_dual_n=zeros(1,n); 
Fup_dual_n=zeros(1,n); 
var_edis_Fup_dual_n=zeros(1,n); 
for j=1:n-1
    edis_dual_n(j) = mean(edis_Fup_dual_sort(1,(j-1)*n_count+1:j*n_count));
    Fup_dual_n(j)= mean(edis_Fup_dual_sort(2,(j-1)*n_count+1:j*n_count));
    var_edis_Fup_dual_n(j) = mean((edis_Fup_dual_sort(2,(j-1)*n_count+1:j*n_count)-Fup_dual_n(j)).^2)^0.5;
end
edis_dual_n(n) = mean(edis_Fup_dual_sort(1,(n-1)*n_count+1:end));
Fup_dual_n(n)= mean(edis_Fup_dual_sort(2,(n-1)*n_count+1:end));
var_edis_Fup_dual_n(n) = mean((edis_Fup_dual_sort(2,(n-1)*n_count+1:end)-Fup_dual_n(n)).^2)^0.5;


edis_Fup_forw = [edis_forw; Fup_forw];
edis_Fup_forw_sort = sortrows(edis_Fup_forw')';
N = length(edis_forw);
n = floor(N/n_count);
edis_forw_n=zeros(1,n); 
Fup_forw_n=zeros(1,n); 
var_edis_Fup_forw_n=zeros(1,n); 
for j=1:n-1
    edis_forw_n(j) = mean(edis_Fup_forw_sort(1,(j-1)*n_count+1:j*n_count));
    Fup_forw_n(j)= mean(edis_Fup_forw_sort(2,(j-1)*n_count+1:j*n_count));
    var_edis_Fup_forw_n(j) = mean((edis_Fup_forw_sort(2,(j-1)*n_count+1:j*n_count)-Fup_forw_n(j)).^2)^0.5;
end
edis_forw_n(n) = mean(edis_Fup_forw_sort(1,(n-1)*n_count+1:end));
Fup_forw_n(n)= mean(edis_Fup_forw_sort(2,(n-1)*n_count+1:end));
var_edis_Fup_forw_n(n) = mean((edis_Fup_forw_sort(2,(n-1)*n_count+1:end)-Fup_forw_n(n)).^2)^0.5;

edis_Fup_inv = [edis_inv; Fup_inv];
edis_Fup_inv_sort = sortrows(edis_Fup_inv')';
N = length(edis_inv);
n = floor(N/n_count);
edis_inv_n=zeros(1,n); 
Fup_inv_n=zeros(1,n); 
var_edis_Fup_inv_n=zeros(1,n); 
for j=1:n-1
    edis_inv_n(j) = mean(edis_Fup_inv_sort(1,(j-1)*n_count+1:j*n_count));
    Fup_inv_n(j)= mean(edis_Fup_inv_sort(2,(j-1)*n_count+1:j*n_count));
    var_edis_Fup_inv_n(j) = mean((edis_Fup_inv_sort(2,(j-1)*n_count+1:j*n_count)-Fup_inv_n(j)).^2)^0.5;
end
edis_inv_n(n) = mean(edis_Fup_inv_sort(1,(n-1)*n_count+1:end));
Fup_inv_n(n)= mean(edis_Fup_inv_sort(2,(n-1)*n_count+1:end));
var_edis_Fup_inv_n(n) = mean((edis_Fup_inv_sort(2,(n-1)*n_count+1:end)-Fup_inv_n(n)).^2)^0.5;

edis_Fup_ndual = [edis_ndual; Fup_ndual];
edis_Fup_ndual_sort = sortrows(edis_Fup_ndual')';
N = length(edis_ndual);
n = floor(N/n_count);
edis_ndual_n=zeros(1,n); 
Fup_ndual_n=zeros(1,n); 
var_edis_Fup_ndual_n=zeros(1,n); 
for j=1:n-1
    edis_ndual_n(j) = mean(edis_Fup_ndual_sort(1,(j-1)*n_count+1:j*n_count));
    Fup_ndual_n(j)= mean(edis_Fup_ndual_sort(2,(j-1)*n_count+1:j*n_count));
    var_edis_Fup_ndual_n(j) = mean((edis_Fup_ndual_sort(2,(j-1)*n_count+1:j*n_count)-Fup_ndual_n(j)).^2)^0.5;
end
edis_ndual_n(n) = mean(edis_Fup_ndual_sort(1,(n-1)*n_count+1:end));
Fup_ndual_n(n)= mean(edis_Fup_ndual_sort(2,(n-1)*n_count+1:end));
var_edis_Fup_ndual_n(n) = mean((edis_Fup_ndual_sort(2,(n-1)*n_count+1:end)-Fup_ndual_n(n)).^2)^0.5;

figure
plot(-edis_dual_n,Fup_dual_n,'o','linewidth',1.3)
hold on
plot(-edis_forw_n,Fup_forw_n,'o','linewidth',1.3)
plot(-edis_inv_n,Fup_inv_n,'o','linewidth',1.3)
plot(-edis_ndual_n,Fup_ndual_n,'o','linewidth',1.3)
patch([-edis_dual_n flip(-edis_dual_n)], [Fup_dual_n-var_edis_Fup_dual_n flip(Fup_dual_n+var_edis_Fup_dual_n)],[0 0.4470 0.7410], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([-edis_forw_n flip(-edis_forw_n)], [Fup_forw_n-var_edis_Fup_forw_n flip(Fup_forw_n+var_edis_Fup_forw_n)],[0.8500 0.3250 0.0980], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([-edis_inv_n flip(-edis_inv_n)], [Fup_inv_n-var_edis_Fup_inv_n flip(Fup_inv_n+var_edis_Fup_inv_n)],[0.9290 0.6940 0.1250], 'FaceAlpha',0.2, 'EdgeColor','none')
patch([-edis_ndual_n flip(-edis_ndual_n)], [Fup_ndual_n-var_edis_Fup_ndual_n flip(Fup_ndual_n+var_edis_Fup_ndual_n)],[0.4940 0.1840 0.5560], 'FaceAlpha',0.2, 'EdgeColor','none')
leg=legend('dual (divergent)','forward','inverse','dual (convergent)');
set(leg,'interpreter','latex')
xlabel('$$\epsilon_{diss}$$','interpreter','latex')
ylabel('$$F_{up}$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e_Fup_dual = [einj_dual+edis_dual; Fup_dual];
e_Fup_dual_sort = sortrows(e_Fup_dual')';
N = length(einj_dual);
n = floor(N/n_count);
e_dual_n=zeros(1,n); 
Fup_dual_n=zeros(1,n); 
var_e_Fup_dual_n=zeros(1,n); 
for j=1:n-1
    e_dual_n(j) = mean(e_Fup_dual_sort(1,(j-1)*n_count+1:j*n_count));
    Fup_dual_n(j)= mean(e_Fup_dual_sort(2,(j-1)*n_count+1:j*n_count));
    var_e_Fup_dual_n(j) = mean((e_Fup_dual_sort(2,(j-1)*n_count+1:j*n_count)-Fup_dual_n(j)).^2)^0.5;
end
e_dual_n(n) = mean(e_Fup_dual_sort(1,(n-1)*n_count+1:end));
Fup_dual_n(n)= mean(e_Fup_dual_sort(2,(n-1)*n_count+1:end));
var_e_Fup_dual_n(n) = mean((e_Fup_dual_sort(2,(n-1)*n_count+1:end)-Fup_dual_n(n)).^2)^0.5;

figure
plot(1e9*e_dual_n,1e9*Fup_dual_n,'o','linewidth',1.3)
hold on
patch(1e9*[e_dual_n flip(e_dual_n)], 1e9*[Fup_dual_n-var_e_Fup_dual_n flip(Fup_dual_n+var_e_Fup_dual_n)],[0 0.4470 0.7410], 'FaceAlpha',0.2, 'EdgeColor','none')
leg=legend('$$F_{large}$$');
set(leg,'interpreter','latex')
xlabel('$$\epsilon_{inj}-\epsilon_{dis}$$','interpreter','latex')
ylabel('$$F$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')

e_Fdown_dual = [einj_dual+edis_dual; Fdown_dual];
e_Fdown_dual_sort = sortrows(e_Fdown_dual')';
N = length(einj_dual);
n = floor(N/n_count);
e_dual_n=zeros(1,n); 
Fdown_dual_n=zeros(1,n); 
var_e_Fdown_dual_n=zeros(1,n); 
for j=1:n-1
    e_dual_n(j) = mean(e_Fdown_dual_sort(1,(j-1)*n_count+1:j*n_count));
    Fdown_dual_n(j)= mean(e_Fdown_dual_sort(2,(j-1)*n_count+1:j*n_count));
    var_e_Fdown_dual_n(j) = mean((e_Fdown_dual_sort(2,(j-1)*n_count+1:j*n_count)-Fdown_dual_n(j)).^2)^0.5;
end
e_dual_n(n) = mean(e_Fdown_dual_sort(1,(n-1)*n_count+1:end));
Fdown_dual_n(n)= mean(e_Fdown_dual_sort(2,(n-1)*n_count+1:end));
var_e_Fdown_dual_n(n) = mean((e_Fdown_dual_sort(2,(n-1)*n_count+1:end)-Fdown_dual_n(n)).^2)^0.5;

figure
plot(1e9*e_dual_n,1e9*Fdown_dual_n,'o','linewidth',1.3)
hold on
patch([1e9*e_dual_n flip(1e9*e_dual_n)], 1e9*[Fdown_dual_n-var_e_Fdown_dual_n flip(Fdown_dual_n+var_e_Fdown_dual_n)],[0.8500 0.3250 0.0980], 'FaceAlpha',0.2, 'EdgeColor','none')
leg=legend('$$F_{small}$$');
set(leg,'interpreter','latex')
xlabel('$$\epsilon_{inj}-\epsilon_{dis}$$','interpreter','latex')
ylabel('$$F$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')


figure
plot(1e9*e_dual_n,Fup_dual_n./Fdown_dual_n,'o','linewidth',1.3)
hold on
plot(1e9*e_dual_n,ones(1,length(e_dual_n)),'k--','linewidth',1.3)
patch([1e9*e_dual_n flip(1e9*e_dual_n)], [Fup_dual_n./Fdown_dual_n.*(1-var_e_Fup_dual_n./Fup_dual_n-var_e_Fdown_dual_n./Fdown_dual_n) flip(Fup_dual_n./Fdown_dual_n.*(1+var_e_Fup_dual_n./Fup_dual_n+var_e_Fdown_dual_n./Fdown_dual_n))],[0 0.4470 0.7410], 'FaceAlpha',0.2, 'EdgeColor','none')
% leg=legend('$$F_{small}$$');
% set(leg,'interpreter','latex')
xlabel('$$\epsilon_{inj}-\epsilon_{dis}\,(Jkg^{-1}s^{-1})$$','interpreter','latex')
ylabel('$$F_{large}/F_{small}$$','interpreter','latex')
set(gca,'FontSize',14,'FontName','Times')