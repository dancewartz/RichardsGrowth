%% Defining parameters
clear; close all; clc;
mu = 1; 
a = mu/1000; 
D = 0.5; 
sigma = sqrt(2);
cD = 2*D/sigma^2;
N = 2^20; 
y0 = mu/a/5000;
Tg = 0.1/mu;
Te = 3/D;
num_cycles = 175;
ts = (Te+Tg)*[1:num_cycles]; 
measure = Te + Tg; 
measure_ind = 1;

ymax = mu/a/2;

dt = 1e-2/sigma^2;
sdt = sqrt(dt);
tse = dt:dt:Te;

y = y0*ones(1,N);
m = y0;

%% Simulation loop
for j = 1:num_cycles
    j
    for k = 1:length(tse)
        dW = sdt*randn(size(y));
        y = exp(-sigma^2/2*dt + sigma*dW).*y; % Seascape Noise - Geometric Brownian Motion
        y = y + dt*D*(m - y); % Diffusion
    end
    
    y(y>ymax) = []; % Cutoff
    
    y = exp(mu*Tg)*y*mu./(-a*y + a*exp(mu*Tg)*y + mu); % growth
  
    m1(j) = mean(y);
    m2(j) = mean(y.^2);
    m = mean(y);
end
%% Plotting
figure
set(gca,'fontsize',12)
hold on
xlabel('t')
ylabel('m_1')
plot(ts,m1,'.','MarkerSize',20)
% plot(ts,m1,'LineWidth',3)

figure
legend
hold on
set(gca,'xscale','log','yscale','log')
plot(m1,m2,'.','MarkerSize',20)
p = polyfit(log(m1(ts>measure)),log(m2(ts>measure)),1)
plot(m1,exp(polyval(p,log(m1))),'--','LineWidth',3,'DisplayName',['\gamma = ' num2str(p(1))])

figure(1)
% odefunc = @(t,x) mu*x - a*exp(polyval(p,log(x)));
odefunc = @(t,x) Tg/(Te+Tg)*mu*x.*(1 - (x/mean(m1(end - round(end/10):end))).^cD);
[tr,mr] = ode45(odefunc,[0 ts(end)],y0);
plot(tr,mr,'--','LineWidth',3,'DisplayName','Fitted Richards')

% writematrix([ts' m1'],'FigData/TwoPhase.csv')
% T = table(mu,a,D,sigma,Te,Tg,N,y0,ymax,num_cycles);
% writetable(T,'FigData/TwoPhaseParams.csv')

% saveas(gca,'Figures/seasonalD0.2.png')
% 
% figure(2)
% saveas(gca,'Figures/m1m2D0.2.png')