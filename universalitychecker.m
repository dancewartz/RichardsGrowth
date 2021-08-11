%% Defining parameters
clear; close all; clc;
mu = 1; 
a = mu/1000; 
K = mu/a;
D = 0.6; 
sigma = sqrt(2);
cD = 2*D/sigma^2;
N = 2^20; 
y0 = mu/a/4000;
ymax = K/5;
Tg = 0.1/mu;
Te = 3/D;
dt = 0.01/D;
sdt = sqrt(dt);
tse = dt:dt:Te;
num_cycles = 175;
ts = (Te+Tg)*[1:num_cycles]; 
measure = Tg + Te;
measure_ind = 1;

y = y0*ones(1,N);
m1 = NaN*ones(1,num_cycles);
m2 = NaN*ones(1,num_cycles);
m = y0;

CDFinv = @(x,m) cD*m./gammaincinv(x,1+cD,'upper');
growthOdeFunc = @(t,y) mu*y.*(1-y/K).*(1 + 0.2*y/K + 0.3*(y/K).^2);

%% Simulation Loop
for j = 1:num_cycles
    j
    for k = 1:length(tse)
        dW = sdt*randn(size(y));
        y = exp(-sigma^2/2*dt + sigma*dW).*y; % Seascape Noise, Geometric Brownian Motion
        y = y + dt*D*(m - y); % Diffusion
    end
    
    y = y(y<ymax); % Cutoff
    
    [~,y_growth] = ode45(growthOdeFunc,[0,Tg],y); % Growth, uses RK4
    
    y = y_growth(end,:);
    
    m = mean(y);
    m1(j) = m;
    m2(j) = mean(y.^2);
end

%% Plotting
figure
set(gca,'fontsize',15)
hold on
xlabel('Time')
ylabel('\langle y \rangle')
plot(ts,m1,'.','MarkerSize',20)

figure
hold on
set(gca,'xscale','log','yscale','log')
plot(m1,m2,'.','MarkerSize',20)
p = polyfit(log(m1(ts>measure)),log(m2(ts>measure)),1)
plot(m1,exp(polyval(p,log(m1))),'--','LineWidth',4)

figure(1)
% odefunc = @(t,x) mu*x - a*exp(polyval(p,log(x)));
odefunc = @(t,x) Tg/(Te + Tg)*mu*x.*(1 - (x/mean(m1(end))).^cD);
[tr,mr] = ode45(odefunc,[0 ts(end)],y0);
plot(tr,mr,'--','LineWidth',3,'DisplayName','Fitted Richards')
saveas(gca,'Figures/universal.png')

m = y0;
test = CDFinv(rand(1,N),m);
test = test(test<ymax);
[~,logedges] = histcounts(log10(test));
[counts,edges] = histcounts(test,10.^logedges,'Normalization','pdf');
midpoints = 1/2*(edges(1:end-1) + edges(2:end));
theo_dist = @(y) (2*D/sigma^2*m)^(2*D/sigma^2 + 1)/gamma(2*D/sigma^2 + 1)*y.^(-2-2*D/sigma^2).*exp(-2*D/sigma^2*m./y);
figure
hold on
legend
xlabel('y')
ylabel('P(y)')
set(gca,'xscale','log','yscale','log');
plot(midpoints,counts,'.','markersize',20,'DisplayName','P(y,0)')
plot(midpoints,theo_dist(midpoints),'--','LineWidth',3,'DisplayName','Theory')

% writematrix([ts' m1'],'FigData/Universal.csv')
% T = table(mu,a,D,sigma,Te,Tg,N,y0,ymax,num_cycles);
% writetable(T,'FigData/UniversalParams.csv')

