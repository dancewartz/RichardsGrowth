%% Defining Parameters
clear; close all; clc;
mu = 1; 
a = mu; 
D = 1; 
sigma = 2;
N = 2^10; 
y0 = mu/a/500; 
T = 12/mu;
measure = 4/mu;
dt = 1e-2*min([1 1/mu 1/D]);
% dt = 1e-2;
y = y0*ones(N); 
ts = dt:dt:T; 
sdt = sqrt(dt);
m1 = NaN*ones(size(ts));
m2 = NaN*ones(size(ts));
m = y0; 

%% Simulation Loop
for j = 1:length(ts)
    dW = sdt*randn(N,N);
    y = y + D*dt*(circshift(y,1,1) + circshift(y,-1,1) + circshift(y,1,2) + circshift(y,-1,2) - 4*y); % Diffusion
    y = exp((mu-sigma^2/2)*dt + sigma*dW).*y; % Geometric Brownian Motion
    y = y./(1+a*dt*y); % Nonlinearity
    
    m = mean(y,'all');
    m1(j) = m;
    m2(j) = mean(y.^2,'all');
    
end

%% Plotting
subplot(2,1,1,'FontSize',12)
hold on
legend('location','northwest')
% set(gca,'fontsize',15)
xlabel('t')
ylabel('m_1')
plot(ts,m1,'LineWidth',3,'DisplayName','Simulation')

subplot(2,1,2,'fontsize',12,'xscale','log','yscale','log')
hold on
legend('location','northwest')
% set(gca,'fontsize',15,'xscale','log','yscale','log')
xlabel('m_1')
ylabel('m_2')
plot(m1,m2,'.','markersize',20,'DisplayName','Simulation')
p = polyfit(log(m1(ts>measure)),log(m2(ts>measure)),1)
plot(m1,exp(polyval(p,log(m1))),'--','LineWidth',3,'DIsplayName',['\gamma = ' num2str(p(1))])
% plot(m1,m2(round(measure/dt))/m1(round(measure/dt))^2*m1.^2,'k--','LineWidth',3,'DisplayName','Logistic')


subplot(2,1,1)
odefunc = @(t,x) mu*x - a*exp(polyval(p,log(x)));
[tr,mr] = ode45(odefunc,[measure T],m1(round(measure/dt)));
plot(tr,mr,'--','LineWidth',3,'DisplayName','Fitted Richards')

% odefunc = @(t,x) mu*x - mu*m1(end)^(-1)*x.^2;
% [tl,ml] = ode45(odefunc,[measure T],m1(round(measure/dt)));
% plot(tl,ml,'k--','LineWidth',3,'DisplayName','Fitted Logistic')

% saveas(gca,'Figures/twodim.png')
% writematrix([ts' m1' m2'],'FigData/TwoDim.csv')
% Tab = table(mu,a,D,sigma,N,y0,dt,T);
% writetable(Tab,'FigData/TwoDimParams.csv')