% HEXAGONAL DOMAIN: combination logistic and Runge-Kutta with KSA in 2D
clear all

% Parameters

% lattice dimension
dim = 2;

% half length of the side of square
L = 100;

Pm = 1; % transition rate per unit time of moving to another lattice site

Pp = 0.5; % proliferation rate per unit time of giving rise to another agent

Pd = 0.1; % death rate per unit time

TimeReductionFactor = 1.5; % when death rate is increase, this has to be increased as well


S = load('gillespie2D_Pp05_Pd01.mat'); % load the data for a particular choice of Pp and Pd, note that L=100 is fixed

% figure 
% 
% plot(S.save_time,S.save_density,'LineWidth', 3);set(gca,'FontSize',14); % number of cells in the system
% %h_legend = legend('1st component','2nd component','Total');
% %set(h_legend,'FontSize',14)
% xlabel('rescaled time','FontSize',14)
% ylabel('rescaled density','FontSize',14)



%%%%%%%%%%%%%%%%
% solution to the corresponding logistic growth equation

%cA0 = Q_initial/L^2; % initial density

cA0_bar = (Pp)/(Pp - Pd)*S.save_cA0;

f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 



for i = 1 : length(S.save_time)
       y_log(i) = f(S.save_time(i)); 
end

%%% plot logistic solution 

% hold on
% plot(S.save_time,y_log,'LineWidth', 3)
% xlabel('rescaled time','FontSize',14)
% ylabel('rescaled density','FontSize',14)


%%%%%%%%%%%%%%%%%%%%%%%%%
% adding KSA


%cA0 = Q_initial/L^2; % initial density of cells

dt = 0.01; % time step

N = round(S.save_store_time(length(1:round(length(S.save_store_time)/TimeReductionFactor)))/dt); % number of simulations, chosen to be the same as for Gillespie

% Additional parameters


y(:,1) = S.save_cA0^2 * ones(L+1,1); % initial density of two-point distribution functions 
y(1,1) = S.save_cA0; % initial condition for 1-point distribution function


t = 0; % start time

dt = 0.01; % time step

% 4th order Runge Kutta

for i = 2 : N
    
    t = t + dt;

    k1 = dynamics (t,y(:,i-1),Pm,Pp,Pd,L);
    k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,Pm,Pp,Pd,L);
    k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,Pm,Pp,Pd,L);
    k4 = dynamics (t+dt, y(:,i-1) + dt*k3,Pm,Pp,Pd,L);
    y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  
    
    store_time_Runge(i) = t;
end


rescaled_time = (Pp-Pd)*store_time_Runge; % rescale time to allow parameter comparison

y_rescaled = Pp/(Pp-Pd) * y;


hold on
plot(rescaled_time,y_rescaled(1,:),'LineWidth', 3)
h_legend = legend('Gillespie 2D','Logistic growth','KSA - square','KSA - hexagonal');
set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)
title(['L = ' num2str(L) ', P_m = ' num2str(Pm) ', P_p = ' num2str(Pp) ', P_d = ' num2str(Pd) ' '],'FontSize',14)
set(gca,'linewidth',3)
set(gca,'FontWeight','bold')
set(gca,'FontSize',24)






function deriv = dynamics(t,y,Pm,Pp,Pd,L)
    
 
     
    deriv = zeros(L+1,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*1/3* (2*y(2) + 3 * y(3) - 5* y(2)) - 2 * Pd*y(2) + ...
        Pp * 1/3 * (y(1) - y(2)) + 1/3 * Pp * (2*y(2) + 3*y(3)) *(y(1) - y(2))^2 /(y(1)^2*(1-y(1)))); % distance 1
    
    
    for n=3:L
    
        deriv(n) = 1/3 * Pm * (y(n-1) + 2 * y(n) +3 * y(n+1)- 6*y(n)) - ...
            2 * Pd * y(n) + 1/3 * Pp * (y(1) - y(n))* (y(1) - y(2))*(y(n-1)+2*y(n)+3*y(n+1))/(y(1)^2*(1-y(1))); % 
    
    end
    
    deriv(L) = 1/3 * Pm * (y(L-1) + 2 * y(L) +3 * y(2)- 6*y(L)) - ...
            2 * Pd * y(L) + 1/3 * Pp * (y(1) - y(L))* (y(1) - y(2))*(y(L-1)+2*y(n)+3*y(2))/(y(1)^2*(1-y(1))); % Periodic boundary conditions
    
    
end



