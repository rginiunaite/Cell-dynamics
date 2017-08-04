% combination of Gillespie (40 realisations), logistic and Runge-Kutta with power-2 approximation in 1D

clear all

% Parameters

% lattice dimension
dim = 1;

% length of the domain
L = 100;

Pm = 1; % transition rate per unit time of moving to another lattice site

Pp = 1; % proliferation rate per unit time of giving rise to another agent

Pd = 0; % death rate per unit time

TimeReductionFactor = 2; % when death rate is increase, this has to be increased as well


S = load('gillespie1D.mat'); % load the data for a particular choice of Pp and Pd, note that L=100 is fixed

figure 

plot(S.save_time,S.save_density,'LineWidth', 3);set(gca,'FontSize',14); % number of cells in the system
%h_legend = legend('1st component','2nd component','Total');
%set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)


%%%%%%%%%%%%%%%%
% solution to the corresponding logistic growth equation

%cA0 = Q_initial/L; % initial density
cA0_bar = (Pp - Pd)/(Pp)*S.save_cA0;

f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 

for i = 1 : length(S.save_time)
       y(i) = f(S.save_time(i)); 
end

%%% plot logistic solution 

hold on
plot(S.save_time,y,'LineWidth', 3)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)

dt = 0.01; % time step

N = round(S.save_store_time(length(1:round(length(S.save_store_time)/TimeReductionFactor)))/dt); % number of simulations, chosen to be the same as for Gillespie

% Initial conditions 

y  = zeros (L+1,N);

y(:,1) = S.save_cA0^2 * ones(L+1,1); % initial density of two-point distribution functions 
y(1,1) = S.save_cA0; % initial condition for 1-point distribution function


t = 0; % start time


for i = 2 : N
    
    t = t + dt;

    k1 = dynamics (t,y(:,i-1),L,Pm,Pp,Pd);
    k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,L,Pm,Pp,Pd);
    k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,L,Pm,Pp,Pd);
    k4 = dynamics (t+dt, y(:,i-1) + dt*k3,L,Pm,Pp,Pd);
    y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  
    
    store_time_Runge(i) = t;
end


rescaled_time = (Pp-Pd)*store_time_Runge; % rescale time to allow parameter comparison

y_rescaled = (Pp-Pd)/Pp * y;


hold on
plot(rescaled_time,y_rescaled(1,:),'LineWidth', 3)
h_legend = legend('Gillespie 1D','Logistic growth','Power 2');
set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)
title(['L = ' num2str(L) ', P_m = ' num2str(Pm) ', P_p = ' num2str(Pp) ', P_d = ' num2str(Pd) ' '],'FontSize',14)
set(gca,'linewidth',3)
set(gca,'FontWeight','bold')
set(gca,'FontSize',24)


function deriv = dynamics(t,y,L,Pm,Pp,Pd)
    

    deriv = zeros(L+1,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*(y(3) - y(2)) - 2* Pd * y(2) + ...
        Pp * (y(1) - y(2)) + Pp * (0.5*(2*y(3)*(y(1) - y(2))/y(1) + (y(1)-y(2))^2/(1-y(1))) -y(1)^2*(1-y(1)))); % dynamics of two-point distribution function (distance 1)
    
 % dynamics of all two-point functions from distance 2 to L-1
    for i = 3:L
        deriv(i) = (Pm*(y(i-1) + y(i+1) - 2* y(i)) - 2 * Pd * y(i) +...
            Pp * 0.5*((y(1)-y(i))*(y(i+1) + y(i-1))/y(1) + 2 * (y(1)-y(i))*(y(1)-y(2))/(1-y(1)) + (y(i-1)+y(i+1))*(y(1)-y(2))/y(1) - y(1)^2*(1-y(1))));
    end

    deriv(L+1) = (Pm*(y(L) + y(2) - 2* y(L+1)) - 2 * Pd * y(L+1) +...
            Pp * 0.5*((y(1)-y(L+1))*(y(1) + y(L))/y(1) + 2 * (y(1)-y(L+1))*(y(1)-y(2))/(1-y(1)) + (y(L)+y(1))*(y(1)-y(2))/y(1) - y(1)^2*(1-y(1)))); % for the last one, since we have
        % periodic boundary condition L+1 neighbour is 1 neighbour, i.e. y(2)
    
end

