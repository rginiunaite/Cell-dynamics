% 1d ODE system for one-point distribution function coupled with two-point
% distribution functions for all distances
% Runge-Kutta method
clear all

%%%%%%%%%%%%%%
%%%% Parameters inside the function
 L = 100; % length of the domain,change inside the function too!!!!!!!!!!!
 
 Pm = 1; % transition rate per unit time of moving to another lattice site
 
 Pp = 0.5; % proliferation rate per unit time of giving rise to another agent
 
 Pd = 0.1; % death rate per unit time


cells = zeros(L,1); % array that will store the sites, this is just for 1d

% initially choose randomly if each state is occupied or no

% choose a random number for each state
rand_initial = rand(L,1);

for i = 1:L      
   % if the random number is greater than 0.5, then assume that initially
   % there is a cell at that site. 1 corresponds to that cell
   if rand_initial(i)>0.5
       cells(i) = 1; 
   end
end


Q = sum(cells); % store the number of cells at different times

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is for logistic growth

Q_initial = Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cA0 = Q/L; % initial density of cells


N = 10000; % number of simulations

% Initial conditions 

y  = zeros (L+1,N);

y(:,1) = cA0^2 * ones(L+1,1); % initial density of two-point distribution functions 
y(1,1) = cA0; % initial condition for 1-point distribution function


t = 0; % start time

dt = 0.01; % time step


for i = 2 : N
    
    t = t + dt;

    k1 = dynamics (t,y(:,i-1));
    k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1);
    k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2);
    k4 = dynamics (t+dt, y(:,i-1) + dt*k3);
    y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  
    
    store_time(i) = t;
end


rescaled_time = (Pp-Pd)*store_time; % rescale time to allow parameter comparison

y_rescaled = (Pp-Pd)/Pp * y;


figure
plot(rescaled_time,y_rescaled(1,:),'LineWidth', 2)
%h_legend = legend('Gillespie 1D','Logistic growth','Total');
%set(h_legend,'FontSize',14)
xlabel('time','FontSize',14)
ylabel('density','FontSize',14)




%%%%%%%%%%%%%%%%
% solution to the corresponding logistic growth equation

cA0 = Q_initial/L; % initial density
cA0_bar = (Pp - Pd)/(Pp)*cA0;

f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 

for i = 1 : length(store_time)
    %store_time(i) = (Pp - Pd) * store_time(i); % rescale time
    y_logistic(i) = f(rescaled_time(i)); 
end

%%% plot logistic solution 

hold on
plot(rescaled_time,y_logistic,'LineWidth', 2)
h_legend = legend('KSA','Logistic growth');
set(h_legend,'FontSize',14)
xlabel('time','FontSize',14)
ylabel('density','FontSize',14)




function deriv = dynamics(t,y)
    
    L = 100; % length of the domain

    Pm = 1; % transition rate per unit time of moving to another lattice site

    Pp = 0.5; % proliferation rate per unit time of giving rise to another agent

    Pd = 0.1; % death rate per unit time


    deriv = zeros(L+1,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*(y(3) - y(2)) - 2* Pd * y(2) + ...
        Pp * (y(1) - y(2)) + Pp * (y(3)*(y(3)*(y(1)-y(2))^2)/(y(1)^2*(1-y(1))))); % dynamics of two-point distribution function (distance 1)
    
 % dynamics of all two-point functions from distance 2 to L-1
    for i = 3:L
        deriv(i) = (Pm*(y(i-1) + y(i+1) - 2* y(i)) - 2 * Pd * y(i) +...
            Pp * ((y(1)-y(i))*(y(1)-y(2))*(y(i-1)+y(i+1)))/(y(1)^2*(1-y(1))));
    end

    deriv(L+1) = (Pm*(y(L) + y(2) - 2* y(L+1)) - 2 * Pd * y(L+1) +...
            Pp * ((y(1)-y(L+1))*(y(1)-y(2))*(y(L)+y(2)))/(y(1)^2*(1-y(1)))); % for the last one, since we have
        % periodic boundary condition L+1 neighbour is 1 neighbour, i.e. y(2)
    
end