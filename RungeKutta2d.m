% 2d ODE system for one-point distribution function coupled with two-point
% distribution functions for all distances
% Runge-Kutta method
clear all

%%%%%%%%%%%%%%
%%%% Parameters inside the function


% lattice dimension
dim = 2;

% length of the side of square
L = 10;

max_rad = sqrt(2*(L/2)^2); % maximal radius that would cover all the lattice

Pm = 1; % transition rate per unit time of moving to another lattice site

Pp = 0.5; % proliferation rate per unit time of giving rise to another agent

Pd = 0.1; % death rate per unit time


cells = zeros(L,L); % array that will store the sites

% initially choose randomly if each state is occupied or no

% choose a random number for each state, only 5% occupied
rand_initial = rand(L,L)/1.9;

for i = 1:L      
   % if the random number is greater than 0.5, then assume that initially
   % there is a cell at that site. 1 corresponds to that cell
   for j = 1:L
      if rand_initial(i,j)>0.5
        cells(i,j) = 1; 
      end
   end
end



Q = sum(sum(cells)); % store the number of cells at different times


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is for logistic growth

Q_initial = Q;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cA0 = Q/L; % initial density of cells

N = 10000; % number of simulations

% Initial conditions 

dr = 0.5; % radially approximated at this distance (after distance 5), i.e 5, 5 + dr, 5 + 2dr

numberDer = (round(max_rad) - 5)/dr ; % this is to cound the number of additional derivatives
     
y = zeros(14 + numberDer,1); % 14 fixed derivatives, the rest varying depending on the size of the grid


y(:,1) = cA0^2 * ones(14+numberDer,1); % initial density of two-point distribution functions 
y(1,1) = cA0; % initial condition for 1-point distribution function


t = 0; % start time

dt = 0.01; % time step


for i = 2 : N
    
    t = t + dt;

    k1 = dynamics (t,y(:,i-1),Pm,Pp,Pd,max_rad);
    k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,Pm,Pp,Pd,max_rad);
    k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,Pm,Pp,Pd,max_rad);
    k4 = dynamics (t+dt, y(:,i-1) + dt*k3,Pm,Pp,Pd,max_rad);
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
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)




function deriv = dynamics(t,y,Pm,Pp,Pd,max_rad)
    
    dr = 0.5; % radially approximated at this distance (after distance 5), i.e 5, 5 + dr, 5 + 2dr

    numberDer = (round(max_rad) - 5)/dr ; % this is to cound the number of additional derivatives
     
    deriv = zeros(14 + numberDer,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*0.5* (y(4) + 2 * y(3) - 3* y(2)) - 2 * Pd*y(2) + ...
        Pp * 0.5 * (y(1) - y(2)) + 0.5 * Pp * (y(3) + y(4)) *(y(1) - y(2))^2 /(y(1)^2*(1-y(1))));
    
    deriv(3) = 0.5 * Pm * (y(2) + y(5) + y(2) + y(5) - 4*y(3)) - ...
        2 * Pd * y(3) + 0.5 * Pp * (y(1) - y(3))* (y(1) - y(2))*(y(2)+y(5)+y(2)+ y(5))/(y(1)^2*(1-y(1)));
    
    deriv(4) = 0.5 * Pm * (y(7) + y(2) + y(5) + y(5) - 4*y(4)) - ...
        2 * Pd * y(4) + 0.5 * Pp * (y(1) - y(4))* (y(1) - y(2))*(y(7)+y(2)+y(5)+ y(5))/(y(1)^2*(1-y(1)));
    
    deriv(5) = 0.5 * Pm * (y(8) + y(3) + y(6) + y(4) - 4*y(5)) - ...
        2 * Pd * y(5) + 0.5 * Pp * (y(1) - y(5))* (y(1) - y(2))*(y(8)+y(3)+y(6)+ y(4))/(y(1)^2*(1-y(1)));
    
    deriv(6) = 0.5 * Pm * (y(9) + y(5) + y(9) + y(5) - 4*y(6)) - ...
        2 * Pd * y(6) + 0.5 * Pp * (y(1) - y(6))* (y(1) - y(2))*(y(9)+y(5)+y(9)+ y(5))/(y(1)^2*(1-y(1))); 
    
    deriv(7) = 0.5 * Pm * (y(10) + y(4) + y(8) + y(8) - 4*y(7)) - ...
        2 * Pd * y(7) + 0.5 * Pp * (y(1) - y(7))* (y(1) - y(2))*(y(10)+y(4)+y(8)+ y(8))/(y(1)^2*(1-y(1)));     
    
    deriv(8) = 0.5 * Pm * (y(11) + y(9) + y(5) + y(7) - 4*y(8)) - ...
        2 * Pd * y(8) + 0.5 * Pp * (y(1) - y(8))* (y(1) - y(2))*(y(11)+y(9)+y(5)+ y(7))/(y(1)^2*(1-y(1)));
    
    deriv(9) = 0.5 * Pm * (y(12) + y(8) + y(6) + y(13) - 4*y(9)) - ...
        2 * Pd * y(9) + 0.5 * Pp * (y(1) - y(9))* (y(1) - y(2))*(y(12)+y(8)+y(6)+ y(13))/(y(1)^2*(1-y(1)));   
    
    deriv(10) = 0.5 * Pm * (y(7) + y(14) + y(11) + y(11) - 4*y(10)) - ...
        2 * Pd * y(10) + 0.5 * Pp * (y(1) - y(10))* (y(1) - y(2))*(y(7)+y(14)+y(11)+ y(11))/(y(1)^2*(1-y(1)));
    
    deriv(11) = 0.5 * Pm * (y(13) + y(10) + y(8) + y(15) - 4*y(11)) - ...
        2 * Pd * y(11) + 0.5 * Pp * (y(1) - y(11))* (y(1) - y(2))*(y(13)+y(10)+y(8)+ y(15))/(y(1)^2*(1-y(1)));
    
    deriv(12) = 0.5 * Pm * (y(9) + y(9) + y(14) + y(14) - 4*y(12)) - ...
        2 * Pd * y(12) + 0.5 * Pp * (y(1) - y(12))* (y(1) - y(2))*(y(9)+y(9)+y(14)+ y(14))/(y(1)^2*(1-y(1)));  
    
    deriv(13) = 0.5 * Pm * (y(11) + y(14) + y(9) + y(15) - 4*y(13)) - ...
        2 * Pd * y(13) + 0.5 * Pp * (y(1) - y(13))* (y(1) - y(2))*(y(11)+y(14)+y(9)+ y(15))/(y(1)^2*(1-y(1))); 
    
    deriv(14) = 0.5 * Pm * (y(13) + y(12) + y(16) + y(16) - 4*y(14)) - ...
        2 * Pd * y(14) + 0.5 * Pp * (y(1) - y(14))* (y(1) - y(2))*(y(13)+y(12)+y(16)+ y(16))/(y(1)^2*(1-y(1)));
    
    
    % now the distance is more than 5, so all will be approximated
    % uniformly 
    
   
    
    current_radius = 5;
    
    i = 15; % will start considering regular grid from this derivative
    
    % update derivatives until 1 before maximal radius
    while current_radius +dr < max_rad -dr % this is not correct
        
        current_radius = current_radius + dr; % increase radius by dr
        
        deriv(i) = Pm * (y(i-1)+ y(i+1) - 2 *y(i))/dr^2 + 1/current_radius * (y(i+1) - y(i-1))/(2*dr) ...
            - 2 * Pd * y(i) + Pp * (y(1) - y(i))* (y(1) - y(2))*(y(i+1)+ y(i-1))/(y(1)^2*(1-y(1)));
        i = i+1;
    end
    
    last = 14 + numberDer;
    %periodic boundary conditions
    deriv(last) =  Pm * (y(last-1)+ y(2) - 2 *y(last))/dr^2 + 1/current_radius * (y(2) - y(last-1))/(2*dr) ...
            - 2 * Pd * y(last) + Pp * (y(1) - y(last))* (y(1) - y(2))*(y(2)+ y(last-1))/(y(1)^2*(1-y(1)));
    
end






