% combination logistic and Runge-Kutta with power 2 approximation in 2D


% Parameters

% lattice dimension
dim = 2;

% length of the side of square
L = 100;

Pm = 1; % transition rate per unit time of moving to another lattice site

Pp = 0.5; % proliferation rate per unit time of giving rise to another agent

Pd = 0; % death rate per unit time

TimeReductionFactor = 4; % when death rate is increase, this has to be increased as well

S = load('gillespie2D_Pp05_Pd0.mat'); % load the data for a particular choice of Pp and Pd, note that L=100 is fixed

figure 

plot(S.save_time,S.save_density,'LineWidth', 3);set(gca,'FontSize',14); % number of cells in the system
%h_legend = legend('1st component','2nd component','Total');
%set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)



%%%%%%%%%%%%%%%%
% solution to the corresponding logistic growth equation

%cA0 = Q_initial/L^2; % initial density

cA0_bar = (Pp - Pd)/(Pp)*S.save_cA0;

f = @(t) cA0_bar * exp(t)/(1+cA0_bar*(exp(t)-1)); % solution to the equation 

for i = 1 : length(S.save_time)
       y_log(i) = f(S.save_time(i)); 
end

%%% plot logistic solution 

hold on
plot(S.save_time,y_log,'LineWidth', 3)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)


%%%%%%%%%%%%%%%%%%%%%%%%%
% adding KSA


%cA0 = Q_initial/L^2; % initial density of cells

dt = 0.01; % time step

N = round(S.save_store_time(length(1:round(length(S.save_store_time)/TimeReductionFactor)))/dt); % number of simulations, chosen to be the same as for Gillespie

% Additional parameters

max_rad = sqrt(2*(L/2)^2); % maximal radius that would cover all the lattice
% Initial conditions 

dr = 0.5; % radially approximated at this distance (after distance 5), i.e 5, 5 + dr, 5 + 2dr

numberDer = (round(max_rad) - 5)/dr ; % this is to cound the number of additional derivatives
     
y = zeros(14 + numberDer,1); % 14 fixed derivatives, the rest varying depending on the size of the grid

y(:,1) = S.save_cA0^2 * ones(14+numberDer,1); % initial density of two-point distribution functions 
y(1,1) = S.save_cA0; % initial condition for 1-point distribution function


t = 0; % start time

dt = 0.01; % time step

% 4th order Runge Kutta

for i = 2 : N
    
    t = t + dt;

    k1 = dynamics (t,y(:,i-1),Pm,Pp,Pd,max_rad,dr);
    k2 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k1,Pm,Pp,Pd,max_rad,dr);
    k3 = dynamics (t+dt/2,y(:,i-1)+dt/2 * k2,Pm,Pp,Pd,max_rad,dr);
    k4 = dynamics (t+dt, y(:,i-1) + dt*k3,Pm,Pp,Pd,max_rad,dr);
    y(:,i) = y(:,i-1) + dt/6 * (k1 + 2*k2 + 2*k3 + k4);  
    
    store_time_Runge(i) = t;
end


rescaled_time = (Pp-Pd)*store_time_Runge; % rescale time to allow parameter comparison

y_rescaled = (Pp-Pd)/Pp * y;


hold on
plot(rescaled_time,y_rescaled(1,:),'LineWidth', 3)
h_legend = legend('Gillespie 2D','Logistic growth','Power 2');
set(h_legend,'FontSize',14)
xlabel('rescaled time','FontSize',14)
ylabel('rescaled density','FontSize',14)
title(['L = ' num2str(L) ', P_m = ' num2str(Pm) ', P_p = ' num2str(Pp) ', P_d = ' num2str(Pd) ' '],'FontSize',14)
set(gca,'linewidth',3)
set(gca,'FontWeight','bold')
set(gca,'FontSize',24)






function deriv = dynamics(t,y,Pm,Pp,Pd,max_rad,dr)
    
    

    numberDer = (round(max_rad) - 5)/dr ; % this is to cound the number of additional derivatives
     
    deriv = zeros(14 + numberDer,1);

    deriv(1) = (Pp*(y(1)-y(2)) - Pd *y(1)); % dynamics of one-point distribution function
    
    deriv(2) = (Pm*0.5* (y(4) + 2 * y(3) - 3* y(2)) - 2 * Pd*y(2) + ...
        Pp * 0.5 * ((y(1)-y(2))*(y(4)+2*y(3))/y(1) + 3 * (y(1)-y(2))^2/(1-y(1)) - 3*y(1)^2*(1-y(1)))); % distance 1
    
    deriv(3) = 0.5 * Pm * (y(2) + y(5) + y(2) + y(5) - 4*y(3)) - ...
        2 * Pd * y(3) + 0.5 * Pp * ((y(2) + y(5) + y(2) + y(5))*(2*y(1) - y(2) - y(3))/(y(1)) + 4 * (y(1)-y(3))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1))); % distance srt(2)
    
    deriv(4) = 0.5 * Pm * (y(7) + y(2) + y(5) + y(5) - 4*y(4)) - ...
        2 * Pd * y(4) + 0.5 * Pp *((y(7) + y(2) + y(5) + y(5))*(2*y(1) - y(2) - y(4))/(y(1)) + 4 * (y(1)-y(4))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1))); % distance 2
    
    deriv(5) = 0.5 * Pm * (y(8) + y(3) + y(6) + y(4) - 4*y(5)) - ...
        2 * Pd * y(5) + 0.5 * Pp * ((y(8) + y(3) + y(6) + y(4))*(2*y(1) - y(2) - y(5))/(y(1)) + 4 * (y(1)-y(5))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1))); % distance sqrt(5)
    
    deriv(6) = 0.5 * Pm * (y(9) + y(5) + y(9) + y(5) - 4*y(6)) - ...
        2 * Pd * y(6) + 0.5 * Pp * ((y(9) + y(5) + y(9) + y(5))*(2*y(1) - y(2) - y(6))/(y(1)) + 4 * (y(1)-y(6))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1)));  % distance sqrt(8)
    
    deriv(7) = 0.5 * Pm * (y(10) + y(4) + y(8) + y(8) - 4*y(7)) - ...
        2 * Pd * y(7) + 0.5 * Pp * ((y(10) + y(4) + y(8) + y(8))*(2*y(1) - y(2) - y(7))/(y(1)) + 4 * (y(1)-y(7))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1)));  % distance 3     
    
    deriv(8) = 0.5 * Pm * (y(11) + y(9) + y(5) + y(7) - 4*y(8)) - ...
        2 * Pd * y(8) + 0.5 * Pp * ((y(11) + y(9) + y(5) + y(7))*(2*y(1) - y(2) - y(8))/(y(1)) + 4 * (y(1)-y(8))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1))); 
    
    deriv(9) = 0.5 * Pm * (y(12) + y(8) + y(6) + y(13) - 4*y(9)) - ...
        2 * Pd * y(9) + 0.5 * Pp * ((y(12) + y(8) + y(6) + y(13))*(2*y(1) - y(2) - y(9))/(y(1)) + 4 * (y(1)-y(9))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1))); 
    
    deriv(10) = 0.5 * Pm * (y(7) + y(14) + y(11) + y(11) - 4*y(10)) - ...
        2 * Pd * y(10) + 0.5 * Pp * ((y(7) + y(14) + y(11) + y(11))*(2*y(1) - y(2) - y(10))/(y(1)) + 4 * (y(1)-y(10))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1))); 
    
    deriv(11) = 0.5 * Pm * (y(13) + y(10) + y(8) + y(15) - 4*y(11)) - ...
        2 * Pd * y(11) + 0.5 * Pp * ((y(13) + y(10) + y(8) + y(15))*(2*y(1) - y(2) - y(11))/(y(1)) + 4 * (y(1)-y(11))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1))); 
        
    deriv(12) = 0.5 * Pm * (y(9) + y(9) + y(14) + y(14) - 4*y(12)) - ...
        2 * Pd * y(12) + 0.5 * Pp * ((y(9) + y(9) + y(14) + y(14))*(2*y(1) - y(2) - y(12))/(y(1)) + 4 * (y(1)-y(12))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1))); 
    
    deriv(13) = 0.5 * Pm * (y(11) + y(14) + y(9) + y(15) - 4*y(13)) - ...
        2 * Pd * y(13) + 0.5 * Pp * ((y(11) + y(14) + y(9) + y(15))*(2*y(1) - y(2) - y(13))/(y(1)) + 4 * (y(1)-y(13))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1)));   % distance sqrt(20)
    
    
    %%%% If dr =0.5 fixed, I will not vary it, then
    
    deriv(14) = 0.5 * Pm * (y(13) + y(12) + y(16) + y(16) - 4*y(14)) - ...
        2 * Pd * y(14) + 0.5 * Pp * ((y(13) + y(12) + y(16) + y(16))*(2*y(1) - y(2) - y(14))/(y(1)) + 4 * (y(1)-y(14))*(y(1)-y(2))/(1-y(1)) - 4 *y(1)^2*(1-y(1))); % distance 5
    
    
    
    
    %%%% be careful because actually y(16) depend on dr    FALSE!!!!!!!!!
%     if sqrt(32) < 5 + dr && sqrt(34) < 5 +dr
%     deriv(14) = 0.5 * Pm * (y(13) + y(12) + y(15) + y(15) - 4*y(14)) - ...
%         2 * Pd * y(14) + 0.5 * Pp * (y(1) - y(14))* (y(1) - y(2))*(y(13)+y(12)+y(15)+ y(15))/(y(1)^2*(1-y(1)));
%     elseif  sqrt(32) > 5 +dr
%     deriv(14) = 0.5 * Pm * (y(13) + y(12) + y(16) + y(16) - 4*y(14)) - ...
%         2 * Pd * y(14) + 0.5 * Pp * (y(1) - y(14))* (y(1) - y(2))*(y(13)+y(12)+y(16)+ y(16))/(y(1)^2*(1-y(1)));
%     else 
%     deriv(14) = 0.5 * Pm * (y(13) + y(12) + y(15) + y(16) - 4*y(14)) - ...
%         2 * Pd * y(14) + 0.5 * Pp * (y(1) - y(14))* (y(1) - y(2))*(y(13)+y(12)+y(15)+ y(16))/(y(1)^2*(1-y(1)));       
%     end       
    % now the distance is more than 5, so all will be approximated
    % uniformly 
    
   
    
    current_radius = 5;
    
    i = 15; % will start considering regular grid from this derivative
    
    % update derivatives until 1 before maximal radius
    while current_radius +dr < max_rad -dr % this is very approximate
        
        current_radius = current_radius + dr; % increase radius by dr
        
        deriv(i) = Pm * (y(i-1)+ y(i+1) - 2 *y(i))/dr^2 + 1/current_radius * (y(i+1) - y(i-1))/(2*dr) ...
            - 2 * Pd * y(i) + Pp * ((y(i+1) + y(i-1) )*(2*y(1) - y(2) - y(i))/(y(1)) + 2 * (y(1)-y(i))*(y(1)-y(2))/(1-y(1)) - 2 *y(1)^2*(1-y(1)));
        
        
        i = i+1;
    end
    
    last = 14 + numberDer;
    %periodic boundary conditions
    deriv(last) =  Pm * (y(last-1)+ y(2) - 2 *y(last))/dr^2 + 1/current_radius * (y(2) - y(last-1))/(2*dr) ...
            - 2 * Pd * y(last) + Pp * ((y(2) + y(last-1) )*(2*y(1) - y(2) - y(last))/(y(1)) + 2 * (y(1)-y(last))*(y(1)-y(2))/(1-y(1)) - 2 *y(1)^2*(1-y(1)));
   
    
end



