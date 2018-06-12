%% Written by Sean McQuade and Zheming An on July 7th 2017.  
% Model attributed to Benedetto Piccoli and Kwangwon Lee. This code was
% implemented to produce 280 trajectories in the paper titled: 
%"Experimental and Mathematical Analyses Relating Circadian Period 
% and Phase of Entrainment in Neurospora crassa "

%%parameters

global C_Z C_P epsilon tau

% C_Z is the "entrainment strength", C_P is the "alignment strength".
C_Z_list = linspace(0.28,0.28,1);                     % List of values we will use for C_Z
C_P_list = linspace(0.16,0.16,1);                     % List of values we will use for C_P
epsilon = pi/24;                                      % Entrainment window size (radius)
alpha_N = -pi/4;                                      % Initial angle of endogenous clock
alpha_P = pi/6;                                       % Initial angle of light-sensitive protein
tau = 26;                                             % Endogenous period of the biological clock
number_of_trials = 1;           % This value is useful only when we sample other parameters randomly.
T_interval = 0.01;              % Time step
T_initial=T_interval;           % Initial time
T_final =600;                    % Final time
number_of_timesteps = T_final/T_interval;

%% preallocate
circadian_clocks = zeros(number_of_timesteps,3);      % Matrix to store simulation data of three oscillators
phase_of_entrainment = zeros(floor(T_final/24),1);    % Matrix to store phase displacement

for c1 = 1:1
    C_Z = C_Z_list(c1);
for c2 = 1:1
    C_P = C_P_list(c2);
for k=1:number_of_trials
    %The initial condition of three oscillators      
    theta_N_initial = alpha_N;
    theta_P_initial = alpha_P;
    initial_conditions = [0, theta_N_initial,theta_P_initial];
    circadian_clocks(1,:) = initial_conditions;
    
    for i = 1:number_of_timesteps-1
%% The Runge-Kutta method
%k1 = h*f(y) , k2 = h*f(y+1/2*k1), k3 = h*f(y+1/2*k2) , k4 = h*f(y+k3)
%new y = y + 1/6(k1)+1/3(k2)+1/3(k3)+1/6(k4)
%step2: k2 is found using the derivative found using k1 to advance the
%state with a 1/2 step in the k1 direction.  from these "k1_states" we find
%the derivative and call it k2.
%step3: k3 is the derivative of y after you step 1/2 in the k2 direction.
%step4: k4 is the derivative after stepping in the k3 direction
%otherwise, new y = old y + 1/6(k1 +2*(k2) + 2*(k3) + k4) 
        
        k1= Derivative_Runge_Kutta_4th_order_method( circadian_clocks(i,:) );
        k2_clocks = circadian_clocks(i,:) + 0.5*T_interval*k1';

        k2= Derivative_Runge_Kutta_4th_order_method(k2_clocks);
        k3_clocks = circadian_clocks(i,:) + 0.5*T_interval*k2';

        k3= Derivative_Runge_Kutta_4th_order_method(k3_clocks);
        k4_clocks = circadian_clocks(i,:) + T_interval*k3';

        k4= Derivative_Runge_Kutta_4th_order_method(k4_clocks);
        circadian_clocks(i+1,:) =  circadian_clocks(i,:) + T_interval*0.1667*(k1' + 2*k2' + 2*k3' + k4');

% Calulate the sinusoidal values of phase angle for three oscillators
        zeitgeber_signal = sin(circadian_clocks(:,1));
        circadian_signal = sin(circadian_clocks(:,2));
        p_signal = sin(circadian_clocks(:,3));          
    end
%% find the time values corresponding to the peaks of the signal
        [n_max, n_time]= findpeaks(circadian_signal,'MinPeakHeight',0.99);
        [z_max, z_time]= findpeaks(zeitgeber_signal,'MinPeakHeight',0.99);
        days = size(z_time,1);
                
        for j = 1:days
        [phase_abs_value, phase_index] = min(abs(z_time(j) - n_time));
        phase_of_entrainment(j) = (z_time(j) - n_time(phase_index)) .* T_interval;
        end
        
%% Plot the signals for three oscillators and phase of entrainment vs time.      
        figure; %signals representing the three oscillators
        plot(T_initial:T_interval:T_final,zeitgeber_signal,T_initial:T_interval:T_final,circadian_signal,T_initial:T_interval:T_final,p_signal)
        set(gca,'xaxisLocation', 'bottom')
        legend('zeitgeber','endogenous clock','light-sensitive protein')
        title_string = sprintf('C_Z = %0.2f, C_P = %0.2f, ', C_Z, C_P);
        title(title_string);
        xlabel('time')
        ylabel('signal')
        set(gca,'fontsize',18);
        
        figure; %plots the phase of entrainment changing over time
        plot(z_time*T_interval,phase_of_entrainment,'LineWidth',2)
        title2 = sprintf('Phase of entrainment, Period = %0.2f', tau );
        title(title2);
        xlabel('time')
        ylabel('phase of entrainment')
        set(gca,'xaxisLocation', 'bottom')
        set(gca,'fontsize',18);

        hold on

end
end
end
