function [ derivative ] = Derivative_Runge_Kutta_4th_order_method(y)

global C_Z C_P epsilon tau 

%rename the variables from the single input vector.
theta_Z = y(1);      % theta_Z is the phase angle of zeitgeber signal
theta_N = y(2);      % theta_N is the phase angle of endogenous clock signal
theta_P = y(3);      % theta_P is the phase angle of light-sensitive protein signal

%% Using modular function to normalize phase angles of three oscillators.
% We use the normalized values to define the correct phase angles. 
Normalized_theta_Z = mod(theta_Z,2*pi);
Normalized_theta_P = mod(theta_P,2*pi);
Normalized_theta_N = mod(theta_N,2*pi);

%% Function to dsecribe entrainment happens during the dusk (phase angle equal to pi), and C_Z is the "entrainment strength"
if abs(Normalized_theta_Z - pi) < epsilon       % epsilon Entrainment window size (radius)
    Psi_1 = C_Z;
else
    Psi_1 = 0;
end

%% Functions to describe the direction of the angle alignment. 
 % vector_1 defines how the phase angle entrainment happens between light-sensitive protein and zeitgeber signal.  The term is a 
if abs(Normalized_theta_Z - Normalized_theta_P) < pi     
        vector_1 = Normalized_theta_Z-Normalized_theta_P;   
else
    if Normalized_theta_Z - Normalized_theta_P > 0
        vector_1 = Normalized_theta_Z - Normalized_theta_P - 2*pi;
    else
        vector_1 = Normalized_theta_Z - Normalized_theta_P + 2*pi;
    end
end    
 % vector_2 defines how the phase angle alignment happens between endogenous clock and light-sensitive protein.
if abs(Normalized_theta_P - Normalized_theta_N) < pi     
        vector_2 = Normalized_theta_P-Normalized_theta_N;   
else
    if Normalized_theta_P - Normalized_theta_N > 0
        vector_2 = Normalized_theta_P - Normalized_theta_N - 2*pi;
    else
        vector_2 = Normalized_theta_P - Normalized_theta_N + 2*pi;
    end
end    


dtheta_Z = 2*pi/24;                              % The derivative of the phase angle of Zeitgeber signal
dtheta_N = (2*pi)/tau +C_P*(vector_2 + 2*pi/3);  % The derivative of the phase angle of endogenous clock signal, C_P is the "alignment strength"
dtheta_P = (2*pi)/tau+ Psi_1*vector_1;           % The derivative of the phase angle of light-sensitive protein signal

derivative = [dtheta_Z, dtheta_N, dtheta_P]';
end

