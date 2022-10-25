function [correl_R_L_speed, correl_ang_fwd_speed, nb_cycle, ang_speed, fwd_speed] = Neural_oscillator(exhaust_rate, reciprocal_inhibition_strength, time, display)


% Step base model of two simple 'steady state' neurons inhibiting each other.
%The neurons grow an internal self feedback (R_NF and L_NF) that
%lead each of them independantly to stabilise their firing rate (R and L) towards a steady state (default = 0.5) 


%display = 1; % if display = 1; the code will display the graph coded at the
%bottom

%% default parameters (the commented ones can be stated when calling the function)
% time = 1000; % number of steps
% exhaust_rate = 0.01; %rate at which the internal feedback is modulated
% reciprocal_inhibition_strength = 1; % synpatic weight of both reciprocal inhibitory synapses.
steady_state = 0.5; % background firing rate of the neuron on which they would stabilise in absence of external modulation
noise = 0.01; % add a random value (drawn from a normal distribution of std = noise) to the firing rate of each neuron at each time step 

%% declaration of cumulative variables

R=[]; % actual neuron firing rate
R_NF = []; % internal feedback of neuron
L=[]; % actual neuron firing rate
L_NF = []; % internal feedback of neuron

%% initial conditions

R(1) = 0.5;
L(1) = 0.6;
R_NF(1) = 0;
L_NF(1) = 0;

%% run the simulation

for i = 1:time
    
    %calculate novel firing rate as sum of self activation, internal feedback, and simultaneous reciprocal lateral
    %inhibition, + noise
    
    R(i+1) = R(i) - R_NF(i)  - (L(i) * reciprocal_inhibition_strength) + randn(1)*noise;
    L(i+1) = L(i) - L_NF(i)  - (R(i) * reciprocal_inhibition_strength) + randn(1)*noise;
    
  
    %Set hard bondaries to the neurons firing rate
    if R(i+1) < 0; R(i+1) = 0; end
    if R(i+1) > 2; R(i+1) = 2; end
    if L(i+1) < 0; L(i+1) = 0; end
    if L(i+1) > 2; L(i+1) = 2; end

    
    %update the internal feedback value for next step
    R_NF(i+1) = R_NF(i) + (R(i+1) - (R_NF(i) + steady_state))*exhaust_rate;
    L_NF(i+1) = L_NF(i) + (L(i+1) - (L_NF(i) + steady_state))*exhaust_rate;
    
%     R_NF(i+1) = R_NF(i) + (R(i+1) -  steady_state)*exhaust_rate; % Works also with pacemaker neurons
%     L_NF(i+1) = L_NF(i) + (L(i+1) -  steady_state)*exhaust_rate;
    
end

%% Calculate the second order variables resulting from the simulation

    ang_speed = R-L;
    fwd_speed = R+L;
    
    % Check whether L and R are anticorrelated (peak in antiphase)
    r_mat = corrcoef(R,L); % Calculate matrice of correlation
    correl_R_L_speed = r_mat(2,1);

    % Check whether angular and forward speed are correlated (peak in the same time)
    r_mat = corrcoef(abs(ang_speed'),fwd_speed'); % Calculate matrice of correlation
    correl_ang_fwd_speed = r_mat(2,1);

    % Check number of oscillations cycle (proxy for 1/frequency)
    a =  ang_speed (2:end) .* ang_speed(1:end-1);
    f = find(a<0); % get events index when ang_speed crosses 0
    nb_cycle = length(f) /2;


%% Plot the figures

if display == 1

    subplot(4,1,1); 
    plot([R;L]');
    title ('LAL neurons')

    subplot(4,1,2); 
    plot(ang_speed,'g'); hold on,
    plot([0,time+1],[0,0],'k--')
    title('ang speed (L - R)')

    subplot(4,1,3);
    plot(fwd_speed,'k'); hold on,
    title('fwd speed (L + R)')
    
    %check steadiness (if crossess0 events are regular)
    subplot(4,1,4);
    plot(f,'o');
    title('steadiness of turning reversal');
    xlabel ('event number');
    ylabel ('time when event occured')
    
    
end


   