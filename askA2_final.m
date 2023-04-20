%% Define the chemical reaction network
% Here, we consider the Michaelis-Menten reaction with enzyme kinetics:
% E + S <-> ES -> E + P
% with parameters k1, k2, and k3
clear; close all; clc

k1=2; k1r=1; k2=1.5;
cs0=8; ce0=4; ces0=0; cp0=0;

constants=[k1 k1r k2];
c0=[ce0 cs0 ces0 cp0];

% Define the simulation parameters
tmin=0; tmax=6;
steps=1000;
tspan=linspace(tmin,tmax,steps);

% Solve ode for concentrations: [E],[S],[ES],[P]
[ti1,ci]=ode89(@(t,c)odeA2(t,c,constants),tspan,c0);

% Plotting the results
figure(1)
sgtitle('ode89')
subplot(2,2,1);
plot(ti1,ci(:,1),'r','LineWidth',1.5)
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[E]')
grid on; box on; axis tight
subplot(2,2,2);
plot(ti1,ci(:,2),'g','LineWidth',1.5)
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[S]')
grid on; box on; axis tight
subplot(2,2,3);
plot(ti1,ci(:,3),'b','LineWidth',1.5)
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[ES]')
grid on; box on; axis tight
subplot(2,2,4);
plot(ti1,ci(:,4),'y','LineWidth',1.5)
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[P]')
grid on; box on; axis tight




%% Gillespies algorithm for Michaelis-Menten kinetics

selection = 2;

switch selection
    case 1
        VL = 1e-17;
    case 2 
        VL = 1e-16;
    case 3
        VL = 1e-15;
end

ce0=ce0*1e-6; cs0=cs0*1e-6; ces0=ces0*1e-6; cp0=ce0*1e-6;
NA=6.022e23;

constants=[k1*1e6/(NA*VL) k1r k2];

% Define the simulation parameters
tmin=0; 
maxiter=100000; % maximum iterations in each simulation trial

ntrials=30; % number of simulation trials

% Run the simulation trials
for i=1:ntrials

     % Define the initial conditions for the current trial
    ti2(1)=tmin;
    X1(1)=round(ce0*(NA*VL)); % #molecules E
    X2(1)=round(cs0*(NA*VL)); % #molecules S
    X3(1)=round(ces0*(NA*VL)); % #molecules ES
    X4(1)=round(cp0*(NA*VL)); % #molecules P
    
    xout{i}(1)=X1(1);
    xout{ntrials+i}(1)=X2(1);
    xout{2*ntrials+i}(1)=X3(1);
    xout{3*ntrials+i}(1)=X4(1);
    tout{i}(1)=ti2(1);

    % Main loop
    iter=1;
    while true
        
        % Compute the reaction rates
        a=[X1(end)*X2(end)*constants(1) X3(end)*constants(2) X3(end)*constants(3)];
        a0=sum(a);

        % Generate the two random numbers needed for the algorithm
        r=rand(1,2);
        
        % Compute the time increment
        dt=(1/a0)*log(1/r(1));

        % Update the time
        ti2(end+1)=ti2(end)+dt;
        
        % Determine which reaction occurs
        if r(2)*a0<a(1)
            % Reaction 1: E + S -> ES
            m=1;
        elseif r(2)*a0 < (a(1)+a(2))
            % Reaction 2: ES -> E + S
            m=2;
        else
            % Reaction 3: ES -> E + P
            m=3;
        end
        
        % Update the #molecules of reactants and products
        switch m
        case 1
            X1(end+1) = X1(end)-1;
            X2(end+1) = X2(end)-1;
            X3(end+1) = X3(end)+1;
            X4(end+1) = X4(end);
        case 2 
            X1(end+1) = X1(end)+1;
            X2(end+1) = X2(end)+1;
            X3(end+1) = X3(end)-1;
            X4(end+1) = X4(end);
        case 3
            X1(end+1) = X1(end)+1;
            X2(end+1) = X2(end);
            X3(end+1) = X3(end)-1;
            X4(end+1) = X4(end)+1;
        end
        
        % Store the results
        flag = isinf(ti2(end));
        if flag==0
            tout{i}(iter)=ti2(end);
            xout{i}(iter)=X1(iter);
            xout{ntrials+i}(iter)=X2(iter);
            xout{2*ntrials+i}(iter)=X3(iter);
            xout{3*ntrials+i}(iter)=X4(iter);
        end

        % Termination criteria
        if ti2(iter)>10
            break;end
        % Update the main counter
        iter=iter+1;

        if iter>maxiter
            break;end
    end %end inner loop
    
    figure(2)
    subplot(2,2,1);
    plot(ti2,X1/(NA*VL)*1e6,'r','LineWidth',0.01)
    hold on
    subplot(2,2,2);
    plot(ti2,X2/(NA*VL)*1e6,'g','LineWidth',0.01)
    hold on
    subplot(2,2,3);
    plot(ti2,X3/(NA*VL)*1e6,'b','LineWidth',0.01)
    hold on
    subplot(2,2,4);
    plot(ti2,X4/(NA*VL)*1e6,'y','LineWidth',0.01)
    hold on
    

    clear ti2 X1 X2 X3 X4 iter a

end %end outer loop

% Loop over the simulation runs
for j=1:ntrials
    % Interpolate the data to a common time vector
    tmax = max(cellfun(@max,tout));
    tnew = linspace(tmin, tmax, steps);
    Enew{j} = interp1(tout{j},xout{j},tnew);
    Snew{j} = interp1(tout{j},xout{ntrials+j}, tnew);
    ESnew{j} = interp1(tout{j},xout{2*ntrials+j}, tnew);
    Pnew{j} = interp1(tout{j},xout{3*ntrials+j}, tnew);
    
    % Compute the running averages
    if j == 1
        tavg = tnew;
        Eavg = Enew{1};
        Savg = Snew{1};
        ESavg = ESnew{1};
        Pavg = Pnew{1};
    else
        tavg = (tavg*(j-1) + tnew)/j;
        Eavg = (Eavg*(j-1) + Enew{j})/j;
        Savg = (Savg*(j-1) + Snew{j})/j;
        ESavg = (ESavg*(j-1) + ESnew{j})/j;
        Pavg = (Pavg*(j-1) + Pnew{j})/j;
    end
end

% Converting #molecules to concentrations
Eavg = Eavg/(NA*VL)*1e6;
Savg = Savg/(NA*VL)*1e6;
ESavg = ESavg/(NA*VL)*1e6;
Pavg = Pavg/(NA*VL)*1e6;

% Plot the average curve
figure(2)
sgtitle('gillespie');
subplot(2,2,1);
plot(tavg,Eavg,'--k','LineWidth',1.2)
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[E]')
grid on; box on; axis tight
subplot(2,2,2);
plot(tavg,Savg,'--k','LineWidth',1.2)
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[S]')
grid on; box on; axis tight
subplot(2,2,3);
plot(tavg,ESavg,'--k','LineWidth',1.2)
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[ES]')
grid on; box on; axis tight
subplot(2,2,4);
plot(tavg,Pavg,'--k','LineWidth',1.2)
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[P]')
grid on; box on; axis tight

% Comparison of gillespie's algorithm with ode solver
figure(3)
sgtitle(" gillespie's algorithm vs ode89 ")
subplot(2,2,1);
plot(ti1,ci(:,1),'r','LineWidth',1.5)
hold on
plot(tavg,Eavg,'--k','LineWidth',1.2)
hold off
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[E]'); legend('ode89','gillespie')
grid on; box on; axis tight
subplot(2,2,2);
plot(ti1,ci(:,2),'g','LineWidth',1.5)
hold on
plot(tavg,Savg,'--k','LineWidth',1.2)
hold off
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[S]'); legend('ode89','gillespie')
grid on; box on; axis tight
subplot(2,2,3);
plot(ti1,ci(:,3),'b','LineWidth',1.5)
hold on
plot(tavg,ESavg,'--k','LineWidth',1.2)
hold off
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[ES]'); legend('ode89','gillespie')
grid on; box on; axis tight
subplot(2,2,4);
plot(ti1,ci(:,4),'y','LineWidth',1.5)
hold on
plot(tavg,Pavg,'--k','LineWidth',1.2)
hold off
xlabel('t (s)'); ylabel('concenctration (μM)'); title('[P]'); legend('ode89','gillespie')
grid on; box on; axis tight



% Relative error
r_E = abs(Eavg'-ci(:,1))./ci(:,1)*100;
r_S = abs(Savg'-ci(:,2))./ci(:,2)*100;
r_ES = abs(ESavg'-ci(:,3))./ci(:,3)*100;
r_P = abs(Pavg'-ci(:,4))./ci(:,4)*100;

figure(4)
subplot(2,2,1);
plot(tavg,r_E)
xlabel('t (s)'); ylabel('Percentage error'); title('[E]')
grid on; box on; axis tight
subplot(2,2,2);
plot(tavg,r_S)
xlabel('t (s)'); ylabel('Percentage error'); title('[S]')
grid on; box on; axis tight
subplot(2,2,3);
plot(tavg,r_ES)
xlabel('t (s)'); ylabel('Percentage error'); title('[ES]')
grid on; box on; axis tight
subplot(2,2,4);
plot(tavg,r_P)
xlabel('t (s)'); ylabel('Percentage error'); title('[P]')
grid on; box on; axis tight


