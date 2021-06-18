%% SCRIPT FOR THE PROJECT
%{
This is the script which tries to reproduce the results of the paper and
then tries to extend the vision with other pathways.
%}

%% INITIAL VALUES

% To make molar concentrations equal to abundances, we set the volume as
% the inverse of the Avogadro number
avogadro = 6.02214076e23;
volume = 1/avogadro;

% 1 grammo = 6.022 x 10^23 Da = 1 avogadro Da

% INACTIVE FIBRINOGEN
% initial concentration of Fibrinogen = 2.5 mg/mL
% Fibrinogen molecular mass: ~340 kDa = 340000 Da (NB diversa dala molar
% mass)
% Quindi sapendo a quanto corrisponde il rapporto Dalton e grammi e sapendo
% che una molecola di fibrinogeno pesa circa 340 kDa, possiamo trovare
% quante molecole di fibrinogeno ci sono inizialmente in un mL di soluzione
% calcoliamo la massa di una molecola di fibrinogeno e facciamo 2.5/ans
mass_fib = 340000/(avogadro*1e3);
% initial amount of fibrinogen
fbni = 2.5/mass_fib;
fbni = 1000;

% THROMBIN
% initial concentration of Thrombin = 0.75 units/mL
% Thrombin Molecular mass: 37.4 kDa
% Specific Activity: >= 2,000 NIH units/mg
% 1 NIH unit = 0.324 +- 0.073 microgrammi
m_tr1 = 0.324 * 0.75; % conversion unit-microgrammi * units nostre
m_tr2 = 37400/(avogadro*1e6);
% initial amount of fibrinogen
thb = m_tr1/m_tr2;
thb = 800;

% ACTIVE FIBRINOGEN
fbna = 0;

% FIBRIN MATRIX
fm = 100;

% COMPLEX InatcivateFIBRINOGEN-THROMBIN
c0 = 0;

% COMPLEX AcivateFIBRINOGEN-THROMBIN
c1 = 0;

% COMPLEX AcivateFIBRINOGEN-COMPLEX C1
c2 = 0;


%% ODE SYSTEM

% create the vectors of ODEs
% fbna(1), fm(2), thb(3), fbni(4), c0(5), c1(6), c2(7)
% Z = zeros(7,1);
% create the vectors of concentrations
% fbna(1), fm(2), thb(3), fbni(4), c0(5), c1(6), c2(7)
Xinit = zeros(7,1);
% create vectors of rates
% k+(1), k-(2), k(3), k1+(4), k1-(5), k2+(6), 
% k2-(7), k3+(8), k3-(9), k4+(10), k4-(11)
r = zeros(11,1);

% Initial Rates
r(1) = 0.031;
r(2) = 0;
r(3) = 1931;
r(4) = 0.49;
r(5) = 250;
r(6) = 298;
r(7) = 0.26;
r(8) = 2.6;
r(9) = 26.9;
r(10) = 15.5;
r(11) = 375;

% Initial Amount
Xinit(1) = fbna;
Xinit(2) = fm;
Xinit(3) = thb;
Xinit(4) = fbni;
Xinit(5) = c0;
Xinit(6) = c1;
Xinit(7) = c2;

% ODEs
% fbna derivatives
%{
Z(1) = -r(4)*Xinit(1)*Xinit(3) + r(5)*Xinit(6) - r(8)*Xinit(4)*Xinit(6) + r(9)*Xinit(7) + r(3)*Xinit(5);
Z(2) = r(6)*Xinit(6) - r(7)*Xinit(3)*Xinit(2) + r(10)*Xinit(7) - r(11)*Xinit(6)*Xinit(2);
Z(3) = -r(4)*Xinit(1)*Xinit(3) + r(5)*Xinit(6) + r(6)*Xinit(6) - r(7)*Xinit(3)*Xinit(2) ...
    - r(1)*Xinit(4)*Xinit(3) + r(2)*Xinit(5) + r(3)*Xinit(5);
Z(4) = -r(1)*Xinit(4) + r(2)*Xinit(5);
Z(5) = r(1)*Xinit(4)*Xinit(3) -  r(2)*Xinit(5) - r(3)*Xinit(5);
Z(6) = r(4)*Xinit(1)*Xinit(3) - r(5)*Xinit(6) - r(6)*Xinit(6) + r(7)*Xinit(3)*Xinit(2) ...
    + r(9)*Xinit(7) - r(8)*Xinit(1)*Xinit(6) + r(10)*Xinit(7) - r(11)*Xinit(6)*Xinit(2);
Z(7) = r(8)*Xinit(1)*Xinit(6) - r(9)*Xinit(7) + r(11)*Xinit(6)*Xinit(2) - r(10)*Xinit(7);
%}


%% ODEs Simulation with ODE15S

%integration bounds
tin=0;
tfin=10;

%integrate the model
[T,Y] = ode15s(@(t,X) enzReact(t,X,r), tin:0.05:tfin, Xinit);


% plot the simulation results
figure(1)
plot(T,Y, 'LineWidth', 1.5)
% fbna(1), fm(2), thb(3), fbni(4), c0(5), c1(6), c2(7)
legend({'FBNa', 'FM','THB' , 'FBNi' , 'C0', 'C1', 'C2'})

% To get the same colors to the data and simulations

%set colors as Matlab default

colorMap = lines(6);

figure(1)
plot(T,Y(:,1), 'LineWidth', 1.5, 'color', colorMap(1,:))
hold on 
for i=2:5
    plot(T,Y(:,i), 'LineWidth', 1.5, 'color', colorMap(i,:))
end
hold off





