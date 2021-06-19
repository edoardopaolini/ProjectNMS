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
%fbni = 2.5/mass_fib;
fbni = 1000;

% THROMBIN
% initial concentration of Thrombin = 0.75 units/mL
% Thrombin Molecular mass: 37.4 kDa
% Specific Activity: >= 2,000 NIH units/mg
% 1 NIH unit = 0.324 +- 0.073 microgrammi
m_tr1 = 0.324 * 0.75; % conversion unit-microgrammi * units nostre
m_tr2 = 37400/(avogadro*1e6);
% initial amount of fibrinogen
%thb = m_tr1/m_tr2;
thb = 300;

% ACTIVE FIBRINOGEN
fbna = 0;

% FIBRIN MATRIX
fm = 50;

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
r(1) = 0.0284;
r(2) = 0;
r(3) = 2064;
r(4) = 0.484;
r(5) = 290;
r(6) = 300;
r(7) = 0.286;
r(8) = 2.64;
r(9) = 22.7;
r(10) = 15.5;
r(11) = 366;

%r = r./(max(r));

% Initial Amount
Xinit(1) = fbna;
Xinit(2) = fm;
Xinit(3) = thb;
Xinit(4) = fbni;
Xinit(5) = c0;
Xinit(6) = c1;
Xinit(7) = c2;

mi = max(Xinit);
%Xinit = Xinit
choice = input('Vuoi i dati normalizzati? [y/n]: ');
if choice=='y'
    Xinit = Xinit./mi;
end

% ODEs
% fbna derivatives
%{
Z(1) = -r(4)*Xinit(1)*Xinit(3) + r(5)*Xinit(6) - r(8)*Xinit(1)*Xinit(6) + r(9)*Xinit(7) + r(3)*Xinit(5);
Z(2) = r(6)*Xinit(6) - r(7)*Xinit(3)*Xinit(2) + r(10)*Xinit(7) - r(11)*Xinit(6)*Xinit(2);
Z(3) = -r(4)*Xinit(1)*Xinit(3) + r(5)*Xinit(6) + r(6)*Xinit(6) - r(7)*Xinit(3)*Xinit(2) ...
    - r(1)*Xinit(4)*Xinit(3) + r(2)*Xinit(5) + r(3)*Xinit(5);
Z(4) = -r(1)*Xinit(4)*Xinit(3) + r(2)*Xinit(5);
Z(5) = r(1)*Xinit(4)*Xinit(3) -  r(2)*Xinit(5) - r(3)*Xinit(5);
Z(6) = r(4)*Xinit(1)*Xinit(3) - r(5)*Xinit(6) - r(6)*Xinit(6) + r(7)*Xinit(3)*Xinit(2) ...
    + r(9)*Xinit(7) - r(8)*Xinit(1)*Xinit(6) + r(10)*Xinit(7) - r(11)*Xinit(6)*Xinit(2);
Z(7) = r(8)*Xinit(1)*Xinit(6) - r(9)*Xinit(7) + r(11)*Xinit(6)*Xinit(2) - r(10)*Xinit(7);
%}


%% ODEs Simulation with ODE15S

%integration bounds
tin=0;
tfin=255;

%integrate the model
[T,Y] = ode15s(@(t,X) enzReact(t,X,r), tin:1:tfin, Xinit);

% To get the same colors to the data and simulations

%set colors as Matlab default

colorMap = lines(7);

if choice=='y'
    Y = Y.*mi;
end

figure(1)
for i=1:7
    plot(T,Y(:,i), 'LineWidth', 1.5, 'color', colorMap(i,:));
    hold on
end
legend({'FBNa', 'FM','THB' , 'FBNi' , 'C0', 'C1', 'C2'})
hold off

%% OPTIMIZATION

Tolerance = 1e-15;
if choice=='y'
    Y = Y./mi;
end

% Initial Rates
r2 = zeros(11,1);
r2(1) = 10000;
r2(2) = 100;
r2(3) = 100;
r2(4) = 100;
r2(5) = 100;
r2(6) = 100;
r2(7) = 1;
r2(8) = 1;
r2(9) = 100;
r2(10) = 1;
r2(11) = 1;

ms              = MultiStart('UseParallel','always','Display','iter','StartPointsToRun','bounds','TolFun',Tolerance,'TolX',Tolerance);

nOptRuns        = 12;
x0              = r2;
lb              = zeros(size(x0));
ub              = 10*ones(size(x0));
problem         = createOptimProblem('lsqnonlin','objective',@(par) myobjfun(par, 0:1:255, Y, 0, 255 , Xinit), ...
                'x0' ,x0 ,'lb',lb,'ub',ub);

% problem         = createOptimProblem('lsqnonlin','objective',@(par) myobjfun(par, TimeExpDataEnzReact, ExpDataEnzReact, 0, 300 , Xiniz), ...
%                 'x0' ,s2 ,'lb',lb,'ub',ub, 'options', opts); 
% more precise but slower

EstParametersMS = run(ms,problem,nOptRuns);

%% Determinismo con nuovi parametri

%integration bounds
tin=0;
tfin=255;
r = EstParametersMS;

%integrate the model
[T,Y] = ode15s(@(t,X) enzReact(t,X,r), tin:1:tfin, Xinit);

% To get the same colors to the data and simulations

%set colors as Matlab default

colorMap = lines(7);

if choice=='y'
    Y = Y.*mi;
end

figure(2)
for i=1:7
    plot(T,Y(:,i), 'LineWidth', 1.5, 'color', colorMap(i,:));
    hold on
end
legend({'FBNa', 'FM','THB' , 'FBNi' , 'C0', 'C1', 'C2'})
hold off


%% Stocastici

% fi(1), fa(2), t(3), c0(4), c1(5), c2(6), fm(7)
% fbna(1), fm(2), thb(3), fbni(4), c0(5), c1(6), c2(7)

vMinus = [1 0 1 0 0 0 0; 0 0 0 1 0 0 0; 0 0 0 1 0 0 0; 0 1 1 0 0 0 0; ...
    0 0 0 0 1 0 0; 0 0 0 0 1 0 0; 0 0 1 0 0 0 1; 0 1 0 0 1 0 0; ...
    0 0 0 0 0 1 0; 0 0 0 0 0 1 0; 0 0 0 0 1 0 1];

vPlus = [0 0 0 1 0 0 0; 1 0 1 0 0 0 0; 0 1 1 0 0 0 0; 0 0 0 0 1 0 0; ...
    0 1 1 0 0 0 0; 0 0 1 0 0 0 1; 0 0 0 0 1 0 0; 0 0 0 0 0 1 0; ...
    0 1 0 0 1 0 0; 0 0 0 0 1 0 1; 0 0 0 0 0 1 0];

c = r;

initialState = zeros(7,1);
initialState(1,1) = Xinit(4,1);
initialState(2,1) = Xinit(1,1);
initialState(3,1) = Xinit(3,1);
initialState(4,1) = Xinit(5,1);
initialState(5,1) = Xinit(6,1);
initialState(6,1) = Xinit(7,1);
initialState(7,1) = Xinit(2,1);
initialState = initialState';

delta = 0.15;

tMax = 10;

dT = 0.001;

[T,Y] = simRSSA_disc(vMinus,vPlus,c,initialState,delta,tMax,dT);

colorMap = lines(7);

figure(1)
for i=1:7
    plot(T,Y(:,i), 'LineWidth', 1.5, 'color', colorMap(i,:));
    hold on
end
% fi(1), fa(2), t(3), c0(4), c1(5), c2(6), fm(7)
legend({'FBNi', 'FBNa', 'THB', 'C0', 'C1', 'C2', 'FM'})
hold off





