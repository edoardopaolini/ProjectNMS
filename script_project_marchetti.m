%% SCRIPT FOR THE PROJECT
%{
This is the script which tries to reproduce the results of the paper and
then tries to extend the vision with other pathways.
%}

%% CODE

% To make molar concentrations equal to abundances, we set the volume as
% the inverse of the Avogadro number
avogadro = 6.02214076e23;
volume = 1/avogadro;

% 1 grammo = 6.022 x 10^23 Da = 1 avogadro Da

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


% initial concentration of Thrombin = 0.75 units/mL
% Thrombin Molecular mass: 37.4 kDa
% Specific Activity: >= 2,000 NIH units/mg
% 1 NIH unit = 0.324 +- 0.073 microgrammi
m_tr1 = 0.324 * 0.75; % conversion unit-microgrammi * units nostre
m_tr2 = 37400/(avogadro*1e6);
% initial amount of fibrinogen
thb = m_tr1/m_tr2;


