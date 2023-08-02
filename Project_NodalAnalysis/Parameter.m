%% Well properties
tub_length = 25;
KOP = 5000;
 %ft
L_h = 3000; %ft
d_tub = 2.441; %inch
rel_rough = 0.0006; %[]
P1 = 1500;%:-100:1000; %psia
T_th = 150; %Farenheit
T1 = Tfah_Rank(T_th); %Rankine


%% Reservoir Properties
A = 350; %acres
r_e = sqrt((A*43560)/pi);
r_w = 2.441/12;

poro = 0.30;
Sw = 0.05;
h = 30; %ft
k_h = 2; %mD
k_v = 0.4; %mD
P_avg_res = 3600;%:-100:2000; %psia
T_res = 200; %Farenheit
T2 = Tfah_Rank(T_res); %Rankine

%% Fluid Properties
Y_g = 0.71; %[]
m_g1 = 0.02; %cp
M_air = 28.97; %lbmol-ft^3
M_g = Y_g*M_air; %lbmol-ft^3 
R = 10.73; %[10.73 (psi·ft3)/(lbmmol·°R)]
%% discretizing length
%i = (0:250:TVD);
% This log avearge is calculated by hand for 250ft of tubing length, this is the Tempt @125ft
L = tub_length;

%% first guessed q


                            
%% z constants
        A1 = 0.3265;
        A2 = -1.0700;
        A3 = -0.5339;
        A4 = 0.01569;
        A5 = -0.05165;
        A6 = 0.5475;
        A7 = -0.7361;
        A8 = 0.1844;
        A9 = 0.1056;
        A10 = 0.6134;
        A11 = 0.7210;



