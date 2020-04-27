clear; close all; clc

% Load raw data from file "data_buoy44029.mat" 
% http://www.ndbc.noaa.gov/station_page.php?station=44029

% MAT File containts three matrices of data, subdivided by year: 2008, 2009 , 2010
% Each columns is:
% 1     2   3   4   5   6       7       8
% #YY	MM	DD	hh	mm	WDIR	WSPD	GST

load('data_buoy44029.mat')

U_data0=[data_Y2008(:,7);data_Y2009(:,7);data_Y2010(:,7)];

%% Data conversion due to averaging time (8 min to 10 min)
U_data1=U_data0*1.05/1.07;

%% Estimate wind speed at the elevation of the rotor (nacelle) by log-law
z0=0.0002;    % Roughness length in meters - Open sea 
z_rotor=60;   % Elevation of the rotor from sea surface (meters)
z_sensor=4;   % Elevation of the anemometer installed on the buoy (meters)

U_data=U_data1*log(z_rotor/z0)/log(z_sensor/z0);

%% Plotting the histogram of raw data at height of the rotor
hist(U_data,20)

% title, labels, etc.
Ftsize=16;
title('Q2a: Raw data between 2008 and 2010','Fontsize',Ftsize,'FontWeight','bold','Fontname',...
     'Times New Roman')
xlabel('{\itU_i} [m/s]','Fontsize',Ftsize,'FontWeight','bold','Fontname','Times New Roman')
ylabel('Occurrence','Fontsize',Ftsize,'FontWeight','bold','Fontname','Times New Roman')
axset=gca;
set(axset,'Fontname','Times New Roman','Fontsize',Ftsize-2)   
box on


%% DETERMINE "SHAPE" AND "DISPERSION" OF THE WEIBULL DISTRIBUTION REGRESSION METHOD on "Weibull chart"

% Plotting position, Iv, Weibull
U_data_sorted=sort(U_data);     % Data sorting
N=length(U_data);
Iv=(1:1:N)/(N+1);               % Plotting position, Weibull

% Reduced variates, Xwb and Ywb
Ywb=log(-log(1-Iv));
Xwb=log(U_data_sorted);

OK_indices_w=find(abs(Xwb)~=Inf);      % Find finite values log(u) NOT Inf

pwb=polyfit(Xwb(OK_indices_w),Ywb(OK_indices_w)',1);  % Linear regression
Uwb_regr=pwb(1)*Xwb(OK_indices_w)+pwb(2);                % Estimating linear regression line

disp('Q2b: Solution for Weibull shape (k) and scale (alpha), using regression analysis ("Weibull chart")')
k_chart = pwb(1)
alpha_chart=exp(-pwb(2)/pwb(1))

%% Plot CDF of wind speed data

figure
plot(Xwb(OK_indices_w),Ywb(OK_indices_w),'k.',Xwb(OK_indices_w),Uwb_regr,'g--')
title('Q2c: Weibull Fit by Regression on ''Weibull Chart''','Fontsize',Ftsize,'FontWeight','bold','Fontname',...
     'Times New Roman')
xlabel('\itX_w\rm = log(\itu\rm),  \itu\rm: wind speed','Fontsize',Ftsize,'FontWeight',...
       'normal','Fontname','Times New Roman')
ylabel('\itY_w\rm = log[ -log( 1-F\it_u\rm )]','Fontsize',Ftsize,'FontWeight',...
       'normal','Fontname','Times New Roman')
legend('Data (Buoy 44028)','Regression line','location','southeast')
axset=gca;
set(axset,'Fontname','Times New Roman','Fontsize',Ftsize-2)   
box on

%% Calculate Efficiency of wind turbine

% Output Power curve (Pw) vs. wind speed (U_w)
% Pw: output power curve (kiloWatts)
% U_w: wind speed (m/s)

% Interpolated data (from graph)
U_w=[
0
4.00
6.50
7.50
8.85
9.84
10.49
11.00
11.50
12.50
13.75
15.30
16.50
17.50
18.50
19.00
20.00
21.50
22.50
22.50
24.00];

Pw=[
0
25
25
83
263
558
869
1207
1465
1814
2167
2536
2731
2800
2752
2658
2520
30
0
0
0];
    
figure
plot(U_w,Pw,'k*-')
% Tites, labels, etc
title('Q2d: Power Output Curve','Fontsize',Ftsize,'FontWeight',...
       'bold','Fontname','Times New Roman')
xlabel('Wind speed, \itU_w \rm\bf[m/s]','Fontsize',Ftsize,'FontWeight',...
       'bold','Fontname','Times New Roman')
ylabel('Output Power, \itP_w \rm\bf[KWatt]','Fontsize',Ftsize,'FontWeight','bold','Fontname','Times New Roman')
axset=gca;
set(axset,'Fontname','Times New Roman','Fontsize',Ftsize-2)   
box on

Weib_pdf_w=wblpdf(U_w,alpha_chart,k_chart);

figure
plot(U_w,Weib_pdf_w,'k*-')
% Tites, labels, etc
title('Q2d (continued): Weibull PDF vs. \itU_w','Fontsize',Ftsize,'FontWeight','bold','Fontname','Times New Roman')
xlabel('Wind speed, \itU_w \rm\bf[m/s]','Fontsize',Ftsize,'FontWeight','bold','Fontname','Times New Roman')
ylabel('PDF','Fontsize',Ftsize,'FontWeight','bold','Fontname','Times New Roman')
axset=gca;
set(axset,'Fontname','Times New Roman','Fontsize',Ftsize-2)   
box on

% Integrating F(x)dx by means of trapeziodal rule
X=U_w;
F_X=Pw.*Weib_pdf_w;
SumP=0;
for ii=1:length(X)-1
    SumP=SumP+(F_X(ii+1)+F_X(ii))/2*(X(ii+1)-X(ii));
end

disp(' ')
disp('Q2d: Efficiency / productivity of wind turbine (kiloWatts)')
Pw_bar=SumP
