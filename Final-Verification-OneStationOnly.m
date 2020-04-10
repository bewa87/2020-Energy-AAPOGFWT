# Hard coded - for examplification and test purposes on one station only

###
### Step 1: Data Processing For Weibull Distribution
###

pkg load statistics;

fid = fopen('ReadData-03032-ListSylt.txt', 'r');
C   = textscan(fid,'%f32 %f32 %f32 %f32 %f32 %s', 'delimiter', ',');
fclose(fid);

vel_proc       = cell2mat(C(4));
final_vel_proc = sort(vel_proc);
ind_proc       = find(final_vel_proc>0.0);
final_vel_proc = final_vel_proc(ind_proc);

###
### Step 2: Power Curve Modeling
###

##
##  Step 2.1: Preparing Regression
##

log_vel   = [0:0.5:40];
powww01   = zeros(1,6);
pow01     = [powww01 26 73 133 207];
pow02     = [302 416 554 717 907 1126 1375 1652 1985 2282];
powww02   = 3075*ones(1,25);
pow03     = [2585 2821 2997 3050 3067 3074 powww02];
pow04     = zeros(1,30);
log_power = [pow01 pow02 pow03 pow04];

##
##  Step 2.2: Interval for Regression
##

vellog = log_vel(6:29);
powlog = log_power(6:29);

##
##  Step 2.3: Cubic Regression
##

# Case: SAD - Cubic Approximation

x_cubSAD   = fminsearch (@(x) (sum( abs( x(1)*vellog.^3+x(2)*vellog.^2+x(3)*vellog+x(4) - powlog) )),[0;0;0;0]);
y_cubSAD   = x_cubSAD(1)*vellog.^3+x_cubSAD(2)*vellog.^2+x_cubSAD(3)*vellog+x_cubSAD(4);

# Case: SSD - Cubic Approximation

x_cubSSD   = fminsearch (@(x) (sum( ( x(1)*vellog.^3+x(2)*vellog.^2+x(3)*vellog+x(4) - powlog ).^2 )), [0;0;0;0]);
y_cubSSD   = x_cubSSD(1)*vellog.^3+x_cubSSD(2)*vellog.^2+x_cubSSD(3)*vellog+x_cubSSD(4);

# Cubic Comparison Plot

figure(1);
plot(vellog, powlog, "color", "black");
lt01 = title("Wind Speed v vs. Power Output P", "fontsize", 14);
xt01 = xlabel("v", "fontsize", 14);
yt01 = ylabel("P", "fontsize", 14);
hold on
plot(vellog, y_cubSAD, "color", "green");
hold on
plot(vellog, y_cubSSD, "color", "red");

##
##  Step 2.4: Logistic Regression
##

# Case: SAD - Logistic Function

x_logSAD   = fminsearch (@(x) (sum( abs( (x(1))./(x(2)+x(3).*exp(x(4).*vellog+x(5))) - powlog) )),[1;1;1;1;1]);
y_logSAD   = (x_logSAD(1))./(x_logSAD(2)+x_logSAD(3).*exp(x_logSAD(4).*vellog+x_logSAD(5)));

# Case: SSD - Logistic Function

x_logSSD   = fminsearch (@(x) (sum( ( x(1)./(x(2)+x(3).*exp(x(4).*vellog+x(5))) - powlog ).^2 )), [1;1;1;1;1]);
y_logSSD   = (x_logSSD(1))./(x_logSSD(2)+x_logSSD(3).*exp(x_logSSD(4).*vellog+x_logSSD(5)));

# Logistic Comparison Plot

figure(2);
plot(vellog, powlog, "color", "black");
lt02 = title("Wind Speed v vs. Power Output P", "fontsize", 14);
xt02 = xlabel("v", "fontsize", 14);
yt02 = ylabel("P", "fontsize", 14);
hold on
plot(vellog, y_logSAD, "color", "green");
hold on
plot(vellog, y_logSSD, "color", "red");

###
### Step 3: Weibull Distribution (only taken from R - before: Nelder-Mead for
###         Weibull Distribution)
###

##
##  Step 3.1: Parameter Setup 
##

k_R = 2.23;
A_R = 8.06;

##
##  Step 3.2: Relative Frequencies and Nelder-Mead
##

table      = tabulate(final_vel_proc, [0.0:0.1:35.0]);
table(:,5) = table(:,2)./table(350,4);
x_weibull  = fminsearch(@(x) ( sum( ( 0.1*(x(1)/x(2))*(table(:,1)/x(2)).^(x(1)-1).*exp(-(table(:,1)/x(2)).^(x(1))) - table(:,5) ).^2 ) ), [k_R;A_R]);
figure(3);
plot(table(:,1), table(:,5), "color", "black")
lt03 = title("Example 3: Wind Speed Data vs. Probability Density Function", "fontsize", 16);
xt03 = xlabel("v", "fontsize", 14);
yt03 = ylabel("p", "fontsize", 14);
hold on
plot(table(:,1), 0.1*(k_R/A_R)*(table(:,1)/A_R).^(k_R-1).*exp(-(table(:,1)/A_R).^k_R), "color", "blue")
hold on
plot(table(:,1), 0.1*(x_weibull(1)/x_weibull(2))*(table(:,1)/x_weibull(2)).^(x_weibull(1)-1).*exp(-(table(:,1)/x_weibull(2)).^(x_weibull(1))), "color",  "red")

###
### Step 4: Calculating Annual Power Output Generation
###

##
##  Step 4.1: Setup of Scale and Shape Parameters
##

k_O = x_weibull(1);
A_O = x_weibull(2);

##
##  Step 4.2: Preparation
##

vel_av = [0:0.1:25];
p_R    = WeibullPDF(vel_av,A_R,k_R);
p_O    = WeibullPDF(vel_av,A_O,k_O);

##
##  Step 4.3: Annual Power Output Generation
##            With Logistic SSD Regression
##

pow_av = VectorialLogPowerfunction(vel_av,x_logSSD);

# Semi-Empirical Case

pow_semi        = VectorialLogPowerfunction(final_vel_proc,x_logSSD)';
sum_pow_semi    = sum(pow_semi);
year_semi_powW  = sum_pow_semi/(length(pow_semi))*24*365
year_semi_powGW = year_semi_powW/(1000000)

# Case of R-Weibull Parameters

pdfpow_av_R   = pow_av.*p_R;
avhpowout_R   = sum(0.1*pdfpow_av_R);
avypowoutW_R  = avhpowout_R*24*365
avypowoutGW_R = avypowoutW_R/(1000000)

# Case of Octave-Weibull Parameters

pdfpow_av_O   = pow_av.*p_O;
avhpowout_O   = sum(0.1*pdfpow_av_O);
avypowoutW_O  = avhpowout_O*24*365*1000
avypowoutGW_O = avypowoutW_O/(1000000000)