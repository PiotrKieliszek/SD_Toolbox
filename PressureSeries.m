% Shock and Detonation Toolboox
% http://www.galcit.caltech.edu/EDL/public/cantera/html/SD_Toolbox/
%
% This is a demostration of how to vary the initial pressure in a loop for
% a constant volume explosions and the ZND detonation simulations
%
% Using this demo as a guide, users can design their own loops and calulations

clear; clc;
display('Initial Pressure Series')

mech = 'h2air_highT.cti';   
gas = importPhase(mech); 
gas1 = importPhase(mech);
nsp = nSpecies(gas);
P1 = [];

% find Hydrogen nitrogen, and oxygen indices
ih2 = speciesIndex(gas,'H2');
io2  = speciesIndex(gas,'O2');
in2  = speciesIndex(gas,'N2');
x = zeros(nsp, 1);
x(ih2,1) = 2.0;
x(io2,1) = 1.0;
x(in2,1) = 3.76;
T1 = 300;
Po = 100000;
% fig_num & fname are for 'znd' - use '0' for no output
fig_num = 0;
fname = 0;
% plots: 1 = make plots, otherwise no plots
plots = 1;
display('Initial Conditions')
display('   Equivalence Ratio = 1, Temperature = 300 K')

npoints=15;
   disp(['For ', num2str(npoints), ' values of P1'])
for i = 1:npoints
   P1(i) = Po*(0.1 +1.1/npoints*(i-1));
   disp([' ', num2str(i),  ': P1 = ', num2str(P1(i))])
   set(gas,'Temperature',T1,'Pressure',P1(i),'MoleFractions',x);
   
   %%%Constant Volume Explosion Data%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   [cj_speed(i)] = CJspeed(P1(i),T1,x,mech,0); 
   set(gas,'Temperature',T1,'Pressure',P1(i),'MoleFractions',x);
   [gas] = PostShock_fr(cj_speed(i), P1(i), T1, x, mech);
   % SOLVE CONSTANT VOLUME EXPLOSION ODES
   [CVout] = explosion(gas,fig_num);
   ind_time_CV(i) = CVout.ind_time;
   exo_time_CV(i) = CVout.exo_time;   

   %%%%%ZND Detonation Data%%%%%
   % FIND POST SHOCK STATE FOR GIVEN SPEED
   set(gas1, 'T', T1, 'P', P1(i), 'X', x);
   [gas] = PostShock_fr(cj_speed(i), P1(i), T1, x, mech);
   Ts(i) = temperature(gas); %frozen shock temperature   
   Ps(i) = pressure(gas); %frozen shock pressure
   % SOLVE ZND DETONATION ODES
   [ZNDout] = znd(gas,gas1,fig_num,cj_speed(i),fname);
   ind_time_ZND(i) = ZNDout.ind_time_ZND;
   ind_len_ZND(i) = ZNDout.ind_len_ZND;
   exo_time_ZND(i) = ZNDout.exo_time_ZND;
   exo_len_ZND(i) = ZNDout.exo_len_ZND;
   tsteps = size(ZNDout.T,2);
   Tf_ZND(i) = ZNDout.T(tsteps);
   
   %%Calculate CJstate Properties%%%
   [gas] = PostShock_eq(cj_speed(i),P1(i),T1,x, mech);
   T2(i) = temperature(gas);
   P2(i) = pressure(gas);
   rho2(i) = density(gas);

   %Approximate the effective activation energy using finite differences
    Ta = Ts(i)*(1.02);
    set(gas, 'T', Ta, 'P', Ps(i), 'X', x);
    [CVout1] = explosion(gas,0);
    Tb = Ts(i)*(0.98);
    set(gas, 'T', Tb, 'P', Ps(i), 'X', x);
    [CVout2] = explosion(gas,0); 
    taua = CVout1.ind_time;
    taub = CVout2.ind_time;
    % Approximate effective activation energy for CV explosion
    if(taua==0 || taub==0)
        theta_effective_CV(i) = 0;
    else
        theta_effective_CV(i) = 1/Ts(i)*((log(taua)-log(taub))/((1/Ta)-(1/Tb))); 
    end
      
    [gas] = PostShock_fr(cj_speed(i)*1.02,P1(i),T1,x, mech);
    Ta= temperature(gas);
    [ZNDout1] = znd(gas,gas1,fig_num,cj_speed(i)*1.02,fname);
    [gas] = PostShock_fr(cj_speed(i)*0.98,P1(i),T1,x, mech);
    Tb = temperature(gas);
    [ZNDout2] = znd(gas,gas1,fig_num,cj_speed(i)*0.98,fname);
    ind_lena = ZNDout1.ind_len_ZND;
    ind_lenb = ZNDout2.ind_len_ZND;
    % Approximate effective activation energy for ZND Detonation
    if(ind_lena==0 || ind_lenb==0)
        theta_effective_ZND(i) = 0;
    else
        theta_effective_ZND(i) = 1/Ts(i)*((log(ind_lena)-log(ind_lenb))/((1/Ta)-(1/Tb))); 
    end
end

if(plots==1)
    % make plots
    close all;
    fontsize=12;
    figure;
    subplot(2,2,1);
    plot(P1(:),T2(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Temperature (K)','FontSize',fontsize);
    title('Post CJ State Temperature','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,2);
    plot(P1(:),P2(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Pressure (Pa)','FontSize',fontsize);
    title('Post CJ State Pressure','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,3);
    plot(P1(:),rho2(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Density (kg/m^3)','FontSize',fontsize);
    title('Post CJ State Density','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,4);
    plot(P1(:),cj_speed(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Ucj (m/s)','FontSize',fontsize);
    title('CJ Speed','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    %%% Plots for the Induction Zone (times and lengths)
    figure;
    subplot(2,2,1);
    plot(P1(:),ind_time_CV(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('tau_{CV_i} (s)','FontSize',fontsize);
    title('Induction time for CV explosion','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,2);
    plot(P1(:),ind_time_ZND(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('tau_{ZND_i} (s)','FontSize',fontsize);
    title('Induction time for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,4);
    plot(P1(:),ind_len_ZND(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Delta_{ZND_i} (m)','FontSize',fontsize);
    title('Induction length for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    %%% Plots for the Exothermic Zone (times and lengths)
    figure;
    subplot(2,2,1);
    plot(P1(:),exo_time_CV(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('tau_{CV_e} (s)','FontSize',fontsize);
    title('Exothermic time for CV explosion','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(2,2,2);
    plot(P1(:),exo_time_ZND(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('tau_{ZND_e} (s)','FontSize',fontsize);
    title('Exothermic time for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);
    
    subplot(2,2,4);
    plot(P1(:),exo_len_ZND(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Delta_{ZND_e} (m)','FontSize',fontsize);
    title('Exothermic length for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    %Ts and Tf for ZND detonation
    figure;
    subplot(1,2,1);
    plot(P1(:),Ts(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('T_s (K)','FontSize',fontsize);
    title('Frozen shock Temperature for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    subplot(1,2,2);
    plot(P1(:),Tf_ZND(:),'k');
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('T_f (K)','FontSize',fontsize);
    title('Final Temperature for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    %Approximation of the effective activation energy for CV explosion
    figure;
    plot(P1, theta_effective_CV);
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Theta_{CV} (J)','FontSize',fontsize);
    title('Effective activation energy for CV explosion','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);

    %Approximation of the effective activation energy for ZND detonation
    figure;
    plot(P1, theta_effective_ZND);
    xlabel('Initial Pressure (Pa)','FontSize',fontsize);
    ylabel('Theta_{ZND} (J)','FontSize',fontsize);
    title('Effective activation energy for ZND detonation','FontSize',fontsize);
    set(gca,'FontSize',12,'LineWidth',2);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE OUTPUT TEXT FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fn = ['Pressure_Series.plt'];
d = date;
P = P1/oneatm;
	
fid = fopen(fn, 'w');
fprintf(fid, '# ZND DETONATION STRUCTURE CALCULATION AND CONSTANT VOLUME EXPLOSION\n');
fprintf(fid, '# CALCULATION RUN ON %s\n\n', d);
	
fprintf(fid, '# INITIAL CONDITIONS\n');
fprintf(fid, '# TEMPERATURE (K) %4.1f\n', T1);
%fprintf(fid, '# PRESSURE (ATM) %2.1f\n', P);
fprintf(fid, '# DENSITY (KG/M^3) %1.4e\n', density(gas1));
q = 'H2:2 O2:1 N2:3.76';    
fprintf(fid, ['# SPECIES MOLE FRACTIONS: ' q '\n']);

%fprintf(fid, '# REACTION ZONE STRUCTURE\n\n');

fprintf(fid, '# THE OUTPUT DATA COLUMNS ARE:\n');
fprintf(fid, 'Variables = "Pressure state 1", "Temperature state 2", "Pressure state 2", "density state 2", "CJ Speed", "Temperature Post Shock", "Pressure Post Shock", "Induction time CV", "Exothermic time CV",  "Effective Activation Energy CV ",  "Effective Activation Energy ZND ", "Induction  time ZND", "Induction length ZND", "Exothermic time ZND", "Exothermic length ZND", "Final Temperature ZND"\n');

z = [P1; T2; P2; rho2; cj_speed; Ts; Ps; ind_time_CV; exo_time_CV; theta_effective_CV;  theta_effective_ZND; ind_time_ZND; ind_len_ZND; exo_time_ZND; exo_len_ZND; Tf_ZND];
fprintf(fid, '%14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e \t %14.5e  \t %14.5e  \t %14.5e  \t %14.5e  \t %14.5e  \t %14.5e\t %14.5e \n', z);

fclose(fid);
