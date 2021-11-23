% Traveltime inversion using raytracing and Gauss-Newton Method; R. Drews.

% Description:
%----------------------------------------------------------------------------------------------
% This code exemplifies how to invert radar-wide angle data using raytracing and Gauss-Newton. 
% The shape of the density profile is prescribed as a function of:
%
% density(depth) = 910 - A*exp(-r*depth)
%
% and the parameters A,r are inverted for together with reflector depths. This is a synthetic examples,  
% it can quickly be adapted for real data (real-case example should appear soon, contact me if not).
% Refer to the following paper for details (citations are appreciated):
%
% "Density anomaly in ice-shelf channels as seen by optical televiewing and wide-angle radar"
%  R. Drews, J. Brown, K. Matsuoka, E. Witrant, M. Philippe, B. Hubbard, F. Pattyn; 
%  submitted to The Cryosphere Discussions
%
%
% Code has been tested on Matlab 2011b.
%----------------------------------------------------------------------------------------------


% Requirements: 
% ---------------------------------------------------------------------------------------------
% The raytracing forward model is based on 
% Margrave, G. F.: Numerical Methods of Exploration Seismology with algorithms in Matlab, CREWES
% Toolbox Version: 1006 and requires the following files in a local folder "raytracing": 
%
% drawray.m  flipy.m  shootray.m  sphdiv.m  surround.m  traceray_pp.m
%
% These can be downloaded here: http://www.codeforge.com/article/174599
%----------------------------------------------------------------------------------------------
 

% To Do:
%---------------------------------------------------------------------------------------------
% 1. Get Graphical User Interface to make it more accessible
% 2. Use ASCII txt files for input data and unify scripts for real and syntethic cases
% 3. Ramping of Layers (to simulate systematic bias, has currently hardcoded shape of slope.)
% 3. Work on better documentation.
% 4. Transfer everything to Python (GUI as well?) 

% Last Changes: 28.09.2015 
%----------------------------------------------------------------------------------------------

% How to Use: The script is documented in file and hopefully somewhat self-explanatory. 



% Ok Let's Go!
%----------------------------------------------------------------------------------------------
clear all;close all;				%Get rid of everything that happened before.
addpath('./raytracing/')			%Add functions of Raytracing forward model


% Define Recording Geometry (CMP, WARR,..)
%----------------------------------------------------------------------------------------------------------
maxDepth=100;Nz=0.25;                           %Maximum Depth and depth discretization of background medium for Raytracing
zsrc=0;zrec=0;                                  %Depth of source and receivers (=0 when at  surface)
xoff = 0:Nz:30;                                 %Location of Geophones from minium to maximum offset 


%What to solve for? Currently implemented are (no error checking):
%SolveFor_A=1;SolveFor_Depth=1;SolveFor_r=1
%SolveFor_A=0;SolveFor_Depth=1;SolveFor_r=1 	a_true should be aguess; (DEFAULT)
%SolveFor_A=0;SolveFor_Depth=1;SolveFor_r=0	a_true should be aguess and r_true should be r_guess
%SolveFor_A=0;SolveFor_Depth=0;SolveFor_r=1     a_true should be aguess and reflector_depths_guess should be reflector_depth
%----------------------------------------------------------------------------------------------------------
SolveFor_A=1;SolveFor_Depth=1;SolveFor_r=1;


% Define Parameters of Forward Model and Inversion (careful, no strict error checking here. Make sure that 
% that if a reflector at reflector_depths=[] is addeed/removed, to also adjust sizes of reflector_depth_guess,
% InternalLayersForConversion and ramp_layers accordingly. 
%----------------------------------------------------------------------------------------------------------
atrue=460;rtrue=0.033;                          %True values for parameter A,r in Depth-Density Model 
aguess=460;rguess=0.015;                        %Initial Guess for parameters A,r in Depth-Density Model
reflector_depths = [5 15 30 50 70];             %True reflector depths
reflector_depths_guess = [10 22 38 41 85];      %Guessed reflector depth
InternalLayersForInversion = [1 0 1 1 1];       %Length determins number of layers, 1 inversion, 0 control; Don't put 0 everywhere.

sigma_i = 8e-8;		                      	%Estimated Uncertainty on Time measurement (in C_v matrix)
depth_sigma = 10;                               %Estimated Uncertainty for Reflector Depths (in m)
surface_density_sigma = 30;                     %Uncertainty of surface density in depth-density model, used for damping)
r_density_sigma = 0.01;                         %Uncertainty for constant r in depth-density model, used for damping)
MaxInvIt = 10;                                  %>1: Maximum number of iterations in adjusting model parameters. Usually only few are needed.
					        %Important is that J does not change anymore (see outuput plots)

                                                %Default Parameters for Raytracing:
caprad=0.1;itermax=20;                          %Capture Radius, maxnumber of iterations
pfan=-1;optflag=1;pflag=1;dflag=0;              %Default ray fan, and various flags


% Add some Noise and Bias to the synthetic data. Choose 0 for ideal cases
%----------------------------------------------------------------------------------------------------------
NoiseScaling =0.1.*1e-08;              %Add normally distributed noise to all reflectors
ramp_layers = [0.0 0 0 0 0]*1e-8;      %Pretend that layers are not only offset in depth, but also slanted (slope hardcoded..)
                                       %add float (e.g. [1.1 -2 0 0 -]to systematically perturb individual reflectors.


%Constants
%----------------------------------------------------------------------------------------------------------
c=3e8;vice=1.68e8;rho_ice=910;rho_air=2;
konst = (c/vice-1)/rho_ice;

% Regularization ?
%----------------------------------------------------------------------------------------------------------
lambdaNL=0e1;                          %this is the Lagrange Multiplier Lambda, choose 0 for ideal cases

% Do You want some plots as output? 
%----------------------------------------------------------------------------------------------------------
plot_initial_guess = 1;
plot_J_during_iteration = 1;
plot_final_results = 1;

%!!!! No more Setting of Parameters Needed Below Here!!!
%----------------------------------------------------------------------------------------------------------

error_checking


%Set weighting matrices and model parameters as function of what is inverted for
%----------------------------------------------------------------------------------------------------------
NumberOfLayersForInversion=length(find(InternalLayersForInversion==1));
if  (SolveFor_A==1 & SolveFor_Depth==1 & SolveFor_r==1)
    sigma_j = [surface_density_sigma r_density_sigma (1:NumberOfLayersForInversion)*0+depth_sigma];
    Cm = diag(sigma_j.^2);
    NumberOfModelParameters = NumberOfLayersForInversion +2 ;
elseif (SolveFor_A==0 & SolveFor_Depth==1 & SolveFor_r==1)
    sigma_j = [r_density_sigma (1:NumberOfLayersForInversion)*0+depth_sigma];
    Cm = diag(sigma_j.^2);
    NumberOfModelParameters = NumberOfLayersForInversion +1;
elseif (SolveFor_A==0 & SolveFor_Depth==1 & SolveFor_r==0)
    sigma_j = [(1:NumberOfLayersForInversion)*0+depth_sigma];
    Cm = diag(sigma_j.^2);
    NumberOfModelParameters = NumberOfLayersForInversion;
elseif (SolveFor_A==0 & SolveFor_Depth==0 & SolveFor_r==1)
    sigma_j = [r_density_sigma];
    Cm = diag(sigma_j.^2);
    NumberOfModelParameters = 1;
end

% Raycolors for plotting, and depth vector
%----------------------------------------------------------------------------------------------------------
zp=0:Nz:maxDepth;
raycolors = ['r','k','g','c','m','b','y'];
xoffall=xoff;


%True and Guessed Density/Velocity profiles
%------------------------------------------
rhoz_true = 910 -atrue*exp(-rtrue*zp);
rhoz_guess = 910 -aguess*exp(-rguess*zp);
vztrue = c./(rhoz_true*konst+1);vzguess = c./(rhoz_guess*konst+1);
vzguess = c./(rhoz_guess*konst+1);

%101. Create Data
%--------------------------------------------------------------------------
verbose=0;
[DataLayers,ControlLayers] = CreateSyntheticData(InternalLayersForInversion,vztrue,zp,zsrc,zrec,reflector_depths,reflector_depths_guess,NoiseScaling,ramp_layers,xoff,xoffall,caprad,pfan,itermax,optflag,pflag,verbose);
NumberOfLayersForInversion = length(DataLayers.tdata);
NumberOfLayersForControl = length(ControlLayers.tdata);
%102 Done. DataLayers created, corresponding velocity model is vztrue(zp)
%--------------------------------------------------------------------------

%Create some vectors/matrices based on the input provided above
%--------------------------------------------------------------------------
InternalLayersForControl =  InternalLayersForInversion~=1;
it=1;aguess_it(it)=aguess;rguess_it(it)=rguess;
%--------------------------------------------------------------------------

%Initial Guess
%--------------------------------------------------------------------------
[SimulatedLayers,delta_t,NumberOfDataPoints] = InitGuess(DataLayers,vzguess,it,zp,zsrc,zrec,NoiseScaling,xoff,xoffall,caprad,pfan,itermax,optflag,pflag,plot_initial_guess);
%Get misfit of initial guess
delta_t_rms(1) = sqrt(sum(delta_t.^2))/length(delta_t);
%--------------------------------------------------------------------------

%Initialize some variables which will be recording parameters during loops
%--------------------------------------------------------------------------
for k =1 :NumberOfLayersForInversion
    reflector_depths_updated{1}(k)=DataLayers.depth_guess{k};
    reflector_depths_updated_tmp{1}(k)=DataLayers.depth_guess{k};
end
vzupdated{1} = vzguess; 
A = zeros(NumberOfDataPoints,NumberOfModelParameters);
Ct = diag((1:NumberOfDataPoints)*0+sigma_i.^2);
number_of_iterations = 0;
if (SolveFor_A==1 & SolveFor_Depth==1 & SolveFor_r==1)
    modelp =[aguess_it(1);rguess_it(1);reflector_depths_updated{1}(:)];
elseif (SolveFor_A==0 & SolveFor_Depth==1 & SolveFor_r==1)
    modelp =[rguess_it(1);reflector_depths_updated{1}(:)];
elseif (SolveFor_A==0 & SolveFor_Depth==1 & SolveFor_r==0)
    modelp =[reflector_depths_updated{1}(:)];
elseif (SolveFor_A==0 & SolveFor_Depth==0 & SolveFor_r==1)
    modelp =[rguess_it(1)];
end
modelp_0=modelp; %initial guess used in regularization term. 
%--------------------------------------------------------------------------


%Start iteration to find updated parameters
%--------------------------------------------------------------------------
for it = 2:MaxInvIt
    number_of_iterations = number_of_iterations +1;

    %Calculate Partial Derivatives
    %-------------------------------------------
    [A]=GetSensitivityMatrix(zeros(NumberOfDataPoints,NumberOfModelParameters),NumberOfLayersForInversion,SimulatedLayers,SolveFor_A,SolveFor_Depth,SolveFor_r,zp,vzupdated,rguess_it,aguess_it,it);
   
    %Setup for J = 1/2 ||t_mod - t_mes||^2 + 1/2 gm'gm
    %          J = 1/2 dt'*Ct*dt + 1/2gm'*Cm*gm 
    %-------------------------------------------------
    %Case: gm = lambda m'm
    gm = (modelp-modelp_0);
    dgm_dm = eye(length(modelp));
    gradJ = delta_t'*inv(Ct)*A + (lambdaNL*gm'*inv(Cm)*dgm_dm);
    PsiJ = A'*inv(Ct)*A+lambdaNL*inv(Cm);
    update = inv(PsiJ)*gradJ';
    modelp = modelp-update;

    %Parameter Update for this iteration
    %-------------------------------------------
    if (SolveFor_A==1 & SolveFor_Depth==1 & SolveFor_r==1)
        aguess_it(it) = modelp(1);
        rguess_it(it) = modelp(2);
        rhoz_guess =  910 - aguess_it(it)*exp(-rguess_it(it)*zp);
        vzupdated{it} =  c./(rhoz_guess*konst+1);
        for k=1:NumberOfLayersForInversion
            reflector_depths_updated{it}(k) = modelp(2+k);
        end
    elseif (SolveFor_A==0 & SolveFor_Depth==1 & SolveFor_r==1)
        aguess_it(it) = aguess;
        rguess_it(it) = modelp(1);
        rhoz_guess =  910 - aguess*exp(-rguess_it(it)*zp);
        vzupdated{it} =  c./(rhoz_guess*konst+1);
        for k=1:NumberOfLayersForInversion
            reflector_depths_updated{it}(k) = modelp(1+k);
        end
    elseif (SolveFor_A==0 & SolveFor_Depth==0 & SolveFor_r==1)
        aguess_it(it) = aguess;
        rguess_it(it) = modelp(1);
        rhoz_guess =  910 - aguess*exp(-rguess_it(it)*zp);
        vzupdated{it} =  c./(rhoz_guess*konst+1);
        for k=1:NumberOfLayersForInversion
            %reflector_depths_updated{it}(k) = DataLayers.vrms_depth{k};
            reflector_depths_updated{it}(k) = DataLayers.depth_guess{k};
        end
    elseif (SolveFor_A==0 & SolveFor_Depth==1 & SolveFor_r==0)
        aguess_it(it) = aguess;
        rguess_it(it) = rguess;
        rhoz_guess =  910 - aguess*exp(-rguess_it(it)*zp);
        vzupdated{it} =  c./(rhoz_guess*konst+1);
        for k=1:NumberOfLayersForInversion
            reflector_depths_updated{it}(k) = modelp(k);
        end
    end
  
        
    %Compute updated rays and distance vector
    delta_t=[];
    for k=1:NumberOfLayersForInversion
        %Simulate final results:
        pfan=-1;optflag=1;pflag=1;dflag=0;% default ray fan, and variousflags
        [tdata,pdata,ldata,rc_data]=traceray_pp(vzupdated{it},zp,zsrc,zrec,reflector_depths_updated{it}(k),DataLayers.xoff{k}',caprad,pfan,itermax,optflag,pflag,dflag,'k');
        %All rays traced?
        if (sum(isinf(tdata))>0)
            display('Not all rays could be traced successfully. Probably initial guess is too far off. Also check capture radius as function of geometry. Or reflectors too shallow for given offset? Stop here.')
            break;
        end
        SimulatedLayers.rc{k,it} = rc_data;
        SimulatedLayers.tdata{k,it} = tdata;
        SimulatedLayers.xoff{k,it} = DataLayers.xoff{k};
        %Evaluate distance error of updated parameters
        delta_t = [delta_t; [SimulatedLayers.tdata{k,it}-DataLayers.tdata{k}']'];
    end
 
  
    delta_t_rms(it) = sqrt(sum(delta_t.^2))/length(delta_t);
   
        
    %Update J
    J(it-1) = 0.5*delta_t'*inv(Ct)*delta_t+0.5*lambdaNL*(modelp'-modelp_0')*inv(Cm)*(modelp-modelp_0);
     
    %Plot J (should decrease each iteration, otherwise something is wrong (e.g. increase lambda))
    if (plot_J_during_iteration)
    	figure(2)
    	subplot(1,2,1)
    	plot(it,J(it-1),'rx'); hold on;
	    xlabel('Number of Iteration')
	    ylabel('J')
    	subplot(1,2,2)
    	plot(1,log(J(1)),'rx');
    	plot(it,log(J(it-1)),'rx'); hold on;
        xlabel('Number of Iteration')
        ylabel('Log(J)')
    end
end % End iterating for Inversion


 
%Do the control with leftover layers if there are any
%---------------------------------------------------------------------------
for k =1:NumberOfLayersForControl
    [tdata,pdata,ldata,rc_data]=traceray_pp(vzupdated{end},zp,zsrc,zrec,ControlLayers.depth{k},ControlLayers.xoff{k}',caprad,pfan,itermax,optflag,pflag,dflag,'k');
    ControlLayers.rc{k} = rc_data;
    ControlLayers.sim_tdata{k} = tdata;
end

%Calculate Initial and Final densities
%---------------------------------------------------------------------------
 rhoz_guess = 910 -aguess_it(1)*exp(-rguess_it(1)*zp);
 HA_Crim_guess = maxDepth*(mean(rhoz_guess)-rho_ice)/(rho_air-rho_ice);
 HA_Crim_true = maxDepth*(mean(rhoz_true)-rho_ice)/(rho_air-rho_ice);
 
if (plot_final_results)
   f3=figure(3);
   set(f3, 'Position',  [300, 150, 1200, 400]); 
   subplot(1,5,1)
   plot(vztrue,zp,'g-');hold on;
   for iterations=1:number_of_iterations
    plot(vzupdated{iterations},zp,'k--');xlabel('Velocity (m/s)');ylabel('Depth (m)');
   end
   xlim([1.67e8,1.99e8]);legend('True','Inverted','Location','SouthEast');flipy;
   subplot(1,5,2)
   plot(1:number_of_iterations+1,log10(delta_t_rms),'r-x'); hold on;
   xlabel('Number of Iteration');ylabel('Log (\Delta t_{rms}) (s)');title('Traveltime residuals')
   subplot(1,5,3)
   plot(1:number_of_iterations+1,aguess_it(:)-atrue,strcat(raycolors(k),'-x'),'markersize',13); hold on;
   xlabel('Number of Iteration');ylabel('Update Parameter: A (kg/m^3) ');title('Update for surface density')
   subplot(1,5,4)
   plot(1:number_of_iterations+1,rguess_it(:)-rtrue,strcat(raycolors(k),'-x'),'markersize',13); hold on;
   xlabel('Number of Iteration');ylabel('Update Parameter: r (m^{-1}) ');title('Update for densification length')
   subplot(1,5,5)
   for k=1:NumberOfLayersForInversion
       for it=1:number_of_iterations+1
        plot(it,reflector_depths_updated{it}(k)-DataLayers.depth{k},strcat(raycolors(k),'-x'),'markersize',13); hold on;
       end
   end
   xlabel('Number of Iteration');ylabel('Update Reflector Depths (m)'); title('Update of Reflector Depths');
   figure(4);
   for k=1:NumberOfLayersForInversion
        subplot(1,2,1)
        plot(DataLayers.xoff{k},DataLayers.tdata{k},'k-x');grid;hold on;
        plot(SimulatedLayers.xoff{k,number_of_iterations+1},SimulatedLayers.tdata{k,number_of_iterations+1},'r-')
        xlabel('meters');ylabel('seconds');
        subplot(1,2,2)
        plot(DataLayers.xoff{k},DataLayers.tdata{k}-SimulatedLayers.tdata{k,number_of_iterations+1}','kx');hold on;
        xlabel('Offset (m)');ylabel('Traveltime residual (\Delta t_{rms})')
        title('Travel time residuals')
   end
   if (NumberOfLayersForControl>0)
    for k =1:NumberOfLayersForControl
        subplot(1,2,1)
        plot(ControlLayers.xoff{k},ControlLayers.tdata{k},'bx');grid;hold on;
        plot(ControlLayers.xoff{k},ControlLayers.sim_tdata{k},'b-')
        xlabel('meters');ylabel('seconds');flipy;
        subplot(1,2,2)
        plot(ControlLayers.xoff{k},ControlLayers.tdata{k}-ControlLayers.sim_tdata{k}','bo')
    end
   end
   
   title('Final Fit (blue is control, red is inverted)')
   
   figure(5)
   subplot(1,2,1)
   rhoz_final = 910-aguess_it(end)*exp(-rguess_it(end)*zp);
   HA_Crim_inv = maxDepth*(mean(rhoz_final)-rho_ice)/(rho_air-rho_ice);
   plot(rhoz_final,zp,'r');hold on;
   plot(rhoz_true,zp,'g');flipy;
   title('Final and true density')
   legend('Final Density','True Density','Location','southwest')
   subplot(1,2,2)
   plot(rhoz_final-rhoz_true,zp,'r')
   flipy;xlabel('Depth (m)');ylabel('Density (kg/m^3)');
   title('Difference between final and true density')
end
  


%Some Output
display('-------------------------------------------------------')
display('Output from Inversion:')
display('-------------------------------------------------------')
display(strcat('Initial  Depth (m): ',num2str(reflector_depths_guess(find(InternalLayersForInversion==1)))))
display(strcat('Inverted Depth (m): ',regexprep(num2str(reflector_depths_updated{end}(:)'), '\s+', ' ')))
display(strcat('True Depth (m)    : ',regexprep(num2str(reflector_depths(find(InternalLayersForInversion==1))), '\s+', ' ')))
display(strcat('Initial A  : ',num2str(aguess)))
display(strcat('Inverted A : ',num2str(aguess_it(end))))
display(strcat('True A     : ',num2str(atrue)))
display(strcat('Initial r  : ',num2str(rguess)))
display(strcat('Inverted r : ',num2str(rguess_it(end))))
display(strcat('True r     : ',num2str(rtrue)))
display(strcat('Guessed firn-air content  : ',num2str(HA_Crim_guess)))
display(strcat('Inverted firn-air content : ',num2str(HA_Crim_inv)))
display(strcat('True firn-ar content      : ',num2str(HA_Crim_true)))
display(strcat('Mean depth av final velocity: ',num2str(mean(c./(rhoz_final*konst+1)))))

