function [SimulatedLayers,delta_t,NumberOfDataPoints] = InitGuess(DataLayers,vzguess,it,zp,zsrc,zrec,NoiseScaling,xoff,xoffall,caprad,pfan,itermax,optflag,pflag,plot_initial_guess)

delta_t=[];NumberOfDataPoints=0;
raycolors = ['r','k','g','c','m','b','y'];
NumberOfLayersForInversion = length(DataLayers.tdata);
if (plot_initial_guess == 1)	
	dflag=2;hFig=figure(1);subplot(1,2,1);flipy;
	set(hFig, 'Position', [200 200 1500 500])
else
    	dflag=0;
end
for k=1:NumberOfLayersForInversion
    if k<=7
        kol = raycolors(k);
    else
        kol = raycolors(mod(k,7));
    end
    %Raytracing 
    [tdata,pdata,ldata,rc_data]=traceray_pp(vzguess,zp,zsrc,zrec,DataLayers.depth_guess{k},DataLayers.xoff{k}',caprad,pfan,itermax,optflag,pflag,dflag,kol);
    %All rays traced?
    if (sum(isinf(tdata))>0)
       display('Not all rays could be traced successfully. Probably initial guess is too far off. Also check capture radius as function of geometry. Or reflectors too shallow for given offset? Stop here.')
       break;
    end
    SimulatedLayers.rc{k,it} = rc_data;
    SimulatedLayers.tdata{k,it} = tdata;
    delta_t = [delta_t, SimulatedLayers.tdata{k,it}-DataLayers.tdata{k}'];
    NumberOfDataPoints = NumberOfDataPoints + length(rc_data');
end
delta_t = delta_t'; %delta_t as as many lines as measurements
if (plot_initial_guess == 1)
    subplot(1,2,1)
    line(xoffall,zrec*ones(size(xoffall)),'color','b','linestyle','none','marker','v')
    line(0,zsrc,'color','r','linestyle','none','marker','*')
    grid;xlabel('Offset (m)');ylabel('Depth (m)');title('Initial Guess Raytracing');
    for k=1:NumberOfLayersForInversion
        if k<=7
            kol = raycolors(k);
        else
            kol = raycolors(mod(k,7));
        end
	figure(1);subplot(1,2,2);
        plot(DataLayers.xoff{k},DataLayers.tdata{k},strcat(kol,'x'));grid;hold on;
        plot(DataLayers.xoff{k},SimulatedLayers.tdata{k,1},kol)
    end
    xlabel('Offset (m)');ylabel('Traveltime (s)');title('Initial Guess: Curves are Model, crosses are Data');
    flipy;
end
