function [DataLayers,ControlLayers]=CreateSyntheticData(InternalLayersForInversion,vztrue,zp,zsrc,zrec,reflector_depths,reflector_depths_guess,NoiseScaling,ramp_layers,xoff,xoffall,caprad,pfan,itermax,optflag,pflag,verbose)

%102. Create data 
%--------------------------------------------------------------------------

raycolors = ['r','k','c','r','m','b','y']; 
if (verbose == 1)
    %Do plot with rays traced in different colors.
	figure;subplot(2,1,1);flipy;
end
NumberOfLayersForInversion = 0;NumberOfLayersForControl = 0;
for k=1:length(InternalLayersForInversion)
    if k<=7
        kol = raycolors(k);
    else
        kol = raycolors(mod(k,7));
    end
    if (InternalLayersForInversion(k)==1)
        if verbose == 1
            dflag=2;subplot(2,1,1);
        end
        NumberOfLayersForInversion = NumberOfLayersForInversion+1;
        %Raytracing with true velocities
        dflag=2;
        [tdata,pdata,ldata,rc_data]=traceray_pp(vztrue,zp,zsrc,zrec,reflector_depths(k),xoff,caprad,pfan,itermax,optflag,pflag,dflag,kol);
        %Capture if not all rays could be traced
        if (sum(isinf(tdata))>0)
            display('Not all rays could be traced successfully. Probably initial guess is too far off. Also check capture radius as function of geometry. Stop here.')
            break;
        end
     
        %Add Noise
        tdata_id = tdata;
        tdata = tdata + randn(1,length(tdata))*NoiseScaling;
        %Add Ramp (to simulat unaccounted layer-dipping and accumulated positiong error)
        if ramp_layers(k) ~= 0
            lslope = 2*ramp_layers(k)/(xoff(end)-xoff(1));xp=(xoff(end)-xoff(1))/2;[atmp btmp]=min(xoff-xp);yp=0;cp=yp-lslope*xp;
            tramp = lslope*xoff+cp;
            tdata = tdata-tramp;
        end
           
        DataLayers.rc{NumberOfLayersForInversion} = rc_data;
        DataLayers.tdata{NumberOfLayersForInversion} = tdata';
        DataLayers.xoff{NumberOfLayersForInversion} = xoffall';
        DataLayers.depth_guess{NumberOfLayersForInversion} = reflector_depths_guess(k);
        DataLayers.depth{NumberOfLayersForInversion} = reflector_depths(k);
    else
        dflag=0;
        NumberOfLayersForControl = NumberOfLayersForControl+1;
        %Raytracing with true velocities
        [tdata,pdata,ldata,rc_data]=traceray_pp(vztrue,zp,zsrc,zrec,reflector_depths(k),xoff,caprad,pfan,itermax,optflag,pflag,dflag,kol);
        %Add Noise
        tdata_id = tdata;
        tdata = tdata + randn(1,length(tdata))*NoiseScaling;
        %Add Ramp
        if ramp_layers(k) ~= 0
            lslope = 2*ramp_layers(k)/(xoff(end)-xoff(1));xp=(xoff(end)-xoff(1))/2;[atmp btmp]=min(xoff-xp);yp=0;cp=yp-lslope*xp;
            tramp = lslope*xoff+cp;
            tdata = tdata-tramp;
        end 
        ControlLayers.rc{NumberOfLayersForControl} = rc_data;
        ControlLayers.tdata{NumberOfLayersForControl} = tdata';
        ControlLayers.xoff{NumberOfLayersForControl} = xoffall';
        ControlLayers.depth_guess{NumberOfLayersForControl} = reflector_depths_guess(k);
        ControlLayers.depth{NumberOfLayersForControl} = reflector_depths(k);
    end
end
if (verbose == 1)
    line(xoff,zrec*ones(size(xoff)),'color','b','linestyle','none','marker','v')
    line(0,zsrc,'color','r','linestyle','none','marker','*')
    grid;xlabel('meters');ylabel('meters');
    for k=1:NumberOfLayersForInversion
        subplot(2,1,2);plot(xoff,DataLayers.tdata{k},'r');grid;hold on;
        xlabel('meters');ylabel('seconds');
    end
    for k=1:NumberOfLayersForControl
        subplot(2,1,2);plot(xoff,ControlLayers.tdata{k},'b');grid;hold on;
        xlabel('meters');ylabel('seconds');
    end
    flipy;title('Created Data based on known interval velocities')
end
