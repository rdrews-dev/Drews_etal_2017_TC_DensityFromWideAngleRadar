function [S] = GetSensitivityMatrix(A,NumberOfLayersForInversion,SimulatedLayers,SolveFor_A,SolveFor_Depth,SolveFor_r,zp,vzupdated,rguess_it,aguess_it,it)


%Function which calculates all partial derivatives and assembles them in an
%Matrix with M columns (#number of model parameters) and N lines (#number
%of datapoints)
%Partial Derivatives (could be vecotrized if speed is an issue

c=3e8;vice=1.68e8;rho_ice=910;rho_air=2;
konst = (c/vice-1)/rho_ice;
S=A; %initialize Sensitivity Matrix with correct dimensions.

%Get number of Rays for Each Layer
for k=1:NumberOfLayersForInversion
    number_of_rays_perlayer(k) = length(SimulatedLayers.rc{k,it-1}');
end
cumsum_number_of_rays_perlayer = [0 cumsum(number_of_rays_perlayer)];


if (SolveFor_A==1 & SolveFor_Depth==1 & SolveFor_r==1)
    %Get Partial Derivatives for each layer
    for k=1:NumberOfLayersForInversion
              %Loop over each ray in Layer k
              for raynumber = 1:number_of_rays_perlayer(k);
                rc=SimulatedLayers.rc{k,it-1}{raynumber};       %local ray
                dl = (diff(rc(:,1)).^2+diff(rc(:,2)).^2).^0.5;  %local path length
                ind_reflection_local = find(dl==0);             %local reflection id 
                [tmp, ind_reflection_global] = min(abs(zp-rc(ind_reflection_local,2))); %global reflection id 

                %derivative dt_i/da 
                dt_i_da{k}(raynumber,it-1) = -konst/c*sum(exp(-rguess_it(it-1)*rc(2:end,2)).*dl);
                S(cumsum_number_of_rays_perlayer(k)+raynumber,1) = -konst/c*sum(exp(-rguess_it(it-1)*rc(2:end,2)).*dl);

                %derivative dt_i/dr_k
                dt_i_dr{k}(raynumber,it-1) = konst*aguess_it(it-1)/c*sum(rc(2:end,2).*exp(-rguess_it(it-1)*rc(2:end,2)).*dl);
                S(cumsum_number_of_rays_perlayer(k)+raynumber,2) = konst*aguess_it(it-1)/c*sum(rc(2:end,2).*exp(-rguess_it(it-1)*rc(2:end,2)).*dl);

                %derivative dt_i/dl (block matrixes, with delta_mn)
                for ll = 1:NumberOfLayersForInversion
                    if (k == ll )
                        cos_theta = (rc(ind_reflection_local,2)-rc(ind_reflection_local-1,2))/dl(ind_reflection_local-1);
                        dt_i_dl{k}(raynumber,it-1) = 2*cos_theta/vzupdated{it-1}(ind_reflection_global-2);
                        S(cumsum_number_of_rays_perlayer(k)+raynumber,2+ll) = 2*cos_theta/vzupdated{it-1}(ind_reflection_global-2);
                    end
                end
              end
    end
elseif (SolveFor_A==0 & SolveFor_Depth==1 & SolveFor_r==1)
    %Get Partial Derivatives for each layer
    for k=1:NumberOfLayersForInversion
              %Loop over each ray in Layer k
              for raynumber = 1:number_of_rays_perlayer(k);
                rc=SimulatedLayers.rc{k,it-1}{raynumber};       %local ray
                dl = (diff(rc(:,1)).^2+diff(rc(:,2)).^2).^0.5;  %local path length
                ind_reflection_local = find(dl==0);             %local reflection id 
                [tmp, ind_reflection_global] = min(abs(zp-rc(ind_reflection_local,2))); %global reflection id 

                %derivative dt_i/dr_k
                dt_i_dr{k}(raynumber,it-1) = konst*aguess_it(1)/c*sum(rc(2:end,2).*exp(-rguess_it(it-1)*rc(2:end,2)).*dl);
                S(cumsum_number_of_rays_perlayer(k)+raynumber,1) = konst*aguess_it(1)/c*sum(rc(2:end,2).*exp(-rguess_it(it-1)*rc(2:end,2)).*dl);

                %derivative dt_i/dl (block matrixes, with delta_mn)
                for ll = 1:NumberOfLayersForInversion
                    if (k == ll )
                        cos_theta = (rc(ind_reflection_local,2)-rc(ind_reflection_local-1,2))/dl(ind_reflection_local-1);
                        dt_i_dl{k}(raynumber,it-1) = 2*cos_theta/vzupdated{it-1}(ind_reflection_global-2);
                        S(cumsum_number_of_rays_perlayer(k)+raynumber,1+ll) = 2*cos_theta/vzupdated{it-1}(ind_reflection_global-2);
                    end
                end
              end
    end
elseif (SolveFor_A==0 & SolveFor_Depth==1 & SolveFor_r==0)
     %Get Partial Derivatives for each layer
    for k=1:NumberOfLayersForInversion
              %Loop over each ray in Layer k
              for raynumber = 1:number_of_rays_perlayer(k);
                rc=SimulatedLayers.rc{k,it-1}{raynumber};       %local ray
                dl = (diff(rc(:,1)).^2+diff(rc(:,2)).^2).^0.5;  %local path length
                ind_reflection_local = find(dl==0);             %local reflection id 
                [tmp, ind_reflection_global] = min(abs(zp-rc(ind_reflection_local,2))); %global reflection id 

                %derivative dt_i/dl (block matrixes, with delta_mn)
                for ll = 1:NumberOfLayersForInversion
                    if (k == ll )
                        cos_theta = (rc(ind_reflection_local,2)-rc(ind_reflection_local-1,2))/dl(ind_reflection_local-1);
                        dt_i_dl{k}(raynumber,it-1) = 2*cos_theta/vzupdated{it-1}(ind_reflection_global-2);
                        S(cumsum_number_of_rays_perlayer(k)+raynumber,ll) = 2*cos_theta/vzupdated{it-1}(ind_reflection_global-2);
                    end
                end
              end
    end
  elseif (SolveFor_A==0 & SolveFor_Depth==0 & SolveFor_r==1)
    %Get Partial Derivatives for each layer
    for k=1:NumberOfLayersForInversion
              %Loop over each ray in Layer k
              for raynumber = 1:number_of_rays_perlayer(k);
                rc=SimulatedLayers.rc{k,it-1}{raynumber};       %local ray
                dl = (diff(rc(:,1)).^2+diff(rc(:,2)).^2).^0.5;  %local path length
                ind_reflection_local = find(dl==0);             %local reflection id 
                [tmp, ind_reflection_global] = min(abs(zp-rc(ind_reflection_local,2))); %global reflection id 

                %derivative dt_i/dr_k
                dt_i_dr{k}(raynumber,it-1) = konst*aguess_it(1)/c*sum(rc(2:end,2).*exp(-rguess_it(it-1)*rc(2:end,2)).*dl);
                S(cumsum_number_of_rays_perlayer(k)+raynumber,1) = konst*aguess_it(1)/c*sum(rc(2:end,2).*exp(-rguess_it(it-1)*rc(2:end,2)).*dl);
              end
    end
end