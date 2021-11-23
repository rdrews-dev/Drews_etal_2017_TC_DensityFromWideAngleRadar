%Very limited check if input parameters are useful:

% Loose error checking for input parameters
%----------------------------------------------------------------------------------------------------------
if (length(ramp_layers) ~= length(InternalLayersForInversion) || ...
    length(ramp_layers) ~= length(reflector_depths_guess) || ...
    length(ramp_layers) ~= length(reflector_depths) || ...
    length(reflector_depths_guess) ~= length(reflector_depths) || ...
    length(reflector_depths_guess) ~= length(InternalLayersForInversion) || ...
    length(InternalLayersForInversion) ~= length(reflector_depths))
    display('Warning: Ramp_Layers, InternalLayersForInversion,reflector_depths_guess,reflector_depths have not common sizes');
    display('Needs adjustment. Stop here.')
    break;
end

if sum(InternalLayersForInversion)==0
    display('Choose at least one Internal Layer for Inverison!')
    break;
end

if MaxInvIt<2
    display('Number of Iterations should be at least 2. MaxInvIt=1 is reserved for initial guess.')
    break;
end


if SolveFor_A==0 & (atrue ~= aguess)
    display('Warning: A is not inverted for but atrue =! aguess. Be aware.')
end

if SolveFor_r==0 & (rtrue ~= rguess)
    display('Warning: r is not inverted for but rtrue =! rguess. Be aware.')
end


if SolveFor_Depth==0 & (isequal(reflector_depths,reflector_depths))
    display('Warning: Depth is not inverted for but Depth Guess and True Depths are not the same. Be aware.')
end


if (maxDepth < max(reflector_depths) || maxDepth < max(reflector_depths_guess))
   display('Choose maxDepth larger than the Reflector depths (and Reflector depths guesses).')
   break;
end