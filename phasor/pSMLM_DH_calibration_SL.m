function [curveAngle curveDistAngle wobbleMatrix] = pSMLM_DH_calibration_SL(input,zposcali,imagegen,path_storage_mat)
%% Phasor double-helix calibration for SMALLLABS
%Edit from original calibration file:
%Save data as .mat file for reading
[curveAngle curveDistAngle DHmagnratio wobbleMatrix] = pSMLM_DH_calibration(input,zposcali,imagegen);
save(path_storage_mat,'curveAngle','curveDistAngle','wobbleMatrix','DHmagnratio');
end