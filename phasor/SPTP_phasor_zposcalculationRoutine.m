function zposlist = SPTP_phasor_zposcalculationRoutine(distX, distY, curveX, curveY)
%% Calculation of z positions from calibration data and distance in X,Y
%Required input: distX, distY: arrays with distance in X,Y
%               curveX, curveY: calibration curves in X,Y (SPTP_phasor_calibrationRoutine output)

%First remove values if they're not larger than their counterpart (i.e. X
%larger than Y or viceversa)
distXorig = distX;
distYorig = distY;
distX(distXorig<distYorig) = 0;
distY(distYorig<distXorig) = 0;
%also make boolean list of distX/distY - if they have a value, get 1,
%otherwise 0.
distXboolean = zeros(size(distX));
distXboolean(distX>0) = 1;
distYboolean = zeros(size(distY));
distYboolean(distY>0) = 1;
%now calculate zpos from these values. Solve 3rd order polynomial for this
a = curveX.Coefficients{1,1};
b = curveX.Coefficients{2,1};
c = curveX.Coefficients{3,1};
d = curveX.Coefficients{4,1};
y = distX;
zpos_distX = (sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) - 27.*a.^2.*d + ...
    27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)./(3.*2.^(1./3).*a) - (2.^(1./3).*(3.*a.*c - b.^2))./...
    (3.*a.*(sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) - ...
    27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)) - b./(3.*a);

a = curveY.Coefficients{1,1};
b = curveY.Coefficients{2,1};
c = curveY.Coefficients{3,1};
d = curveY.Coefficients{4,1};
y = distY;
zpos_distY = (sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) - 27.*a.^2.*d + ...
    27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)./(3.*2.^(1./3).*a) - (2.^(1./3).*(3.*a.*c - b.^2))./...
    (3.*a.*(sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) - ...
    27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)) - b./(3.*a);

zposlist = distXboolean.*zpos_distX+distYboolean.*zpos_distY;

end