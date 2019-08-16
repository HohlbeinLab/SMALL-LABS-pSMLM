function zpos = phasorcalcZpos(sigma,calibrationFile)
%Read calibration .mat file
% keyboard
load(calibrationFile);
%We now have 2 curves with a*x^3+b*x^2+c*x+d - fits for X/Y and Y/X
%We calculate this same as for SPTP
%First all X/Y
a = astig_phasor_curveXY.a;
b = astig_phasor_curveXY.b;
c = astig_phasor_curveXY.c;
d = astig_phasor_curveXY.d;
y = sigma(:,1)./sigma(:,2);
zpos_XY = (sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) - 27.*a.^2.*d + ...
    27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)./(3.*2.^(1./3).*a) - (2.^(1./3).*(3.*a.*c - b.^2))./...
    (3.*a.*(sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) - ...
    27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)) - b./(3.*a);
%Then all Y/X
a = astig_phasor_curveYX.a;
b = astig_phasor_curveYX.b;
c = astig_phasor_curveYX.c;
d = astig_phasor_curveYX.d;
y = sigma(:,2)./sigma(:,1);
zpos_YX = (sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) - 27.*a.^2.*d + ...
    27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)./(3.*2.^(1./3).*a) - (2.^(1./3).*(3.*a.*c - b.^2))./...
    (3.*a.*(sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) - ...
    27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)) - b./(3.*a);
%Then choose XY or YX
zpos = zpos_XY;
zpos(sigma(:,2)>sigma(:,1)) = zpos_YX(sigma(:,2)>sigma(:,1));
zpos = zpos*1000; %um to nm

%read calibration YAML file
% fileID = fopen(calibrationFile);
% caliData = fscanf(fileID,'%s');
% fclose(fileID);
% %find patterns
% a1pos = strfind(caliData,'a1:');
% a2pos = strfind(caliData,'a2:');
% anglepos = strfind(caliData,'angle:');
% b1pos = strfind(caliData,'b1:');
% b2pos = strfind(caliData,'b2:');
% c1pos = strfind(caliData,'c1:');
% c2pos = strfind(caliData,'c2:');
% d1pos = strfind(caliData,'d1:');
% d2pos = strfind(caliData,'d2:');
% a1 = str2num(caliData((a1pos+3):(a2pos-2)));
% a2 = str2num(caliData((a2pos+3):(anglepos-2)));
% b1 = str2num(caliData((b1pos+3):(b2pos-2)));
% b2 = str2num(caliData((b2pos+3):(c1pos-2)));
% c1 = str2num(caliData((c1pos+3):(c2pos-2)));
% c2 = str2num(caliData((c2pos+3):(d1pos-2)));

%calculate Z from sigma1/2, since this doesn't work via ImageJ/FIJI
%incorporation
%Direct translation from java code
%Loop over all localizations
% for i = 1:size(sigma,1)
%     if sigma(i,1) > sigma(i,2)
%         if c1>c2
%             zpos(i) = -1*sqrt(((sigma(i,1)/sigma(i,2))-b1)/a1)+c1;
%         else
%             zpos(i) = sqrt(((sigma(i,1)/sigma(i,2))-b2)/a2)+c2;
%         end
%     else
%         if c1>c2
%             zpos(i) = sqrt(((sigma(i,2)/sigma(i,1))-b2)/a2)+c2;
%         else
%             zpos(i) = -1*sqrt(((sigma(i,2)/sigma(i,1))-b2)/a2)+c2;
%         end
%     end
%     if isnan(zpos(i))
%         zpos(i) = 0;
%     end
% end
                    


%Below: old method, should be wrong
% for i = 1:size(sigma,1)
%     div = sigma(i,1)/sigma(i,2);
%     part1 = a1*c1*c1*a2*div-2*a1*c1*a2*c2*div-a1*b1+a1*a2*c2*c2*div+a1*b2*div+b1*a2*div-a2*b2*div*div;
%     zpos1 = (-1*sqrt(part1)+a1*c1-a2*c2*div)/(a1-a2*div);
%     zpos2 = (sqrt(part1)+a1*c1-a2*c2*div)/(a1-a2*div);
%     zpos3 = (c1*c1*a2*div+b1-a2*c2*c2*div+b2*div)/(2*a2*div*(c1-c2));
%     %Choose the correct zpos - the one that is closest to 0.
%     %zpos3 is chosen if div lies below the curve - error preventing
%     if(div == a1/a2)
%         zpos(i) = zpos3;
%     elseif (abs(zpos2) < abs(zpos1))
%         zpos(i) = zpos2;
%     else
%         zpos(i) = zpos1;
%         if isnan(zpos(i))
%             zpos(i)=0;
%         end
%     end
% end
% zpos = abs(zpos);
end