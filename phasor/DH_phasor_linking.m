%% Testing method(s) for difficult DH-linking
function [particles, allpossibilities] = DH_phasor_linking(loc_list,curveAngleinput,linking_range,magnratiorange,wobbleMatrix)
%%
% Assumed input: loc_list, consisting of frame, x, y, fwhmx, fwhmy columns
%   linking_range: [min_range max_range] OR a matrix 
%   curveAngle: NonLinearModel matching the angle to z-pos. Output of
%   pSMLM_DH_Calibration. If left empty ([]), no zpos will be calculated.
%   magnitude ratio range: similar layout as wobbleMatrix - information on
%   phasor magnitude ratio that are acceptable. However, I find that this
%   is never worth using (unsure why). If left empty, it's unused.
%   wobbleMatrix: matrix in form of [minzpos maxzpos wobbleX wobbleY] to
%   correct for wobble effects. Is an output of pSMLM_DH_calibration. Can
%   be left empty for no wobble correction.
% Output: particles array
%     columns: frame, Xpos, Ypos, Zpos, rotation, distance
%%
% Update list
% 2019-03-26: Added wobble correction
% 2019-03-25: first version
%%
%set min and max range for lookup of combinations - this should be based on
%the calibration curve, and maybe even dependant on the angle with the
%other spot as well.
%If it is a 1-by-2 matrix, it just has min and max range.
%If it's anything else, it's a calibration curve.
global verbose
if verbose
disp([char(datetime),'   Starting DH phasor'])
end

if ((size(linking_range,1) == 1) && (size(linking_range,2) == 2))
    min_range = linking_range(1);
    max_range = linking_range(2);
    nonlinearmodel = 0;
else %It's a NonLinearModel, which has some great information for us!
    %estimate min and max on min and max used values in calibration, with
    %10% increased size - 10% works better than 0% or 25%.
    min_range = min(linking_range.Variables.y)*.90;%*0.9;
    max_range = max(linking_range.Variables.y)*1.1;%*1.1;
    nonlinearmodel = 1;
end
guess_possible_linking = 0; %0 if only deterministic, 1 if we want to guess unresolvable linkages
abs_error_from_estimated_distance = 1; %was 1; %maximum difference (pixels) between found distance of 2 spots and estimated distance of 2 spots based on their rotation. Unused if just a 1-by-2 matrix is used for min/max linking distance.

particles = zeros(floor(size(loc_list,1)/2),9); %needs to be zeros matrix
curparticle_id = 1;
%% Start analysis
%loop over all frames

%Get a,b,c,d from linking_range if they exist
if ~((size(linking_range,1) == 1) && (size(linking_range,2) == 2))
    t = table2array(linking_range.Coefficients(:,1));
    lr_a = t(1); lr_b = t(2); lr_c = t(3); lr_d = t(4);
end

%Get info from curve-angle calibration if it exists
if ~isempty(curveAngleinput)
    %zpos calculation
    %Readability of coefficients
    a = curveAngleinput.Coefficients{1,1};
    b = curveAngleinput.Coefficients{2,1};
    c = curveAngleinput.Coefficients{3,1};
    d = curveAngleinput.Coefficients{4,1};
end

%Quick frame index list lookup
frame_index_list = accumarray([loc_list(:,1)], (1:size(loc_list,1))', [], @(idxlist) {idxlist});


for frame = 1:max(loc_list(:,1))
    %Get localizations on this frame
%     tt = loc_list(:,1)==frame;
%     loc_array = loc_list(tt,2:5);%loc_array = loc_list(tt,2:3);
%     ts = find(loc_list(:,1)==frame);
%     loc_array = loc_list(ts,2:5);
    t = frame_index_list(frame); tt = t{1};
    loc_array = loc_list(tt,2:5);
%     if ts ~= tt
%         fprintf('Hm.\n');
%     end
    %Add id
    loc_array = [[1:size(loc_array,1)]' loc_array];
    %Calculate all euclidian distances between all points
    if size(loc_array,1)>35 %If its incredibly large, pdist is faster
        EuclidianDistances = squareform(pdist(loc_array(:,2:3)));
    else %else 2 forloops is faster for some reason...
        EuclidianDistances = [];
        for ii = 1:size(loc_array,1)
            for j = 1:size(loc_array,1)
                EuclidianDistances(ii,j) = sqrt((loc_array(ii,2)-loc_array(j,2))^2+(loc_array(ii,3)-loc_array(j,3))^2);
            end
        end
    end
    %Calculate which distances are possible to have a matching pair, currently
    %depending on min_range and max_range
    possibilities = cell(size(loc_array,1),4);
    Possibilities_size = [];
    if size(loc_array,1) > 1
        for ii = 1:size(loc_array,1)
            possibilities{ii,1} = ii;
            for j = 1:size(loc_array,1)
                if EuclidianDistances(ii,j) > min_range
                    if EuclidianDistances(ii,j) < max_range
                        %If linking_range is a 1-by-2 matrix, it just has min and max range.
                        %If it's anything else, it's a calibration curve.
                        %If it just has min and max values, use only this
                        if ~nonlinearmodel
                            %Then check on phasor magnitude ratio
                            if ~isempty(magnratiorange)
                                %get correct ratio (top div by bottom)
                                if loc_array(ii,3) > loc_array(j,3)
                                    magnratiotopdivbot = (loc_array(ii,4)/loc_array(ii,5))/(loc_array(j,4)/loc_array(j,5));
                                else
                                    magnratiotopdivbot = (loc_array(j,4)/loc_array(j,5))/(loc_array(ii,4)/loc_array(ii,5));
                                end
                                %calc zpos
                                tempangle = atan2((loc_array(ii,3)-loc_array(j,3)),(loc_array(ii,2)-loc_array(j,2)));
                                y = tempangle;
                                zpos = real(-((1 + i.*sqrt(3)).*(sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 ...
                                    + 4.*(3.*a.*c - b.^2).^3) - 27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3))./(6.*2.^(1./3).*a) ...
                                    + ((1 - i.*sqrt(3)).*(3.*a.*c - b.^2))./(3.*2.^(2./3).*a.* ...
                                    (sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) ...
                                    - 27.*a.^2.*d + 27.*a.^2..*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)) - b./(3.*a));
                                %see if magn ratio is within bounds of
                                %zpos:
                                %find zpos bin
                                zposBin = max(find(zpos>magnratiorange(:,1)));
                                if (magnratiotopdivbot >= magnratiorange(zposBin,3) & magnratiotopdivbot <= magnratiorange(zposBin,4))
                                    possibilities{ii,2} = [possibilities{ii,2} j];
                                else
                                    fprintf('Excluded based on magnitude ratio frame %.0f\n',frame)
                                end
                            else %if no info about magnitude ratios is present
                                possibilities{ii,2} = [possibilities{ii,2} j];
                            end
                        else %Else, we have more info
                            %Calculate the angle
                            tempangle = atan2((loc_array(ii,3)-loc_array(j,3)),(loc_array(ii,2)-loc_array(j,2)));
                            %Estimate distance of DH on this angle
    %                         estimatedDist = predict(linking_range,tempangle);
                            estimatedDist = (lr_a*tempangle^3+lr_b*tempangle^2+lr_c*tempangle+lr_d);

                            %check if error on distance isn't too large
                            if abs(estimatedDist-EuclidianDistances(ii,j))<=abs_error_from_estimated_distance
                                %Then check on phasor magnitude ratio
                                if ~isempty(magnratiorange)
                                    %get correct ratio (top div by bottom)
                                    if loc_array(ii,3) > loc_array(j,3)
                                        magnratiotopdivbot = (loc_array(ii,4)/loc_array(ii,5))/(loc_array(j,4)/loc_array(j,5));
                                    else
                                        magnratiotopdivbot = (loc_array(j,4)/loc_array(j,5))/(loc_array(ii,4)/loc_array(ii,5));
                                    end
                                    %calc zpos
                                    y = tempangle;
                                    zpos = real(-((1 + i.*sqrt(3)).*(sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 ...
                                        + 4.*(3.*a.*c - b.^2).^3) - 27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3))./(6.*2.^(1./3).*a) ...
                                        + ((1 - i.*sqrt(3)).*(3.*a.*c - b.^2))./(3.*2.^(2./3).*a.* ...
                                        (sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) ...
                                        - 27.*a.^2.*d + 27.*a.^2..*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)) - b./(3.*a));
                                    %see if magn ratio is within bounds of
                                    %zpos:
                                    %find zpos bin
                                    zposBin = max(find(zpos>magnratiorange(:,1)));
                                    if (magnratiotopdivbot >= magnratiorange(zposBin,3) & magnratiotopdivbot <= magnratiorange(zposBin,4))
                                        possibilities{ii,2} = [possibilities{ii,2} j];
                                    else
                                        fprintf('Excluded based on magnitude ratio frame %.0f\n',frame)
                                    end
                                else %if no info about magnitude ratios is present
                                    possibilities{ii,2} = [possibilities{ii,2} j];
                                end
                            end
                        end
                    end
                end
            end
            possibilities{ii,3} = size(possibilities{ii,2},2);
            %set link to 0 for starting
            possibilities{ii,4} = 0;
        end
    end
    %Find single possibilities if they exist
    %First, sort possibilities cell on the size of the possibilities
%     possibilities = sortrows(possibilities,[3 1]);

    %now, have a while loop that continues untill all of the spots have only 1
    %link left
    % keyboard
    stuck = 0;
    while (stuck == 0) %% sum([possibilities{:,3}]) > size(possibilities,1) && 
        %Loop for all the possibilities that are not yet linked
        stuck = 0;
        poss_old = possibilities;
        for ii = find([possibilities{:,3}] > 0)
            %if there's only 1 option which is not yet assigned
            if (possibilities{ii,3} == 1 && possibilities{ii,4} == 0)
                %Set the possibility of the other to this option as well
                possibilities{[possibilities{:,1}] == [possibilities{ii,2}],2} = possibilities{ii,1};
                %Set the linking of this and the other to each other
                %Of the one being checked
                possibilities{ii,4} = possibilities{ii,2};
                %and of its partner
                possibilities{[possibilities{:,1}] == [possibilities{ii,2}],4} = possibilities{ii,1};

                %remove the same option from others
                %For this, check the other arrays
                currPartner = find([possibilities{:,1}] == [possibilities{ii,2}]);
                for k = 1:size(possibilities,1)
                    %check if k is not the current instance
                    if k ~= ii
                        %or its partner
                        if k ~= currPartner
                            %skip this is its empty
                            if ~isempty([possibilities{k,2}])
                                %now, only all others are found
                                %remove index of i
                                possibilities{k,2}(find([possibilities{k,2}] == possibilities{ii,1})) = [];
                                %remove index of partner
                                possibilities{k,2}(find([possibilities{k,2}] == possibilities{ii,2})) = [];
                            end
                        end
                    end
                end

                %Sort possibilities again
                for k = 1:size(possibilities,1)
                    possibilities{k,3} = size(possibilities{k,2},2);
                end
% I've had this line, but removal seems to not harm anything, and speed up code.
%               possibilities = sortrows(possibilities,[3 1]);
                ii = size(possibilities,1)+1; %end the parent loop, go back to the while-loop
            end
        end
        %If nothing changes in a loop, just assume that all other particles are
        %unresolvable.
        if isequal(possibilities, poss_old)
            %if we don't want to guess anything
            if guess_possible_linking == 0
                %Collapse all multiple option possibilities to no option at all
                %(ignored basically)
                for k = 1:size(possibilities,1)
                    if size(possibilities{k,2},2)>1
                        possibilities{k,2} = [];
                    end
                end
                %Get out of the while-loop
                stuck = 1;
            %otherwise, if we want to choose the option that's closest to
            %centre of guesstimate or something...
            else
                %todo
            end
        end
    end
%     possibilities = sortrows(possibilities,[1]);
    
    %% Create final x,y positions
    %First, find only 1 original id of linked particles
    %For this, check if the id of the linked particle is smaller than the
    %id of its linked particle. Only keep those ids.
    idlist = [];
    %Loop over all input particles
    for k = 1:size(possibilities,1)
        %check if it's linked
        if ~isempty(possibilities{k,2})
            %check if it's the smallest
            if possibilities{k,1} < possibilities{k,2}
                idlist = [idlist; [k possibilities{k,2}]];
            end
        end
    end
    %Now, create the particles table with x,y,rotation,distance
    particlesT = zeros(size(idlist,1),6);
    for k = 1:size(idlist,1)
        %xypos of linked particles
        xpos1 = loc_array(idlist(k,1),2);
        xpos2 = loc_array(idlist(k,2),2);
        ypos1 = loc_array(idlist(k,1),3);
        ypos2 = loc_array(idlist(k,2),3);
        %magnitude ratios of linked particles
        if ypos1 > ypos2
            magnratio1 = loc_array(idlist(k,1),4)/loc_array(idlist(k,1),5); % top half
            magnratio2 = loc_array(idlist(k,2),4)/loc_array(idlist(k,2),5); % bottom half
        else
            magnratio1 = loc_array(idlist(k,2),4)/loc_array(idlist(k,2),5); % top half
            magnratio2 = loc_array(idlist(k,1),4)/loc_array(idlist(k,1),5); % bottom half
        end
        %other calculations
        particlesT(k,1) = frame;
        particlesT(k,2) = xpos1/2+xpos2/2; %x pos
        particlesT(k,3) = ypos1/2+ypos2/2; %y pos
        particlesT(k,4) = 0; %z pos
        particlesT(k,5) = mod(atan2((ypos2-ypos1),(xpos2-xpos1)),pi()); %angle = atan2(dy,dx)
        particlesT(k,6) = sqrt((xpos2-xpos1)^2+(ypos2-ypos1)^2);
        particlesT(k,7) = magnratio1;
        particlesT(k,8) = magnratio2;
        particlesT(k,9) = magnratio1/magnratio2;
    end
    particles(curparticle_id:curparticle_id-1+size(particlesT,1),:) = particlesT;
    curparticle_id = curparticle_id+size(particlesT,1);
    
    allpossibilities{frame,:,:} = possibilities;
end
particles(particles(:,1)==0,:) = [];

%% zpos calculation
if ~isempty(curveAngleinput)
    %zpos calculation
    %Readability of coefficients
    a = curveAngleinput.Coefficients{1,1};
    b = curveAngleinput.Coefficients{2,1};
    c = curveAngleinput.Coefficients{3,1};
    d = curveAngleinput.Coefficients{4,1};
    y = particles(:,5); %angle array
    %zpos can be calculated from the following monstrous equation (x^3...)
    particles(:,4) = real(-((1 + i.*sqrt(3)).*(sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 ...
        + 4.*(3.*a.*c - b.^2).^3) - 27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^(1./3))./(6.*2.^(1./3).*a) ...
        + ((1 - i.*sqrt(3)).*(3.*a.*c - b.^2))./(3.*2.^(2./3).*a.* ...
        (sqrt((-27.*a.^2.*d + 27.*a.^2.*y + 9.*a.*b.*c - 2.*b.^3).^2 + 4.*(3.*a.*c - b.^2).^3) ...
        - 27.*a.^2.*d + 27.*a.^2..*y + 9.*a.*b.*c - 2.*b.^3).^(1./3)) - b./(3.*a));
end
%% Wobble correction
%Currently, slow-ish: checks every particle again
%loop over all particles
if ~isempty(wobbleMatrix)
    for ii = 1:size(particles,1)
        wobbleBin = max(find(particles(ii,4)>wobbleMatrix(:,1)));
        if isempty(wobbleBin) %if empty, it's before the first bin of the wobbleMatrix
            %We'll just use the first wobbleMatrix entry
            particles(ii,2:3) = particles(ii,2:3)-wobbleMatrix(1,3:4);
        else
            %Correct for the correct wobbleMatrix bin
            particles(ii,2:3) = particles(ii,2:3)-wobbleMatrix(wobbleBin,3:4);
        end
    end
end
%% end function
end