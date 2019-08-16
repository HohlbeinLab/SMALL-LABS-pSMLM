%% Testing method(s) for SaddlePoint-TetraPod-linking
function [finalarr] = SP_TP_phasor_linking(loc_list,maxdist)
%%
% Assumed input: loc_list, consisting of frame, x, y columns
%   
% Output: particles array
%     columns: frame, Xpos, Ypos, Zpos
%%
% Variables
maxdist_secondaryDir = 1;
maxDistMultiplier = 1.5;
%%
global verbose
if verbose
disp([char(datetime),'   Starting SP/TP phasor'])
end
% keyboard
%% Loop over all frames
    finalarr = [];
for curframe = 1:max(loc_list(:,1))
    %Select localizations in current frame
    tarr = loc_list(loc_list(:,1)==curframe,:);
    %add ones (availability for linking) to all locs
    framearr = [tarr ones(size(tarr,1),1)];
%     framearr = [1 5 5.11 1; 1 12 5.2 1; 1 22 5.1 1;];%Testing stuck
%     framearr = [1 5 5.11 1; 1 12 5.2 1; 1 22 5.1 1; 1 22 15.1 1]; %shouldnt get stuck
    %Pair all localizations in the frame
    clear DistanceX DistanceY validcombinations availablearr
    for ii = 1:size(framearr,1)
        for j = 1:size(framearr,1)
%             if ii > j
                DistanceX(ii,j) = abs(framearr(ii,2)-framearr(j,2));
                DistanceY(ii,j) = abs(framearr(ii,3)-framearr(j,3));
                availablearr(ii,j) = 1;
                if ii == j
                    availablearr(ii,j) = 1;
                end
%             end
        end
    end
    allassigned = 0;
    finalmidpointsforctPhasor = [];
    stuck = false;
    while ((allassigned == 0) && (size(framearr,1)>0))
        %Get combinations that are valid
        validcombinations = ((DistanceX<maxdist_secondaryDir).*(DistanceY<maxdist.*maxDistMultiplier)+...
            (DistanceY<maxdist_secondaryDir).*(DistanceX<maxdist.*maxDistMultiplier)).*availablearr;
        for ii = 1:size(validcombinations,1)
            validcombinations(ii,ii) = 0;
        end
        %loop over columns
        availablearrstart = availablearr;
        for combcol = 1:size(validcombinations,2)
            if sum(validcombinations(:,combcol)) == 1
                % Only one combination, mix them, and set them as
                % unavailable for next iteration
                correspondingrow = find(validcombinations(:,combcol)==1);
                %Add the middle point between them for phasor
                finalmidpointsforctPhasor = [finalmidpointsforctPhasor;...
                    curframe... 
                    framearr(combcol,2)*.5+framearr(correspondingrow,2)*.5 ...
                    framearr(combcol,3)*.5+framearr(correspondingrow,3)*.5 ...
                    ];
                %Set them as unavailable
                availablearr(:,combcol) = 0;
                availablearr(:,correspondingrow) = 0;
                availablearr(combcol,:) = 0;
                availablearr(correspondingrow,:) = 0;
                break
            elseif sum(validcombinations(:,combcol)) == 0
                %if it has no combinations, check if its not already
                %assigned to another loc
                if sum(availablearr(:,combcol))>0
                    %Add it for phasor
                    finalmidpointsforctPhasor = [finalmidpointsforctPhasor;...
                        curframe framearr(combcol,2) framearr(combcol,3)];
                    %Set it as unavailable
                    availablearr(:,combcol) = 0;
                    availablearr(combcol,:) = 0;
                end
            else
                %only continue if it's stuck for a full loop
                if stuck
                    %get the one with the lowest displacement in the
                    %low-displacement direction
                    mininX = min(min(DistanceX(validcombinations>0)));
                    mininY = min(min(DistanceY(validcombinations>0)));
                    %link the ones with the lowest displacement
                    if mininX < mininY
                        [r,c] = ind2sub(size(DistanceX),min(find(DistanceX==mininX)));
                    else
                        [r,c] = ind2sub(size(DistanceY),min(find(DistanceY==mininY)));
                    end
                    %Add the middle point between them for phasor
                    finalmidpointsforctPhasor = [finalmidpointsforctPhasor;...
                        curframe... 
                        framearr(r,2)*.5+framearr(c,2)*.5 ...
                        framearr(r,3)*.5+framearr(c,3)*.5 ...
                        ];
                    %Set them as unavailable
                    availablearr(:,r) = 0;
                    availablearr(:,c) = 0;
                    availablearr(r,:) = 0;
                    availablearr(c,:) = 0;
                    break
                end
            end
            if sum(sum(availablearr))==0
                allassigned = 1;
            end
        end
        if isequal(availablearr,availablearrstart)
            stuck = true;
        else
            stuck = false;
        end
    end
    clear tarr framearr
    finalarr = [finalarr; finalmidpointsforctPhasor];
end
%% end function
end