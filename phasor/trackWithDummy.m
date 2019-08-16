function tracks = trackWithDummy(pos,trackParams)
% add a dummy track outside the field of view to ensure continuous spacing in time
dummyTrack(:,1:2) = 10000*ones(max(pos(:,3)),2);
dummyTrack(:,3) = 1:max(pos(:,3)); %continuous frame counter

posAndDummy = [pos; dummyTrack];

[B,IX] = sort(posAndDummy(:,3),1); % sort tracks by ascending frame number
posAndDummy = posAndDummy(IX,1:3); % positions without tracking ordered by frame number

tracks = track(posAndDummy,trackParams.maxDisp,trackParams);
% delete dummy track
tracks(tracks(:,1) == 10000,:) = [];

end