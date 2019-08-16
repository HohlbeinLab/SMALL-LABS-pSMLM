%%Loading tiff files fast
function [FinalImage] = TiffLoader(FileTif,varargin)
warning('off')
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
if nargin == 1 %if full movie
    NumberImages=length(InfoImage);
    startframe = 1;
    endframe = NumberImages;
else
    framearr = varargin(1);
    startframe = framearr{1}(1);
    endframe = framearr{1}(2);
end
% FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
FinalImage=zeros(nImage,mImage,endframe-startframe+1,'uint16');
 
TifLink = Tiff(FileTif, 'r');
% for i=1:NumberImages
for i=startframe:endframe
   TifLink.setDirectory(i);
   FinalImage(:,:,i-startframe+1)=TifLink.read();
end
TifLink.close();
warning('on')
end