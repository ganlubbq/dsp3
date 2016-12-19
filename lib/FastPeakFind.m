function  cent=FastPeakFind(d, threshold, filt ,edg, filename)
% Analyzes noisy 2D images , finds x-y positions of peaks to 1 pixel accuracy
% The code was meant to be as fast as possible, so I kept it pretty basic.
% The code assumes that the peaks are relatively sparse, test whether there 
% is too much pile up and set threshold or user defined filter accordingly.   
%
% Inputs:
%   d           The 2D data raw image - assumes a Double\Single-precision
%               floating-point or unit16 array. Please note that the code
%               casts the raw image to uint16, if your image is 8-bit depth,
%               I recommend to change all the places that uint16 is used to
%               uint8 for faster run times. If your image dynamic range is between 0 and 1, I
%               multiplied to fit uint16. This might not be optimal for
%               generic use, so modify according to your needs.
%   threshold   A number between 0 and max(raw_image(:)) to remove  background
%   filt        A filter matrix used to smooth the image. The filter size
%               should correspond the characteristic size of the peaks 
%   edg         A number>1 for skipping the first few and the last few 'edge' pixels
%   filename    In case the user would like to save the peak positions to a file
%
%
% Output:
%   cent        a 2xN matrix of the x and y coordinates of peaks
%
%   Example:
%
%   p=FastPeakFind(image);
%   imagesc(image); hold on
%   plot(p(:,2),p(:,1),'k.'); hold off
%
%   Comments \ improvements are welcomed
%   Nate (nate2718281828@gmail.com)
%   Ver 1.3 , Date: 8/20/2012.

%% defaults
% if (nargin < 1)
%     d=conv2(reshape(single( 255*(rand(1,256*256)>0.9995) ),[256 256]) ,fspecial('gaussian', 15,2),'same')+5*rand(256);
%     imagesc(d)
% end

% if size(d,3)>1 %I added this in case one uses imread (JPG\PNG\...).
%     d=uint16(rgb2gray(d));
% end

% if isfloat(d) %For the case the input image is double, casting to uint16 keeps enough dynamic range while speeds up the code.
%    if max(d(:))<=1
%        d =  uint16( d.*2^16./(max(d(:))));
%    else
%        d = uint16(d);
%    end
% end

if (nargin < 2)
    threshold = (max([min(max(d,[],1))  min(max(d,[],2))])) ;
end

if (nargin < 3)
    filt = (fspecial('gaussian', 7,1)); %if needed modify the filter according to the expected peaks sizes
end

if (nargin < 4)
    edg =3;
end

if (nargin < 5)
    savefileflag = false;
else
    savefileflag = true;
end

%%
if any(d(:))  ; %for the case of non zero raw image
%     d = medfilt2(d,[3,3]);
%     
%     % apply threshold
%     d=d.*uint16(d>threshold);
    
    if any(d(:))  ; %for the case of the image is still non zero
        
%         % smooth image
%         d=conv2(single(d),filt,'same') ;
        
        % Apply again threshold (and change if needed according to SNR)
        d=d.*(d>0.9*threshold);
        
        
        % peak find - using the local maxima approach - 1 pixle resolution
        % d will be noisy on the edges, since no hits are expected there anyway we'll skip 'edge' pixels.
        
        [x y]=find(d(edg:size(d,1)-edg,edg:size(d,2)-edg));
        
        % initialize peak find outputs
        cent=[[];[]];
        x=x+edg;
        y=y+edg;
        for j=1:length(y)
            if (d(x(j),y(j))>=d(x(j)-1,y(j)-1 )) &&...
                    (d(x(j),y(j))>d(x(j)-1,y(j))) &&...
                    (d(x(j),y(j))>=d(x(j)-1,y(j)+1)) &&...
                    (d(x(j),y(j))>d(x(j),y(j)-1)) && ...
                    (d(x(j),y(j))>d(x(j),y(j)+1)) && ...
                    (d(x(j),y(j))>=d(x(j)+1,y(j)-1)) && ...
                    (d(x(j),y(j))>d(x(j)+1,y(j))) && ...
                    (d(x(j),y(j))>=d(x(j)+1,y(j)+1));

%All these alternatives were slower...
%if all(reshape( d(x(j),y(j))>=d(x(j)-1:x(j)+1,y(j)-1:y(j)+1),9,1))
%if  d(x(j),y(j)) == max(max(d((x(j)-1):(x(j)+1),(y(j)-1):(y(j)+1))))
%if  d(x(j),y(j))  == max(reshape(d(x(j),y(j))  >=  d(x(j)-1:x(j)+1,y(j)-1:y(j)+1),9,1)) 
 
 
               cent= [cent ; [x(j),y(j)]];
                
            end
        end
        
        if savefileflag
            dlmwrite([filename '.txt'],[cent],   '-append', ...
                'roffset', 0,   'delimiter', '\t', 'newline', 'unix');
        end
        
    else % in case image after threshold is all zeros
        cent=[];
        return
    end
    
else % in case raw image is all zeros (dead event)
    cent=[];
    return
end

if (nargin < 1); colormap(bone);hold on; plot(cent(:,2),cent(:,1),'r+');hold off; end

return
