function [imDPC,idx,imTop,imBot] = applyDPC_Gen(imSet,iNA,rNA,thStep)
%%
%applyDPC_Gen function
%
%Purpose: This function generates DPC images for an arbitrary number of
%axes based on the theta step size given. This value should be in radians
%
%Inputs: imSet - MxNxP 3D matrix containing images for generating DPC
%                images.
%        angles - Px2 matrix containing the illumination angle for each
%                 image in third dimension.
%        rNA - 2x1 vector containing minimum and maximum acceptable
%              illumination NAs for calculating DPC image. Has the form:
%                   [(minimum value) (maximum Value)] ex) [0, 0.2];
%        thStep - variable describing angular position for axis of asymmetry when generating DPC images.
%                 value must be in radians for proper operation of code.
%
%Outputs: imDPC - MxN DPC image using arbitrarily defined axis of
%                 asymmetry.
%         idx - Matrix containing index values for which images in imSet
%               were used for DPC image generation.
%         imTop - MxN image containing summation of all images used when
%                 generating first half of DPC image.
%         imBot - MxN image containing summation of all images used for
%                 generating bottom half of DPC image.
%
%-------------------------------------------------------------------------%


%Calculate radial NA
iNA(:,3) = sqrt(iNA(:,1).^2 + iNA(:,2).^2);

%Calculate angle for illuminations
ang = round(atan2(iNA(:,1),iNA(:,2)),2);

% Calculate average angular stepsize for annular illumination grid
for q = 1:length(ang)-1
    dtheta(q) = ang(q+1) - ang(q);
end
dtheta = abs(median(dtheta));

%Number of DPC image sets to generate
nDPC = round(pi / thStep);

%Preallocate DPC image stack
imDPC = zeros(size(imSet,1),size(imSet,2),nDPC);

%Generate image stacks for saving top image and bottom image for generating
%DPC images
imTop = imDPC;
imBot = imDPC; 

%Preallocate index matrix tracking images used for each DPC measurement
idx = zeros(size(imSet,3),2*nDPC);
cnt = 1;
%Iterate through number of DPC images to obtain
for k = 1:size(imSet,3)

    %Check what illumination NAs to use
    if(iNA(k,3) >= rNA(1) && iNA(k,3) <= rNA(2))
        %Check if angle is in fourth quadrant
        for q = 1:nDPC
            %Checks if image uses illumination from upper or lower quadrant
            if(sign(ang(k)) == 1)
                val = round((q-1)*thStep,2);
                if((ang(k)-val) > dtheta/2 && (ang(k)-val) < (pi - dtheta/2))
%                     if(q == 1 && ang(k) >= pi - dtheta/2 && ang(k) <= pi + dtheta/2)
%                     else
                        idx(k,2*q-1) = k;
                        imTop(:,:,q) = imTop(:,:,q) + imSet(:,:,k);
%                     end
                elseif((ang(k)-val) < -dtheta/2)
                    idx(k,2*q) = k;
                    imBot(:,:,q) = imBot(:,:,q) + imSet(:,:,k);                   
                end
            elseif(sign(ang(k)) == -1)
                val = round(pi - thStep*(q-1),2);
                if((abs(ang(k)) - val) > dtheta/2 && (ang(k)-val) < (pi - dtheta/2))
                    idx(k,2*q-1) = k;
                    imTop(:,:,q) = imTop(:,:,q) + imSet(:,:,k);                    
                elseif((abs(ang(k)) - val) < -dtheta/2)
                    idx(k,2*q) = k;
                    imBot(:,:,q) = imBot(:,:,q) + imSet(:,:,k);  
                end
            elseif(sign(ang(k)) == 0)
                val = round(3.14 - thStep*(q-1),2);
                %Ignores any zero-value angles for T-B case
                if((q-1)*thStep ~= 0)
                    if(k ~= (size(imSet,3)+1)/2)
                        if(val > 3.14)
                            idx(k,2*q-1) = k;
                            imTop(:,:,q) = imTop(:,:,q) + imSet(:,:,k);     
                        elseif(val < 3.14)
                            idx(k,2*q) = k;
                            imBot(:,:,q) = imBot(:,:,q) + imSet(:,:,k);
                        end
                    end
                end
                
            end
        end
        
        
        
    end %End of illumination NA selection loop
end %End of imSet for loop

%Generate Total Brightfield for each DPC image
totBF = imTop + imBot;

%Generate DPC images
imDPC = (imTop - imBot)./totBF;


end %End of 