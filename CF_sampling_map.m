function [ re_sum_CF_w ] = CF_sampling_map( CF_size, vdist, CF_ecc, CF_polarangle, index_cfconnection, radius, pRF_ecc, pRF_pol, pRF_size)
% CF_coverage creates a map of the cortical sampling based on the
% connective fields between a source and target ROI and the pRF estimates.
% This function requires the following inputr
% CF_size=connective field size estimates for all the voxels
% Vdist= distance in cortical surface between the target and the source for all the voxels
% CF_ecc= eccentricity estimate obtained with connective field model for all the voxels
% CF_polarangle= polar angle (phase) estimate with connective field model for all the voxels
% index_cfconnection=indext of the connection voxel source-target estimate with connective field model for all the voxels
% radius= radius of visual stimulation
% pRF_ecc= pRF eccentricity of the source ROI
% pRF_pol=pRF polar angle of the source ROI
% pRF_size=pRF size of the source ROI

% 1. reconstruction of the connective field
for n_voxel=1:size(CF_size,2)
    for u=1:size(index_cfconnection,2)
        n=index_cfconnection(n_voxel);
        CF(:,n_voxel)=single(exp(-1.*((vdist(:,n).^2)./(2.* CF_size(n_voxel).^2))));
    end
end
                      
% 2. reconstruction of the pRF sampling map

size_map=50;
n = 1;
matSize = 51;

[x, y] = meshgrid(0:radius/size_map:radius);
x=x(:);y=y(:);

for n_voxel=1:size(pRF_size,2)
    x0=abs(pRF_ecc(n_voxel).*cos(pRF_pol(n_voxel)));
    y0=abs(pRF_ecc(n_voxel).*sin(pRF_pol(n_voxel)));
    pRF(:,n_voxel) = exp(((x-x0).^2+(y-y0).^2)/(-2*pRF_size(n_voxel).^2));
end

sumV1=sum(pRF,2);

re_sumV1= reshape(sumV1,matSize,[],n);

% figure,
% imagesc(flipud(re_sumV1)) %Vizualize pRF sampling map


% 3. Calculate CF sampling map

CF_w=pRF*CF;
sum_CF_w=sum(CF_w,2);
re_sum_CF_w=reshape(sum_CF_w,matSize,[],n);
re_sum_CF_w=flipud(re_sum_CF_w);

%figure, imagesc(re_sum_CF_w)  %Vizualize CF sampling map

end