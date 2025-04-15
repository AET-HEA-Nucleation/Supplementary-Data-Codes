clear
rec = importdata('Recon_ARSP.mat');
recon=rec.reconstruction;

tracing=importdata('otsu0.9NAatom_L1_Kmean.mat');
model=tracing.curr_model;
atom=tracing.labels;

% model=temp_model-10;
% atom=ones(1,size(model,2));

dimz=size(recon,3);

P=[1 2 3];
recon=permute(recon,P);
model=model(P,:);

% save CT2view dimz recon model atom
figure(1)
% scatter3(model(1,atom==1),model(2,atom==1),model(3,atom==1),'go','markerfacecolor','g')
% hold on
scatter3(model(1,atom~=1),model(2,atom~=1),model(3,atom~=1),'ko','markerfacecolor','r')
hold off
axis image

% ind=model(1,:)>50 & model(1,:)<200 & model(2,:)>50 & model(2,:)<200 & model(2,:)>50 & model(2,:)<200;
% sum(ind)
%% View species

dz=1;
clim1=[400*dz/3 (dz+0.5)*max(recon,[],'all')];
for i=dz+1:dz:dimz-dz
    slice=sum(recon(:,:,i+(-dz:dz)),3);
%     ind1=(model(3,:)<(i+dz) & model(3,:)>=(i-dz)) & atom==1;
    ind2=(model(3,:)<(i+dz) & model(3,:)>=(i-dz)) & atom~=1;
    
    figure(2)
    imagesc(slice,clim1)
    colormap gray
    axis image
    hold on
%     scatter(model(2,ind1),model(1,ind1),'go','markerfacecolor','g')
    scatter(model(2,ind2),model(1,ind2),'ro','markerfacecolor','r')
    hold off
    pause
end
%% compare tracing
model2=importdata('../model_refined_HEA_181027_t1_yk55_yk_0829.mat')/0.347+151+[0.2677 0.1054 -0.4869]'/0.347;

dz=1;
clim1=[000*dz/3 (dz+0.5)*max(recon,[],'all')];
for i=dz+1:dz:dimz-dz
    slice=sum(recon(:,:,i+(-dz:dz)),3);
    ind1=(model2(3,:)<(i+dz) & model2(3,:)>=(i-dz));
    ind2=(model(3,:)<(i+dz) & model(3,:)>=(i-dz)) & atom~=1;
    
    figure(2)
    imagesc(slice,clim1)
    colormap gray
    axis image
    hold on
    scatter(model2(2,ind1),model2(1,ind1),'go')
    scatter(model(2,ind2),model(1,ind2),'ro')
    hold off
    pause
end