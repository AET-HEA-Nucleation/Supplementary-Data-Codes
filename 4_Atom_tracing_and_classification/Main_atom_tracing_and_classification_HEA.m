%% Main Atom Tracing and classification
% code to perform polynomial tracing on reconstruction volume
% reconstruction volume should be upsampled by 3*3*3 linear interpolation
% each local maximum is fitted with 9*9*9 voxels (3*3*3 before interpolation)
% by 4th order polynomial equation to get the position
% Then, Non-atoms, typically originating from weak noise in the reconstruction, 
% were removed by K-means clustering of the integrated intensities around
% potential atoms. The final species are determined by applying local
% classification to account for potential intensity fluctuations in the 3D
% volume.
% Note: In the maintext, the atomic model were manually checked against the reconstruction slice 
% by slice to correct for misidentified or unidentified atoms, which typically 
% constituted a small fraction (<1%) of the total atom population in each
% HEA nanoparticle. This may result in differences in the output of this
% code.
addpath('./src/');

% Select the index of HEA dataset to proceed
ii = 4;

% load reconstructed 3D volume
Recon_filename=['../3_Final_reconstruction_volume/HEA_' num2str(ii) '_particle_volume.mat'];
Dsetvol = importdata(Recon_filename);

% Th: intensity threshold for the local maxima pixel
% local maxima with intensity less than this value will not be traced
% because they are way too weak to become actual atoms
MaxIter = 14;   CritIter = 7;   Th = 1;

% Set the experiemental pixel size for each dataset
if ii==11 || ii==15
    Res = 0.469/3;
else
    Res = 0.347/3;
end

% numpeak: maxmium number of peaks to trace
ratio=(Res/0.347*3*size(Dsetvol,1)/260)^3;  numpeak=100000*ratio;
saveInterval = 1000;    minDist = 2.0 / Res;    SearchRad = 3;
  
ourputstring = ['./Output/Atom_tracing_all_peaks_HEA_' num2str(ii) '_nanoparticle.mat'];

BoxSize0=3; %box size used for average when sorting peaks
BoxSize1=9; %box size used to find maxima
BoxSize2=9; %box size used box for fitting off-center gauss
BoxSize3=7; %used to compute visualization matrix

% upsampling the reconstruction matrix by 3*3*3 by linear interpolation
% better to run at super conputer since the size of the interpolated
% volume will be larger than 16G
xx = (1:size(Dsetvol,1)) - round((size(Dsetvol,1)+1)/2);
yy = (1:size(Dsetvol,2)) - round((size(Dsetvol,2)+1)/2);
zz = (1:size(Dsetvol,3)) - round((size(Dsetvol,3)+1)/2);

xxi = ((3*xx(1)):(xx(end)*3))/3;
yyi = ((3*yy(1)):(yy(end)*3))/3;
zzi = ((3*zz(1)):(zz(end)*3))/3;

xxi = xxi(3:end);
yyi = yyi(3:end);
zzi = zzi(3:end);

[Y,X,Z] = meshgrid(yy,xx,zz);
[Yi,Xi,Zi] = meshgrid(yyi,xxi,zzi);

Dsetvol = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'spline',0);
clear Xi Yi Zi
FinalVol = My_paddzero(Dsetvol,size(Dsetvol)+20);
FinalVol_single = single(FinalVol);

% get polynomial power array
fitCoeff = [];
for i=0:4
    for j=0:4
        for k=0:4
            if i+j+k <= 4                
                if max([i j k]) == 4
                    fitCoeff(end+1,:) = [i j k -1];                
                else
                    fitCoeff(end+1,:) = [i j k 0];                
                end
            end
        end
    end
end

% get the local maxima from the reconstruction volume
se = strel3d(3);
dilatedBW = imdilate(FinalVol,se);
maxPos = find(FinalVol==dilatedBW & FinalVol>Th);
maxVals = FinalVol(maxPos);
[~,sortInd] = sort(maxVals,'descend');
maxNum = min(numpeak,length(sortInd));
maxPos = maxPos(sortInd(1:maxNum));

fprintf(1,'numpeak = %d \n',length(maxPos));

maxXYZ = zeros(length(maxPos),3);
for i=1:length(maxPos)
    [xx,yy,zz] = ind2sub(size(FinalVol),maxPos(i));
    maxXYZ(i,:) = [xx yy zz];  
end
clear Dsetvol dilatedBW

% initialize the parameters
Q = 0.5;    Alpha = 1;
cropHalfSize = SearchRad;
[X,Y,Z] = ndgrid(-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize,-cropHalfSize:cropHalfSize);
SphereInd = find(X.^2+Y.^2+Z.^2 <=(SearchRad+0.5)^2);
XYZdata.X = X(SphereInd);
XYZdata.Y = Y(SphereInd);
XYZdata.Z = Z(SphereInd);

Orders = fitCoeff(:,1:3);
PosArr = zeros(size(maxXYZ));
TotPosArr = zeros(size(maxXYZ));

exitFlagArr = zeros(1, size(maxXYZ,1));
CoeffArr = repmat(fitCoeff(:,4),[1 size(maxXYZ,1)]);

% perform the main tracing loop
for i=1:size(maxXYZ,1)
    endFlag = 0;
    consecAccum = 0;
    iterNum = 0;
    while ~endFlag    
        iterNum = iterNum + 1;
        if iterNum>MaxIter
          exitFlagArr(i) = -4;
          endFlag = 1;
        end
        cropXind = maxXYZ(i,1) + (-cropHalfSize:cropHalfSize);
        cropYind = maxXYZ(i,2) + (-cropHalfSize:cropHalfSize);
        cropZind = maxXYZ(i,3) + (-cropHalfSize:cropHalfSize);

        cropVol = FinalVol(cropXind,cropYind,cropZind);

        Pos = PosArr(i,:);
        GaussWeight = exp(-1*Alpha*( (X(SphereInd)-Pos(1)).^2 + (Y(SphereInd)-Pos(2)).^2 + (Z(SphereInd)-Pos(3)).^2 ) / cropHalfSize^2 );
        
        fun = @(p,xdata) calculate_3D_polynomial_Rogers(xdata.X,xdata.Y,xdata.Z,Pos,Orders,p).*GaussWeight;

        opts = optimset('Display','off');
        
        [p1,fminres1] = lsqcurvefit(fun,CoeffArr(:,i),XYZdata,cropVol(SphereInd).*GaussWeight,[],[],opts);
        CoeffArr(:,i) = p1;
        
        [dX, dY, dZ] = calc_dX_dY_dZ_Rogers(Orders,CoeffArr(:,i));
        if dX ==-100 && dY == -100 && dZ == -100
            exitFlagArr(i) = -1;
            endFlag = 1;
        else
            maxedShift = max([dX dY dZ],-1*[Q Q Q]);
            minedShift = min(maxedShift,[Q Q Q]);
            PosArr(i,:) = PosArr(i,:) + minedShift;
            if max(abs(PosArr(i,:))) > cropHalfSize
                exitFlagArr(i) = -2;
                endFlag = 1;
            elseif max(abs(minedShift)) < Q
                if consecAccum == CritIter-1
                    goodAtomTotPos = TotPosArr(1:i-1,:);
                    goodAtomTotPos = goodAtomTotPos(exitFlagArr(1:i-1)==0,:);
                    Dist = sqrt(sum((goodAtomTotPos - repmat(PosArr(i,:)+maxXYZ(i,:),[size(goodAtomTotPos,1) 1])).^2,2));
                    if min(Dist) < minDist
                      exitFlagArr(i) = -3;
                    else
                      TotPosArr(i,:) = PosArr(i,:) + maxXYZ(i,:);
                    end
                    endFlag = 1;
                else
                    consecAccum = consecAccum + 1;
                end
            else
                consecAccum = 0;
            end
        end
    end
    fprintf(1,'peak %d, flag %d \n',i,exitFlagArr(i));
    
    if mod(i,saveInterval) == 0    
        save(ourputstring,'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');
    end
end
save(ourputstring,'PosArr','TotPosArr','CoeffArr','Orders','exitFlagArr');

%% Use K-mean method to remove non atom peaks from tracing results

saveprefix=['Atom_tracing_non_atom_removed_HEA_' num2str(ii) '_nanoparticle'];

% load recontrcution volume
rec = importdata(Recon_filename);
rec = double(rec);

% load tracing results
TracingResult = importdata(ourputstring);
atom_pos = TracingResult.TotPosArr(TracingResult.exitFlagArr==0,:)';

atom_pos= (atom_pos / 3) - 2;

b1 = find(atom_pos(1,:)<1 | atom_pos(1,:)>size(rec,1)-1);
b2 = find(atom_pos(2,:)<1 | atom_pos(2,:)>size(rec,2)-1);
b3 = find(atom_pos(3,:)<1 | atom_pos(3,:)>size(rec,3)-1);
bT = union(union(b1,b2),b3);
atom_pos(:,bT) = [];

% generate tight support to mask out the volume outside the sample boundary from the
% reconstruction
tight_support = My_obtain_tight_support_ver1(rec);

curr_model0      = atom_pos;

% set parameters for the K-mean calculation
lnorm = 1;  % error metric order
interp_type = 'cubic';  % upsampling method for the reconstruction
curr_model      = curr_model0;
labels          = ones(1,size(curr_model0,2));
Num_species     = 4;    % number of different species to seperate by K-mean including non-atom
Num_atom        = size(curr_model,2);
O_Ratio = 3;    % up sampling ratio
halfSize = 7/O_Ratio;   % number of pixels for atoms
plothalfSize = 4/O_Ratio;   % number of pixels for atoms for ploting
ds = 1/O_Ratio;

% obtain global intensity histogram
[XX,YY,ZZ] = ndgrid(-halfSize: ds :halfSize, -halfSize: ds :halfSize, -halfSize: ds :halfSize);
SphereInd = find(XX.^2+YY.^2+ZZ.^2 <=(halfSize+0.5*ds)^2);

[XXp,YYp,ZZp] = ndgrid(-plothalfSize: ds :plothalfSize, -plothalfSize: ds :plothalfSize, -plothalfSize: ds :plothalfSize);
SphereInd_plot = find(XXp.^2+YYp.^2+ZZp.^2 <=(plothalfSize+0.5*ds)^2);

SPHyn = 1;
if SPHyn
    useInd      = SphereInd;
    useInd_plot = SphereInd_plot;
else
    useInd = 1:length(XX);
    useInd_plot = 1:length(XXp);
end

% generate points coordinates
YY = YY(useInd); XX = XX(useInd); ZZ = ZZ(useInd);
y_set = zeros(length(XX), Num_atom);
x_set = zeros(length(YY), Num_atom);
z_set = zeros(length(ZZ), Num_atom);

YYp = YYp(useInd_plot); XXp = XXp(useInd_plot); ZZp = ZZp(useInd_plot);
y_set_plot = zeros(length(XXp), Num_atom);
x_set_plot = zeros(length(YYp), Num_atom);
z_set_plot = zeros(length(ZZp), Num_atom);

% interpolations for points
for k=1:Num_atom
    y_set(:,k) = YY + curr_model(2,k);
    x_set(:,k) = XX + curr_model(1,k);
    z_set(:,k) = ZZ + curr_model(3,k);
    
    y_set_plot(:,k) = YYp + curr_model(2,k);
    x_set_plot(:,k) = XXp + curr_model(1,k);
    z_set_plot(:,k) = ZZp + curr_model(3,k);
end
if strcmp(interp_type,'linear')
    points = splinterp3(rec, y_set,x_set,z_set);
else
    points =    interp3(rec, y_set,x_set,z_set, interp_type);
end
points(isnan(points))=0;
points_plot = splinterp3(rec, y_set_plot, x_set_plot, z_set_plot);
points_plot(isnan(points_plot))=0;
% integrate intensity (for 3x3x3 and 5x5x5 voxels) for each traced peak
intensity_integ      = sum(points);
intensity_integ_max  = max(points);
intensity_integ_plot = sum(points_plot);

separate_part = 100;

[hist_inten,cen_integ_total]= hist(intensity_integ,Num_species+1);

ind = cen_integ_total(1)<=intensity_integ & intensity_integ<=cen_integ_total(2);
avg_atom_1=mean(points(:,ind),2);

ind = cen_integ_total(2)<=intensity_integ & intensity_integ<=cen_integ_total(3);
avg_atom_2=mean(points(:,ind),2);

ind = cen_integ_total(3)<=intensity_integ & intensity_integ<=cen_integ_total(4);
avg_atom_3=mean(points(:,ind),2);

ind = cen_integ_total(4)<=intensity_integ & intensity_integ<=cen_integ_total(5);
avg_atom_4=mean(points(:,ind),2);

% atom classification iteration of non-atom and Type 1 2 3
avg_atom = [avg_atom_1,avg_atom_2,avg_atom_3,avg_atom_4];
dist_mat = zeros(Num_species,1);

[hist_inten_plot,cen_integ_total_plot]= hist(intensity_integ_plot,separate_part);

% main K-mean iteration
iter=0;
while true
    iter=iter+1;

    y_up = round(max(hist_inten_plot)/10*12);
    intensity_integ_1_plot = intensity_integ_plot(labels == 1);
    intensity_integ_2_plot = intensity_integ_plot(labels == 2);
    intensity_integ_3_plot = intensity_integ_plot(labels == 3);
    intensity_integ_4_plot = intensity_integ_plot(labels == 4);
    figure(200+2*halfSize*O_Ratio+1)
    clf
    subplot(Num_species+1,1,1);
    hist(intensity_integ_plot,separate_part);
    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
    ylim([0 y_up]);
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    title(sprintf('OR=%d, interp=%s\nboxsize %d',O_Ratio,interp_type,2*halfSize*O_Ratio+1));
    
    subplot(Num_species+1,1,2);
    hist(intensity_integ_1_plot,cen_integ_total_plot);
    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
    ylim([0 y_up]);
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    title(sprintf('%d 1 atoms',sum(labels==1)));
    
    subplot(Num_species+1,1,3);
    hist(intensity_integ_2_plot,cen_integ_total_plot)
    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
    ylim([0 y_up]);
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    title(sprintf('%d 2 atoms',sum(labels==2)));
    
    subplot(Num_species+1,1,4);
    hist(intensity_integ_3_plot,cen_integ_total_plot)
    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
    ylim([0 y_up]);
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    title(sprintf('%d 3 atoms',sum(labels==3)));   

    subplot(Num_species+1,1,5);
    hist(intensity_integ_4_plot,cen_integ_total_plot)
    xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
    ylim([0 y_up]);
    xlabel('integrated intensity (a.u.)');
    ylabel('# atoms');
    title(sprintf('%d 3 atoms',sum(labels==4)));  

    avg_atom1 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);
    avg_atom2 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);
    avg_atom3 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);
    avg_atom4 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);
    
    avg_atom1(useInd) = avg_atom(:,1);
    avg_atom2(useInd) = avg_atom(:,2);
    avg_atom3(useInd) = avg_atom(:,3);
    avg_atom4(useInd) = avg_atom(:,4);
    
    figure(300+2*halfSize*O_Ratio+1);
    img(sum(avg_atom1,3),'non-atom', sum(avg_atom2,3),'atom1',sum(avg_atom3,3),'atom2',sum(avg_atom4,3),'atom3');
    drawnow;

    old_labels = labels;
       
    obj = 0;
    for n = 1:Num_atom
        for k=1:Num_species
            dist_mat(k) = norm(points(:,n) - avg_atom(:,k), lnorm).^lnorm/sum(points(:,n).^lnorm);
        end
        [dist2, idx] = min(dist_mat);        
        labels(n) = idx;
        obj = obj + dist2;
    end
    
    for k=1:Num_species
        avg_atom(:,k) = sum(points(:,labels==k),2)/sum(labels==k);
    end
    avg_atom(isnan(avg_atom))=0;
    avg_atom(:,1)=0; % constraint non atom peak to be zero
    fprintf('%02i. obj = %.3f\n',iter, (obj/Num_atom)^(1/lnorm) );
    
    % if there is no change in the atomic specise classification, break
    if ~any(old_labels-labels), break; end
    
end
% apply support
atomtype=labels;
for i=1:length(atomtype)  %if previous determined noatom is within the tight support, then change back to atom
        Rpos = round(curr_model(:,i));
        if tight_support(Rpos(1),Rpos(2),Rpos(3))>0
            if atomtype(i)==1
                atomtype(i) = 2; % within the support but classified as NA atom
            end
        else
            if atomtype(i)~=1
                atomtype(i) = 1; % outside the support but classified as atom
            end
        end
end
labels=atomtype;
% delete non-atom
curr_model(:,labels==1)=[];
labels(labels==1)=[];
labels = labels-1;
save(['./Output/' saveprefix],'curr_model','labels','intensity_integ_plot','cen_integ_total_plot','intensity_integ_max','hist_inten_plot');

%% perform local k-mean on the summed intensity of 9*9*9 volume around traced atom position

data=importdata(['./Output/' saveprefix '.mat']);
new_model = data.curr_model;
atom = data.labels;

Dsetvol         = importdata(Recon_filename);

xx = (1:size(Dsetvol,1)) - round((size(Dsetvol,1)+1)/2);
yy = (1:size(Dsetvol,2)) - round((size(Dsetvol,2)+1)/2);
zz = (1:size(Dsetvol,3)) - round((size(Dsetvol,3)+1)/2);
xxi = ((3*xx(1)):(xx(end)*3))/3;
yyi = ((3*yy(1)):(yy(end)*3))/3;
zzi = ((3*zz(1)):(zz(end)*3))/3;
xxi = xxi(3:end); yyi = yyi(3:end); zzi = zzi(3:end);
[Y,X,Z] = meshgrid(yy,xx,zz);
[Yi,Xi,Zi] = meshgrid(yyi,xxi,zzi);
Dsetvol = interp3(Y,X,Z,Dsetvol,Yi,Xi,Zi,'cubic',0);
clear Xi Yi Zi

FinalVol = My_paddzero(Dsetvol,size(Dsetvol)+20);
FinalVol_single = single(FinalVol);

% set parameters
% the half size of volume around atom to calculate intensity, halfsize = 3
% means to sum a 7*7*7 box around the atom to get the intensity
classify_info = struct('Num_species', 3,  'halfSize',  3,  'plothalfSize',  1, ...
      'O_Ratio', 1, 'SPHyn',  1,  'PLOT_YN',  1,  'separate_part',  70);

% calculate the atomic position with upsampled volume
new_model_L = (new_model +2).*3;

% apply k-mean classification on the reconstruction
[temp_model, temp_atomtype] = initial_class_kmean_sub(...
    FinalVol_single, new_model_L, classify_info);

classify_info.Radius = 10/Res; % radius is 10A

% when there are 5 iterations with same number of atoms flipped (back and
% forth), the interation will be stopped
[temp_model, local_atomtype] = local_class_kmean_sub(FinalVol_single, temp_model, temp_atomtype, classify_info);

% apply function plot_class_hist() to achieve the histogram information peak_info
% please see the descriptions in subfunction to get more details
[peak_info,intensity_arr,intensity_plot_arr] = ...
    plot_class_hist(FinalVol_single, temp_model, local_atomtype, classify_info);

save(['./Output/Atom_local_classified_HEA_' num2str(ii) '_nanoparticle.mat'],'curr_model','local_atomtype')