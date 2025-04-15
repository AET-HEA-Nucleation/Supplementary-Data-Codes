%% Main Position refinement 
% refinement of the atomic coordinates by minimizing the error between the
% calculated linear projections and the measured projections
clear
clc
addpath('./src/')

% Select the index of HEA dataset to proceed
ii= 4;

% save filename prefix
savename=['HEA_' num2str(ii) '_'];

% experimental pixle size
if ii==11 || ii==15
    Res=0.469;
else
    Res=0.347;
end

% read in files: measured projections, angles and atomic position
projections=importdata(['../1_Measured_data/HEA_' num2str(ii) '/projections.mat']);
angles=importdata(['../1_Measured_data/HEA_' num2str(ii) '/angles.mat']);
cen=round((size(projections,1)+1)/2);
data = importdata(['./Input/Atom_local_classified_HEA_' num2str(ii) '_nanoparticle.mat']);
atoms = data.local_atomtype;
model=(data.curr_model-repmat(cen*ones(3,1),[1,length(atoms)]))*Res;

% process and select input data
pjs = 1:size(projections,3);
angles = angles(pjs,:);
projections = max(projections,0);
projections = projections(2:end,2:end,pjs);
projections=My_paddzero(projections,size(projections)+[50 50 0]);

[N1,N2,num_pj] = size(projections);
halfWidth = 4;  % the cropped bondary size for each atoms
Z_arr = [28 45 78]; % averaged atomic Z number of Type 1 2 3

% initialize refinement parameters
xdata = [];
xdata.Res     = Res;
xdata.Z_arr   = Z_arr;
xdata.halfWidth = halfWidth;
xdata.atoms = atoms;
xdata.model = model;
xdata.angles = angles;

% starting H, B factor and bounds
para0 = [1 1.36 2.58;
    10 10 10];
model_refined = model;

lb=[1, 1, 1;...
    5, 5, 5];
ub=[1, 2, 3;...
    15, 15, 15];
opt = optimset('TolFun',1e-12);

% scan for best initial values
x0 = para0;
x0(1,:)=x0(1,:)/x0(1,1);
xdata.model = model_refined;
xdata.model_ori = model_refined;

% search for optimum B factor by using scanning method
% scan with large steps
num=0;
for h_scan=-0.5:0.2:0.5
    for h_scan2=-1:0.2:1
    for b_scan=-3:0.6:3
        num=num+1;
        x(:,:,num)=x0+[0 h_scan h_scan2;b_scan b_scan b_scan];
    end
    end
end
err_arr=zeros(1,num);
for i=1:num
    [y_pred,para0] = Cal_Bproj_2type(x(:,:,i), xdata, projections);
    err_arr(i)=sum( abs(y_pred(:)-projections(:)) ) / sum(abs(projections(:)));
end
errRscan=[];
[errR,ind]=min(err_arr);
errRscan(end+1)=errR;

% scan with small steps
num=0;
x0=x(:,:,ind);
for h_scan=-0.1:0.04:0.1
    for h_scan2=-0.1:0.04:0.1
    for b_scan=-0.3:0.06:0.3
        num=num+1;
        x(:,:,num)=x0+[0 h_scan h_scan2;b_scan b_scan b_scan];
    end
    end
end
err_arr2=zeros(1,num);
for i=1:num
    [y_pred,para0] = Cal_Bproj_2type(x(:,:,i), xdata, projections);
    err_arr2(i)=sum( abs(y_pred(:)-projections(:)) ) / sum(abs(projections(:)));
end
errRscan=[];
[errR,ind]=min(err_arr2);
errRscan(end+1)=errR;
para0=x(:,:,ind);
[y_pred,~] = Cal_Bproj_2type(para0, xdata, projections);
save(['./Output/' savename 'HB_scan_initial.mat'],'y_pred','para0','errRscan','err_arr','err_arr2')

% repeating the main refinement procedure few times
for jjjj=1:2
x0 = para0;
x0(1,:)=x0(1,:)/x0(1,1);

xdata.model = model_refined;
xdata.model_ori = model_refined;
xdata.projections=projections;

% fit the best H and B factor with above initial values
[para0, ~,~] = lsqcurvefit(@Cal_Bproj_2type2, x0, xdata, projections, lb, ub, opt);
[y_pred,~] = Cal_Bproj_2type(para0, xdata, projections);
save(['./Output/' savename 'HB_Fit_' num2str(jjjj) '.mat'],'y_pred','para0')

% optimized around the above best H and B factor by gradient descend
xdata.projections=[];
xdata.step_sz    = 1;
xdata.iterations = 10;
[y_pred,para0,errR] = gradient_B_2type_difB(para0, xdata, projections);
save(['./Output/' savename 'HB_Fit_itr10_' num2str(jjjj) '.mat'],'y_pred','para0','errR')

% refined atomic coordinates with the best refined B factor by gradient descend
xdata.step_sz    = 1;
xdata.iterations = 50;
[y_pred,para,errR] = gradient_fixHB_XYZ(para0, xdata, projections); model_refined = para(3:5,:);
save(['./Output/' savename 'HB_Fit_xyz_itr50_' num2str(jjjj) '.mat'],'y_pred','para','errR')
end
%% Remove surface atoms by identifying atoms that are on the convex hull of the particle
% To ensure the accuracy of the further analysis on the atomic structures, 
% surface atoms of each HEA nanoparticle were excluded prior to computing structural and chemical orders.

atoms0 = xdata.atoms;
model0 = double(para(3:5,:));

shp = alphaShape(model0(1,:)',model0(2,:)',model0(3,:)',4);
bf = boundaryFacets(shp);
ind = unique(bf(:));
surface = model0(:,ind);
satoms = atoms0(ind);

model = model0; model(:,ind)=[];
atoms = atoms0; atoms(ind)=[];

save(['./Output/Final_atomic_model_' savename 'nanoparticle.mat'],'atoms','model')
