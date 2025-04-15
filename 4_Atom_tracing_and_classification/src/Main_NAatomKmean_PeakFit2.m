%% Classification with known peak centers 
clear
addpath('../../../functions');
inputpath='';
saveprefix='otsu0.9NAatom';
rec = importdata([inputpath 'Recon_ARSP.mat']);
rec=double(rec.reconstruction);

TracingResult = importdata([inputpath 'recon1ysTracing_ARSP3_result.mat']);
atom_pos = TracingResult.TotPosArr(TracingResult.exitFlagArr==0,:)';

atom_pos= (atom_pos / 3) - 2;

b1 = find(atom_pos(1,:)<1 | atom_pos(1,:)>size(rec,1)-1);
b2 = find(atom_pos(2,:)<1 | atom_pos(2,:)>size(rec,2)-1);
b3 = find(atom_pos(3,:)<1 | atom_pos(3,:)>size(rec,3)-1);
bT = union(union(b1,b2),b3);
atom_pos(:,bT) = [];


tight_support = My_obtain_tight_support_ver1(rec);
save support tight_support
% tight_support=importdata('support.mat');
se =strel3d(6);
tight_support  = imdilate(tight_support ,se);  

curr_model0      = atom_pos;

% temp=mod(curr_model0(:),1);
% figure()
% hist(temp,100)


lnorm = 1;
interp_type = 'cubic';
curr_model      = curr_model0;
labels          = ones(1,size(curr_model0,2));
Num_species     = 4;
Num_atom        = size(curr_model,2);
%labels = randi(Num_species, [Num_atom,1]);
O_Ratio = 3;
halfSize = 9/O_Ratio;
plothalfSize = 5/O_Ratio;
ds = 1/O_Ratio;
%% obtain global intensity histogram

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
points(isinf(points))=0;
points_plot = splinterp3(rec, y_set_plot, x_set_plot, z_set_plot);
points_plot(isnan(points_plot))=0;
points_plot(isinf(points_plot))=0;
% integrate intensity (for 3x3x3 and 5x5x5 voxels) for each traced peak
intensity_integ      = sum(points);
intensity_integ_max  = max(points);
% intensity_integ_plot = sqrt(sum(points_plot.^2));
% intensity_integ_plot=intensity_integ;
intensity_integ_plot = sum(points_plot);
% intensity_integ_plot = max(points_plot);

separate_part = 100;
L_forAver = 10;
[hist_inten,cen_integ_total]= hist(intensity_integ,separate_part);

dist_mat = zeros(Num_species,1);

[hist_inten_plot,cen_integ_total_plot]= hist(intensity_integ_plot,separate_part);

% fit histogram with two gaussian
i_guess = [max(hist_inten_plot)/10 cen_integ_total_plot(round(separate_part/6)) cen_integ_total_plot(round(separate_part/10))...
        max(hist_inten_plot) cen_integ_total_plot(round(separate_part/4*2)) cen_integ_total_plot(round(separate_part/10))];
Xdata = cen_integ_total_plot;
Ydata = hist_inten_plot;
[p, fminres, fitresult] = My_two_gaussianfit(Xdata, Ydata,i_guess);  
% [p, fminres, fitresult] = My_two_gaussianfit_constant_ratio(Xdata, Ydata, i_guess);  

% signle gaussian used in fitting
fun = @(p,xdata) p(1)*exp(-((xdata-p(2))/p(3)).^2);
fitresult1=fun(p(1:3),Xdata);
fitresult2=fun(p(4:6),Xdata);

% p([2 5])=[1.8e5 4.8e5];
% p([2 5])=[2.1e5 4.9e5];

ind1=find(abs(intensity_integ_plot-0)<1e4);
ind2=find(abs(intensity_integ_plot-2e5)<1e4);
ind3=find(abs(intensity_integ_plot-3e5)<1e4);
ind4=find(abs(intensity_integ_plot-4e5)<1e4);


avg_atom(:,1) = sum(points(:,ind1),2)/length(ind1);
avg_atom(:,2) = sum(points(:,ind2),2)/length(ind2);
avg_atom(:,3) = sum(points(:,ind3),2)/length(ind3);
avg_atom(:,4) = sum(points(:,ind4),2)/length(ind4);
avg_atom(isnan(avg_atom))=0;

y_up = round(max(hist_inten_plot)/10*12);
iter=0;

while true
    iter=iter+1;
    intensity_integ_1_plot = intensity_integ_plot(labels == 1);
    intensity_integ_2_plot = intensity_integ_plot(labels == 2);
    intensity_integ_3_plot = intensity_integ_plot(labels == 3);
    intensity_integ_4_plot = intensity_integ_plot(labels == 4);

    figure(200+2*halfSize*O_Ratio+1)
    clf
    subplot(Num_species+1,1,1);
    hist(intensity_integ_plot,separate_part);
    hold on
    plot(Xdata, fitresult, 'r-', 'LineWidth',2);
    plot(Xdata, fitresult1, 'g-', 'LineWidth',2);
    plot(Xdata, fitresult2, 'g-', 'LineWidth',2);
    hold off
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
    title(sprintf('%d 4 atoms',sum(labels==4)));    
    
    avg_atom1 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);
    avg_atom2 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);
    avg_atom3 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);
    avg_atom4 = zeros(2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1,2*halfSize*O_Ratio+1);
    
    avg_atom1(useInd) = avg_atom(:,1);
    avg_atom2(useInd) = avg_atom(:,2);
    avg_atom3(useInd) = avg_atom(:,3);
    avg_atom4(useInd) = avg_atom(:,4);
    
    figure(300+2*halfSize*O_Ratio+1);
    clf
    img(sum(avg_atom1,3),'atom1', sum(avg_atom2,3),'atom2',sum(avg_atom3,3),'atom3',sum(avg_atom4,3),'atom4');
    drawnow;
    
    old_labels = labels;
       
    obj = 0;
    for n = 1:Num_atom
        for k=1:Num_species
            dist_mat(k) = norm(points(:,n) - avg_atom(:,k), lnorm).^lnorm/sum(points(:,n).^lnorm);
            if sum(points(:,n).^lnorm)==0
                dist_mat(k)=0;
            end
        end
        [dist2, idx] = min(dist_mat);        
        labels(n) = idx;
        if isinf(dist2)
            n
        end
        obj = obj + dist2;
    end
    
    for k=1:Num_species
        avg_atom(:,k) = sum(points(:,labels==k),2)/sum(labels==k);
    end
    avg_atom(:,1)=0; % constraint nan atom peak to be zero
    avg_atom(isnan(avg_atom))=0;
    fprintf('%02i. obj = %.3f\n',iter, (obj/Num_atom)^(1/lnorm) );
    
    % if there is no change in the atomic specise classification, break
    if ~any(old_labels-labels), break; end
    
end

%% apply support
atomtype=labels;
for i=1:length(atomtype)  %if previous determined noatom is within the tight support, then change back to atom
        Rpos = round(curr_model(:,i));
        if tight_support(Rpos(1),Rpos(2),Rpos(3))>0
            if atomtype(i)==1
                atomtype(i) = 5;
            end
        end
end
labels=atomtype;


save([inputpath sprintf('%s_L%d_Kmean.mat',saveprefix,lnorm)],'curr_model','labels','intensity_integ_plot','cen_integ_total_plot','intensity_integ_max','hist_inten_plot');
%%
% atoms=labels;
% atoms(atoms==5)=1;
% model=curr_model(:,atoms~=1);
% save model model

