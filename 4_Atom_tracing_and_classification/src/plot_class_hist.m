function [peak_info,intensity_arr,intensity_plot_arr] = plot_class_hist(...
    RecVol_padded, temp_model, temp_type, classify_info)
% Yao Yang UCLA, 2019.05.16
% from class7 to classify two type

% two step: first classify non-atoms from the originally method
% (use three type, which already modified)
% then select certain number of atoms around 3 peaks, calculate the
% averaged atom size then move them freely

% Modified with Minh's implement:
% 1. reduced several sub-functions
% 2. add L1 and L2 norm options
% Notice: temp_model should be (3 by N)

if isfield(classify_info,'halfSize')       halfSize = classify_info.halfSize; %#ok<SEPEX>
else halfSize = 1;          %#ok<SEPEX>
end
if isfield(classify_info,'plothalfSize')   plothalfSize = classify_info.plothalfSize; %#ok<SEPEX>
else plothalfSize = 4;      %#ok<SEPEX>
end
if isfield(classify_info,'SPHyn')          SPHyn = classify_info.SPHyn; %#ok<SEPEX>
else SPHyn = 1;             %#ok<SEPEX>
end
if isfield(classify_info,'separate_part')  separate_part = classify_info.separate_part; %#ok<SEPEX>
else separate_part = 100;   %#ok<SEPEX>
end
if isfield(classify_info,'PLOT_YN')        PLOT_YN = classify_info.PLOT_YN; %#ok<SEPEX>
else PLOT_YN = 0;           %#ok<SEPEX>
end

Num_types = numel(unique(temp_type(temp_type>0)));
peak_info = zeros(Num_types+2,separate_part);
%% obtain global intensity histogram

[xX,yY,zZ] = ndgrid(-halfSize:halfSize,-halfSize:halfSize,-halfSize:halfSize);
SphereInd = find(xX.^2+yY.^2+zZ.^2 <=(halfSize+0.5)^2);
[xXp,yYp,zZp] = ndgrid(-plothalfSize:plothalfSize,-plothalfSize:plothalfSize,-plothalfSize:plothalfSize);
SphereInd_plot = find(xXp.^2+yYp.^2+zZp.^2 <=(plothalfSize+0.5)^2);

if SPHyn
  useInd = SphereInd;
  useInd_plot =SphereInd_plot;
else
  useInd = 1:length(xX);
  useInd_plot = 1:length(xXp);
end

intensity_arr = zeros(numel(useInd),size(temp_model,2));
intensity_plot_arr = zeros(numel(useInd_plot),size(temp_model,2));
%% compute initial average Fe atom, Pt atom and non atom

% re-produce array for integrated intensity
intensity_integ = zeros(1,size(temp_model,2));
intensity_integ_plot = zeros(1,size(temp_model,2));

% integrate intensity (for 3x3x3 and 5x5x5 voxels) for each traced peak
for j=1:size(temp_model,2)
    curr_pos  = round(temp_model(:,j));
    box_integ = RecVol_padded(curr_pos(1)-halfSize:curr_pos(1)+halfSize,...
        curr_pos(2)-halfSize:curr_pos(2)+halfSize,curr_pos(3)-halfSize:curr_pos(3)+halfSize);
    intensity_integ(j) = sum(box_integ(useInd));
    intensity_arr(:,j) = box_integ(useInd);
end
for j=1:size(temp_model,2)
    curr_pos  = round(temp_model(:,j));
    box_integ = RecVol_padded(curr_pos(1)-plothalfSize:curr_pos(1)+plothalfSize,...
        curr_pos(2)-plothalfSize:curr_pos(2)+plothalfSize,curr_pos(3)-plothalfSize:curr_pos(3)+plothalfSize);
    intensity_integ_plot(j) = sum(box_integ(useInd_plot));
    intensity_plot_arr(:,j) = box_integ(useInd_plot);
end

    [hist_inten_plot,cen_integ_total_plot]= hist(intensity_integ_plot,separate_part);
    peak_info(1,:) = cen_integ_total_plot;
    peak_info(2,:) = hist_inten_plot;
    y_up = round(max(hist_inten_plot)/10)*12;
    if PLOT_YN
        %     figure(103)
        figure
        clf
        set(gcf,'Position',[1100,0,400,900])
        subplot(Num_types+1,1,1);
        hist(intensity_integ_plot,separate_part);
        xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
        ylim([0 y_up]);
        xlabel('integrated intensity (a.u.)');
        ylabel('# atoms');
        title(sprintf('boxsize %d',halfSize*2+1));
    end
    for i = 1:Num_types
        intensity_integ_sub = intensity_integ_plot(temp_type==i);
        [hist_inten_plot_sub,~]= hist(intensity_integ_sub,cen_integ_total_plot);
        peak_info(2+i,:) = hist_inten_plot_sub;
        if PLOT_YN
            subplot(Num_types+1,1,i+1)
            hist(intensity_integ_sub,cen_integ_total_plot);
            xlabel('integrated intensity (a.u.)');
            ylabel('# atoms');
            title(sprintf('%d Type %i atoms',sum(temp_type==i),i));
            ylim([0 y_up]);
            xlim([0 ceil(max(intensity_integ_plot)/5)*5]);
        end
    end
end