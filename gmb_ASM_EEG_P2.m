
%% ASM taperign and icEEG study: 2nd part
% Second step for result extraction. Required inputs are provided so executing 
% 1st part is not neccessary

% Input:
%   1.  MDL_nTaper & MDL_Taper: Tables containing ROI wise information of
%       ASM effects on icEEG recordings for non-Tapered and Tapered
%       individuals respectivelly'

% Output
%   1.  Figures for the publication corresponding to resutls and
%       supplemenraty sections
%
%   2.  Multiple variables describing the statistical analysis:
%       a)  RG_p_val & RG_ef_sz: P-values and effect sizes studing if there
%           is a statisticaly significant band power reduction (first row) 
%           and its correlation with ASM reduction (second row) considering
%           different cortical segmentations
%       a)  SJ_p_val & SJ_ef_sz: P-values and effect sizes studing if there
%           is a statisticaly significant band power reduction (first row) 
%           and its correlation with ASM reduction (second row) after
%           averaging all RIOs fisrt and all subjects second
%       a)  FB_p_val & FB_ef_sz: P-values and effect sizes studing if there
%           is a statisticaly significant band power reduction (first row) and
%           relative band power reduction (second row).
%
%   3.  hier_mdl: Variable contaitng the multple hierarchical models tested
%       during the study to confirm our findings.

%% Obtain data and adapt it to extract results

% Clear workspace
clear all
clc

% All figures docked
set(0,"DefaultFigureWindowStyle",'docked')

% Load and join tapered and non-tapered patietns data
load([cd,'/data/MDL_Taper.mat'])
load([cd,'/data/MDL_nTaper.mat'])
MDL=[MDL_df;MDL_df_0];

% Lobe Names
Lobe = ["Frontal","Cingulate","Parietal", "Temporal","Occipital",...
    "Amygdala", "Hippocampus", "Caudate"];

% Remove cuadate from the dataset
RM_lb = "Caudate";%
for i=1:length(RM_lb)
    idx = MDL.LOB == RM_lb(i);
    MDL(idx,:) = [];
    MDL.LOB = removecats(MDL.LOB,RM_lb(i));
    Lobe(strcmp(Lobe,RM_lb(i))) = [];
end
MDL.LOB = reordercats(MDL.LOB,Lobe);

% Determine if Lobe corresponds to Cortical structures
Lob_Depth = MDL.LOB=="Hippocampus" | MDL.LOB=="Amygdala";
MDL.Cort = categorical(~Lob_Depth);
MDL.SOZ = categorical(MDL.SOZ);

for i=1:length(Lobe)
    idx = MDL.LOB == Lobe(i);
    MDL.(Lobe(i)) = categorical(idx);
end

% Available Frequency bands
Wave = ["Delta", "Theta", "Alpha", "Beta","Gamma"];

% Avaialable ASM classess
ASM_tp = ["SC", "ML", "GB", "SV", "OT"];%

% Color code
Dt_MkCl= {'#000066','#660000','#662200','#440066','#006600'};
Dt_clr = {'#6666FF','#FF6666','#FF9966','#9966FF','#66FF66'};

%% All subjects: distribution and correlation with tapering strength

% Only Tapered patients: minimum plasma concnetration bellow 1
MDL_cr = MDL(MDL.DPCd~=1,:);

% Figure 3: Corresponding to Figure 3 on the Main text
% Tapering effects on log(δ band power) for all subjects and ROIs. 
% (A) Histogram across all ROIs and subjects of log(δ band power) difference 
% relative to baseline (each dot equates to 3 datapoints). (B) Scatter plot 
% between log(δ band power) difference relative to baseline (y-axis) and the 
% recorded ASM percentage at the point of strongest tapering (black circle 
% in Figure 1B) in each subject (x-axis). Data for all ROIs in a single subject 
% line up on one value of ASM tapering percentage for that subject. 
% For each subject, the underlying data is equivalent to Figure 2D.

% Figure 4: Corresponding to Figure 4 on the Main text
% ASM tapering reduces δ band power in cortical but not sub-cortical regions. 
% Plotting convention is the same as Figure 3. (A) The difference between the 
% log(δ band power) before and during tapering in cortical regions for each 
% subject. (B) The difference between the log(δ band power) before and during 
% tapering in sub-cortical regions for each subject. In our data, only 
% Amygdala and Hippocampus are recorded.

i=1;
Fg = [3,4,4];
R = [4,1,4];
subC_cnt = MDL_cr.Cort=="false";
dX = {1:length(subC_cnt),~subC_cnt,subC_cnt};
C = [1,3,5];
dist = 3;

for k=1:length(Fg)
    figure(Fg(k))
    subplot(2,3,R(k))

    % Get the distridution of Band Power differences
    h = histfit(MDL_cr.(Wave(i))(dX{k}),[],'kernel');
    
    % Curve form of distribution
    CurveX = h(2).XData;
    CurveY = h(2).YData;
    
    % Bar form of distribution ad points
    BarX = h(1).XData;
    BarY = h(1).YData;

    cla;
    hold on;

    DP_Y = [];
    DP_X = [];
    for j=1:length(BarY)
        aux = ceil(BarY(j)/dist);
        DP_Y = [DP_Y,(1:aux)*dist];
        DP_X = [DP_X,zeros(1,aux)+BarX(j)];
    end

    scatter(DP_X,DP_Y,'.','MarkerEdgeColor',Dt_MkCl{C(k)})
    F = fill(CurveX,CurveY,'b');
    F.FaceAlpha = 0.5;
    F.FaceColor = Dt_clr{C(k)};
    F.EdgeColor = Dt_MkCl{C(k)};

    xline(0,"LineStyle",':')

    view([-90 90])
    xlim([-3,2])
    xlabel('log(BPW)')
    ylim([1 60])
    ylabel('N. ROIs')
    title('BPW');

    % Correlation between ASM tapering strenght and Band Power reduction
    subplot(2,3,R(k)+[1,2])

    cla
    scatter (MDL_cr.DPCd(dX{k})*100,MDL_cr.(Wave(i))(dX{k}),"MarkerEdgeColor",Dt_MkCl{C(k)})
    H = lsline;
    H.Color = Dt_MkCl{C(k)};
    yline(0,"LineStyle",':')
    ylim([-3,2])
    xlabel('%')
    title('BPW & DPC')
end

%% Statictical analysis 
% for Band Power reduction, Regional differences and correlation with ASM tpaeinrg strength

% Only Tapered patients: minimum plasma concnetration bellow 1
MDL_cr = MDL(MDL.DPCd~=1,:);

% Adapt data for analysis
I = unique(MDL_cr.ID);
for i=1:length(I)
    A2(i) = nanmean(MDL_cr.Delta(MDL_cr.ID == I(i)));
    B2(i) = nanmean(MDL_cr.DPCd(MDL_cr.ID == I(i)));
end

rdx = isnan(A2);
A2(rdx) = [];
B2(rdx) = [];

for i=1:3
    %-----------------------------------------------------
    % No average across subjects
    switch i
        case 1
            % All ROIs
            A = MDL_cr.Delta;
            B = MDL_cr.SRGt;
        case 2
            % Cortical Structures
            idx = MDL_cr.Cort == "true";
            A = MDL_cr.Delta(idx);
            B = MDL_cr.SRGt(idx);
        case 3
            % Subcortical Structures
            idx = MDL_cr.Cort == "false";
            A = MDL_cr.Delta(idx);
            B = MDL_cr.SRGt(idx);
    end
    
    rdx = isnan(A);
    A(rdx) = [];
    B(rdx) = [];

    % Wilcoxon signed ranked test
    [RG_p_val(1,i),~,stat] = signrank(A);
    RG_ef_sz(1,i) = stat.zval/sqrt(length(A));
    
    % Spearman ranked correlation for DPCd vs Bndpower changes
    [RG_ef_sz(2,i),RG_p_val(2,i)] = corr(A,B,'Type','Spearman');

    %-----------------------------------------------------
    % Averaged across subjects
    if i==1 % Only for all ROIs
        % Wilcoxon signed ranked test
        [SJ_p_val(1),~,stat] = signrank(A2);
        SJ_ef_sz(1) = stat.zval/sqrt(length(A2));

        % Spearman ranked correlation for DPCd vs Bndpower changes
        [SJ_ef_sz(2),SJ_p_val(2)] = corr(A2',B2','Type','Spearman');
    end
    
end

%% Hierarchical modeling for multiple scenarions
% Fixed Effect to apply
fixed = {['DPCd'],...
        ['DPCd*Cort'],...
        ['DPCd*LOB'],...
        ['SOZ*Cort'],...
        ['RSC*Cort'],...
        ['SPK*Cort'],...
        ['DPCd*(',char(strjoin(ASM_tp,'+')),'-1)'],...
        ['(DPCd+SRGt)*Cort'],...
        ['(DPCd+SRGt)*Cort']};

response = 'Delta';

% Mixed effect model fitting
for i=1:length(fixed)
    % Add fixed effect to model formula: Random effect on subjec identifier
    form = [response,' ~1 +',fixed{i},'+(1|ID)'];
    
    % Adapt the dataset for each case
    if i==3 || i==7
        % Only Tapered patients, and cortical structures
        MDL_cr = MDL(MDL.DPCd~=1 & MDL.Cort=="true",:);
        MDL_cr.LOB = removecats(MDL_cr.LOB,["Amygdala","Hippocampus"]);
        MDL_cr.LOB = reordercats(MDL_cr.LOB,...
            ["Occipital","Temporal","Frontal","Cingulate","Parietal"]);

    elseif i==9
        % both Tapered and non-Tapered Subjects
        MDL_cr = MDL;

    else
        % Only Tapered subjects
        MDL_cr = MDL(MDL.DPCd~=1,:);

    end
    
    % Remove non usable data
    rm = isnan(MDL_cr.Delta);
    MDL_cr(rm,:) = [];

    % Fit the model with the data and defined fixed effect
    hier_mdl{i} = fitlme(MDL_cr, form);
end


%% Test Other Frecuecnys and if there is a spectral shift

% Figure 5: Corresponding to Figure S3.1 on the Supplementary material
% Distribution of δ, θ, α, β and γ (Top Row) band power and (Bottom Row) 
% relative band power changes due to tapering for all ROIs and subject

FB_p_val = [];
FB_ef_sz = [];

MDL_cr = MDL(MDL.DPCd~=1,:);
figure(5)
clf

for i=1:length(Wave)
    % Band Power For the specific frecuency band
    A = MDL_cr.([char(Wave(i))]);
    A(isnan(A)) = [];
    
    % Wilcoxon signed ranked test
    [FB_p_val(i,1),~,stat] = signrank(A);
    FB_ef_sz(i,1) = stat.zval/sqrt(length(A));
    
    % Plot the distribution
    subplot(2,length(Wave),i)

    h = histfit(MDL_cr.([char(Wave(i))]),[],'kernel');

    cCurveX = h(2).XData;
    cCurveY = h(2).YData;

    cla;
    hold on;

    Fc = fill(cCurveX,cCurveY,'b');
    Fc.FaceAlpha = 0.5;
    Fc.FaceColor = Dt_clr{1};
    Fc.EdgeColor = Dt_MkCl{1};

    xline(0,"LineStyle",':')

    view([-90 90])
    xlim(3*[-1,1])
    xlabel('log10(band power)')
    ylim([1, 150])
    title(Wave(i));

    % Relative Band Power For the specific frecuency band
    A = MDL_cr.([char(Wave(i)),'_rel']);
    A(isnan(A)) = [];
    
    % Wilcoxon signed ranked test
    [FB_p_val(i,2),~,stat] = signrank(A);
    FB_ef_sz(i,2) = stat.zval/sqrt(length(A));
    
    % Plot the distribution
    subplot(2,length(Wave),length(Wave)+i)

    h = histfit(MDL_cr.([char(Wave(i)),'_rel']),[],'kernel');

    cCurveX = h(2).XData;
    cCurveY = h(2).YData;

    cla;
    hold on;

    Fc = fill(cCurveX,cCurveY,'b');
    Fc.FaceAlpha = 0.5;
    Fc.FaceColor = Dt_clr{1};
    Fc.EdgeColor = Dt_MkCl{1};

    xline(0,"LineStyle",':')

    view([-90 90])
    xlim(0.08*[-1,1])
    xlabel('Relative band power')
    ylim([1, 150])
    ylabel('N. ROIs')

end







