
%% ASM taperign and icEEG study: 1st part
% First step for result extraction and data adaptation code

% Input:
%   1.  iEEG_vs_Tap: Set of .mat files (one for each individual) including 
%       a) the parameters of interest extracted for each ROI during 24h of 
%       Steady-State ASM plasmaconcnetration and minimum concentration; 
%       b) date-time vector for each parameters value; c) ROI information 
%       including name and if each ROI is SOZ, RSC or SPK; d) An independent
%       file (DPC_Times_Data.mat) including timinigs of tapering protocol 
%       and miniumn ASM levles for each individual

% Output
%   1.  Figures for the publication relating methodology and TOD resolved
%       results.
%
%   2.  MDL_Taper: File summarizing the effect of ASM tapering on icEEG for
%       all subjects. This file is required for gmb_ASM_EEG_P2.m, but a
%       version has been provided. To generate a new one, uncomment the last 
%       two lines and run the script.

%% Workspace definition

% Clear workspace
clear all 
clc

% Root Path where data is being stored
path_Orig = cd;

% All figures Docked
set(0,"DefaultFigureWindowStyle",'docked')

% Paths for data of interest
path_TapDT  = [path_Orig,'/data/iEEG_vs_Tap'];    % EEG data at ASM states

% File names to exclude form EEG data folder
Excl_IDs = "DPC_Times_Data";

% Avaialble data
lst = dir(path_TapDT);
aux = cellfun(@(x) x(1),{lst.name});
r_dx = strfind(aux,'.');
lst (r_dx) = [];

% Generated individual identifier list
IDs = string(string(erase({lst.name},'.mat')));
if ~isempty(Excl_IDs)
    IDs(logical(sum(IDs == Excl_IDs',1))) = [];
end

% RGB code for the limits for colorbar
% Blue-White-Red colormap
c_lim_0 = [0,0,1;
    1,  1, 1;
    1,0, 0];

cmap_0 = [linspace(c_lim_0(1,1),c_lim_0(2,1),128),linspace(c_lim_0(2,1),c_lim_0(3,1),128);
    linspace(c_lim_0(1,2),c_lim_0(2,2),128),linspace(c_lim_0(2,2),c_lim_0(3,2),128);
    linspace(c_lim_0(1,3),c_lim_0(2,3),128),linspace(c_lim_0(2,3),c_lim_0(3,3),128)];

% White-blak colormap
c_lim_4 = [1,  1, 1;
    0.4, 0.4, 0.4];

cmap_4 = [linspace(c_lim_4(1,1),c_lim_4(2,1),256);
    linspace(c_lim_4(1,2),c_lim_4(2,2),256);
    linspace(c_lim_4(1,3),c_lim_4(2,3),256)];

% Color for boxplot representations
Dt_MkCl= {'#000066','#660000','#662200','#440066','#006600'};
Dt_clr = {'#6666FF','#FF6666','#FF9966','#9966FF','#66FF66'};

% List of available Lobes that ROIs can be part of
Lobes = ["Frontal","Temporal","Parietal","Occipital","Hippocampus","Amygdala","Cingulate","Caudate"];%
% Lobes to not considere: icEEG on Caudate unrelailable
RM_Lobe = "Caudate";

Lb_cnt = zeros(length(Lobes),2);

% For hourly evoulition Time Templates
T_fin = 1; % Resolution of time templates in hours
TOfD= [0:T_fin:24-T_fin];
TMxT= [-12:T_fin:12-T_fin];

% Generate the ASM Profile
ASM_prf = cell(20, 2);
ASM_prf{1,1} = 'Acetazolamide';
ASM_prf{2,1} = 'Amlodipine';
ASM_prf{3,1} = 'Carbamazapine';
ASM_prf{4,1} = 'Clobazam';
ASM_prf{5,1} = 'Clonazepan';
ASM_prf{6,1} = 'Diazepam';
ASM_prf{7,1} = 'Gabapeptin';
ASM_prf{8,1} = 'Lacosamide';
ASM_prf{9,1} = 'Lamotrigine';
ASM_prf{10,1} = 'Levetiracetam';
ASM_prf{11,1} = 'Lorazepan';
ASM_prf{12,1} = 'Oxacarbazepine';
ASM_prf{13,1} = 'Perampanel';
ASM_prf{14,1} = 'Phenobarbital';
ASM_prf{15,1} = 'Phenytoin';
ASM_prf{16,1} = 'Pregabalin';
ASM_prf{17,1} = 'Sodium Valporate';
ASM_prf{18,1} = 'Topiramate';
ASM_prf{19,1} = 'Trobalt';
ASM_prf{20,1} = 'Zonisamide';
ASM_prf{1,2} = 'OT';%CA
ASM_prf{2,2} = 'OT';%CC
ASM_prf{3,2} = 'SC';
ASM_prf{4,2} = 'GB';
ASM_prf{5,2} = 'GB';
ASM_prf{6,2} = 'GB';
ASM_prf{7,2} = 'OT';%CC
ASM_prf{8,2} = 'SC';
ASM_prf{9,2} = 'SC';
ASM_prf{10,2} = 'SV';
ASM_prf{11,2} = 'GB';
ASM_prf{12,2} = 'SC';
ASM_prf{13,2} = 'OT';%GT
ASM_prf{14,2} = 'GB';
ASM_prf{15,2} = 'SC';
ASM_prf{16,2} = 'OT';%CC
ASM_prf{17,2} = 'ML';
ASM_prf{18,2} = 'ML';
ASM_prf{19,2} = 'OT';%PC
ASM_prf{20,2} = 'ML';
ASM_prf = convertCharsToStrings(ASM_prf);

%% Extract data of interest

load(fullfile(path_TapDT,'DPC_Times_Data.mat'));
TapSrg_Info = TapSrg_Info;

%--------------------------------------------------------------------------
% All the features included
% Extract Drug Plasma Concentration data
DT = load(fullfile(path_TapDT,[char(IDs(1)),'.mat']));
% The included features
Feat_CD = fieldnames(DT.Max_DPC);
Feat_CD(strcmp(Feat_CD,'t')) = [];

%--------------------------------------------------------------------------
% All included ROIs (considering hemispheres)
ROIs = [];
for i=1:length(IDs)
    DT = load(fullfile(path_TapDT,[char(IDs(i)),'.mat']));
    ROIs = [ROIs;DT.ROI_Info.ROI];
end

ROIs = unique(ROIs);

% Initialize Model for Averaged Study
clear MDL_df
MDL_df.ID = [];
MDL_df.ROI = [];
MDL_df.LOB = [];
MDL_df.SOZ = [];
MDL_df.RSC = [];
MDL_df.SPK = [];

for i=1:length(Feat_CD)
    MDL_df.(Feat_CD{i}) = [];
end

% Initialize Model for TOD resolved Study
MDL_tof_dt = MDL_df;
MDL_tof_dt.HOUR = [];
MDL_tmd_dt = MDL_tof_dt;

MDL_data = MDL_df;
MDL_data.State = [];

MDL_df.DPCd = [];
MDL_df.SRGt = [];
MDL_df.TAPt = [];

for i=1:length(IDs)
    %% Obtain icEEG data
    DT = load(fullfile(path_TapDT,[char(IDs(i)),'.mat']));

    % Time variables
    Mx_t = DT.Max_DPC.t;
    Mn_t = DT.Min_DPC.t;

    % From Date-Time to hous
    a = (hour(Mx_t)+(minute(Mx_t)+second(Mx_t)/60)/60);
    b = (hour(Mn_t)+(minute(Mn_t)+second(Mn_t)/60)/60);

    % Same resolution as Time Templates for rearrangement
    tx_o = floor(a/T_fin)*T_fin;
    tn_o = floor(b/T_fin)*T_fin;

    %% Asing 
    [~,idx] = min(abs(a-b'));
    tn_p = floor(datenum(Mn_t-Mn_t(1))*24/T_fin)*T_fin-12;
    tx_p = tn_p(idx);

    MDL_df.ID     = [MDL_df.ID; repmat(categorical(IDs(i)),length(DT.ROI_Info.ROI),1)];
    MDL_tof_dt.ID = [MDL_tof_dt.ID;repmat(categorical(IDs(i)),length(TOfD)*length(DT.ROI_Info.ROI),1)];

    % Add the ROIs
    MDL_df.ROI   = [MDL_df.ROI;categorical(DT.ROI_Info.ROI)];
    MDL_tof_dt.ROI = [MDL_tof_dt.ROI;repmat(categorical(DT.ROI_Info.ROI),length(TOfD),1)];

    % Add the Lobes
    MDL_df.LOB   = [MDL_df.LOB;categorical(DT.ROI_Info.LOB)];
    MDL_tof_dt.LOB = [MDL_tof_dt.LOB;repmat(categorical(DT.ROI_Info.LOB),length(TOfD),1)];

    % Add plasma concentration drop and time from surgery    
    idx = string(TapSrg_Info.Var1) == IDs(i);
    MDL_df.DPCd = [MDL_df.DPCd;TapSrg_Info.Var3(idx)+zeros(size(DT.ROI_Info.LOB))];
    MDL_df.SRGt = [MDL_df.SRGt;TapSrg_Info.Var2(idx)+zeros(size(DT.ROI_Info.LOB))];
    MDL_df.TAPt = [MDL_df.TAPt;TapSrg_Info.Var4(idx)+zeros(size(DT.ROI_Info.LOB))];
    
    % Count Lobes
    for j=1:length(Lobes)
        jdx = string(DT.ROI_Info.LOB) == Lobes(j);
        Lb_cnt(j,1) = Lb_cnt(j,1) + any(jdx); % By subjects
        Lb_cnt(j,2) = Lb_cnt(j,2) + sum(jdx); % By ROIs per subject
    end

    % Add the SOZ & RSC
    MDL_df.SOZ = [MDL_df.SOZ;DT.ROI_Info.SOZ];
    MDL_df.RSC = [MDL_df.RSC;DT.ROI_Info.RSC];
    MDL_df.SPK = [MDL_df.SPK;DT.ROI_Info.SPK];

    MDL_tof_dt.SOZ = [MDL_tof_dt.SOZ;repmat(DT.ROI_Info.SOZ,length(TOfD),1)];
    MDL_tof_dt.RSC = [MDL_tof_dt.RSC;repmat(DT.ROI_Info.RSC,length(TOfD),1)];
    MDL_tof_dt.SPK = [MDL_tof_dt.SPK;repmat(DT.ROI_Info.SPK,length(TOfD),1)];

    % Add the hours of the data
    asdf = repmat(TOfD,length(DT.ROI_Info.RSC),1);
    MDL_tof_dt.HOUR = [ MDL_tof_dt.HOUR;asdf(:)];

    %%

    for j=1:length(Feat_CD)
        %--------------------------------------------------------------------------
        % TOfD data by Subject & ROIs
        aux = [];
        bux = [];
        for k=1:length(TOfD)
            % Index for the Hour @ each segment
            idx_1 = tx_o == TOfD(k);
            idx_2 = tn_o == TOfD(k);

            if any(idx_1) && any(idx_2) % If both have some data
                % Average for that Hour and then difference
                aux(:,k) = nanmean(DT.Min_DPC.(Feat_CD{j})(:,idx_2),2)-...
                    nanmean(DT.Max_DPC.(Feat_CD{j})(:,idx_1),2);
            else
                aux(:,k) = nan(length(DT.ROI_Info.ROI),1);
            end
            bux = [bux;aux(:,k)];

        end

        % Store data for regresion model each ROI individually
        MDL_tof_dt.(Feat_CD{j}) = [ MDL_tof_dt.(Feat_CD{j});bux];

        % Store data for regresion model
        asdf = nanmean(aux,2);%./nanstd(aux')';
        asdf(isinf(asdf)) = nan;
        MDL_df.(Feat_CD{j})     = [ MDL_df.(Feat_CD{j});asdf];

    end
end

%% Add the tapered ASM type
ASM_cnt.CA = 0;
MDL_df.Valp = categorical(zeros(size(MDL_df.Alpha)));
for i=1:length(IDs)
    DT = load(fullfile(path_TapDT,[char(IDs(i)),'.mat']));

    idx = MDL_df.ID == IDs(i);
    
    for j=1:length(DT.Tap_Drg)
        jdx = string(DT.Tap_Drg{j})==ASM_prf(:,1);
        asm_typ = ASM_prf(jdx,2);
        
        % Add the ASM tapering to the dataset
        if ~any(strcmp(fieldnames(MDL_df),asm_typ))
            MDL_df.(asm_typ) = categorical(zeros(size(MDL_df.Alpha)));
        end

        MDL_df.(asm_typ)(idx) = categorical(1);

        % Add Special case: Sodium Valporate (pos: 17)
        if find(jdx)==17
            MDL_df.Valp(idx) = categorical(1);
        end

        % Count it
        if ~any(strcmp(fieldnames(ASM_cnt),asm_typ))
            ASM_cnt.(asm_typ) = 0;
        end
        ASM_cnt.(asm_typ) = ASM_cnt.(asm_typ)+1;
    end
end

% Structures to tables for Model fitting
MDL_df = struct2table(MDL_df);
MDL_tof_dt = struct2table(MDL_tof_dt);

% Reorder the Lobes categorical
Lobes(Lobes==RM_Lobe) = [];

MDL_df(MDL_df.LOB==RM_Lobe,:) = [];
MDL_df.LOB = removecats(MDL_df.LOB,RM_Lobe);
MDL_df.LOB = reordercats(MDL_df.LOB,Lobes);

MDL_tof_dt(MDL_tof_dt.LOB==RM_Lobe,:) = [];
MDL_tof_dt.LOB = removecats(MDL_tof_dt.LOB,RM_Lobe);
MDL_tof_dt.LOB = reordercats(MDL_tof_dt.LOB,Lobes);

%% Example of How difference is extracted

% Figure 1: Corresponding to Figure 2 of the Main text
% Log δ band power difference estimation process for a single subject accounting 
% for time-of-day effects. (A) Baseline-ASM data matrix and (B) Reduced-ASM 
% data matrix, where each row represents a ROI, each column an hour of the 
% day and the grey-scale the log δ band power. (C) Colour-coded difference 
% matrix between reduced-ASM relative to baseline-ASM, where red represent 
% an increased relative to band power, blue a decrease and white no change. 
% (D) Average of (C) over time (one value per ROI) following the same 
% colour-code.

figure(1);

% individual and parameter to show the extraction process
SW_id = '1020';
SW_ft = 'Delta';

aux = [];

% Load the data
DT = load(fullfile(path_TapDT,[char(SW_id),'.mat']));

% Data for the maximum concentration
Mx_t = DT.Max_DPC.t;
Mx_feat = DT.Max_DPC.(SW_ft)(3:end,:);

% Data for the minimum concnetration
Mn_t = DT.Min_DPC.t;
Mn_feat = DT.Min_DPC.(SW_ft)(3:end,:);

% To store the hourlly average
hMx = [];
hMn = [];

tx_o = floor((hour(Mx_t)+(minute(Mx_t)+second(Mx_t)/60)/60)/T_fin)*T_fin;
tn_o = floor((hour(Mn_t)+(minute(Mn_t)+second(Mn_t)/60)/60)/T_fin)*T_fin;

for i=1:length(TOfD)
    idx_1 = tx_o == TOfD(i);
    idx_2 = tn_o == TOfD(i);

    hMx(:,i) = nanmean(Mx_feat(:,idx_1),2);
    hMn(:,i) = nanmean(Mn_feat(:,idx_2),2);
end

% Colormpas where colums are each hour of the day and rows are ROIs
% ---------------------------------------------------------------------
% Maximum/steady-state ASM values
S1 = subplot(3,10,[1:9]);
imagesc(hMx)
yticks([]);
ylabel('ROI')
xticks([1,13,24])
xticklabels({'0','hod','23'})
colorbar;
aux = clim;
title('Pre-Tapering period')

% Minumum ASM values
S2 = subplot(3,10,[1:9]+10);
imagesc(hMn)
xticks([1,13,24])
xticklabels({'0','hod','23'})
yticks([]);
ylabel('ROI')
colorbar;
aux(2,:) = clim;
title('Max-Tapering period')

clim(S1,[min(aux(:)) max(aux(:))])
clim(S2,[min(aux(:)) max(aux(:))])

colormap(S1,cmap_4')
colormap(S2,cmap_4')

% Difference along the TOD
S3 = subplot(3,10,[1:9]+20);
hM_df = hMn-hMx;
hM_df(isnan(hM_df)) = 0;
imagesc(hM_df)
xticks([1,13,24])
xticklabels({'0','hod','23'})
yticks([]);
ylabel('ROI')
colorbar;
aux = clim;
title('TOD matched difference')

% ---------------------------------------------------------------------
% Average Difference along TOD obtainin a single column
S4 = subplot(3,10,30);
imagesc(nanmean(hM_df,2))
yticks([]);
xticks(1);
xticklabels("Avrg.");
aux(2,:) = clim;

clim(S3,max(abs(aux(:)))*[-1 1])
clim(S4,max(abs(aux(:)))*[-1 1])

colormap(S3,cmap_0')
colormap(S4,cmap_0')

%% TOD resolved differences comparing Cortical and Sub-cortical regions

% Figure 2: Corresponding to Figure 5 of the Main text
% Log10(δ band power) difference relative to baseline-ASM resolved along the 
% hours of the day. Solid lines represent the mean across ROIs and subjects, 
% and shadows represent the 95% confidence interval across subjects, for 
% cortical (orange) and sub-cortical (green) regions.

figure(2)
clf

A = [];
B = [];
C = [];

% Feature to sztudy
feat = 'Delta';

% Only Hippocampus and Amygdala are Sub-Cortical
kdx = MDL_tof_dt.LOB == "Hippocampus" |...
            MDL_tof_dt.LOB == "Amygdala";

% Average of all ROIS for each subject for all Lobes and dividing cortical
% and subcortical
for i=1:length(IDs)
    idx = MDL_tof_dt.ID == IDs(i);
    for j=1:length(TOfD)
        jdx = MDL_tof_dt.HOUR == TOfD(j);
        A(i,j) = nanmean(MDL_tof_dt.(feat)(idx&jdx));
        
        B(i,j) = nanmean(MDL_tof_dt.(feat)(idx&jdx&~kdx));
        C(i,j) = nanmean(MDL_tof_dt.(feat)(idx&jdx&kdx));
    end
end

% All ROIs averaged
S1 = subplot(3,1,2);
[P,H] = gmb_ShadowPlot(A,TOfD,Dt_clr(1),0.025);
H.Color = Dt_MkCl{1};

yline(0)
YL = max(abs(ylim));
ylabel('log(BPW)')

xlim([TOfD(1),TOfD(end)])
xticks([0,8,12,20])
xline(8,"LineStyle",":")
xline(20,"LineStyle",":")
xlabel('Hour Of Day')

title('TOD resolved log-band-power')

% Comparing cortical and subcortical ROIS
S2 = subplot(3,1,3);
[p(1),h(1)] = gmb_ShadowPlot(B,TOfD,Dt_clr(3),0.025);
h(1).Color = Dt_MkCl{3};
hold on;
[p(2),h(2)] = gmb_ShadowPlot(C,TOfD,Dt_clr(5),0.025);
h(2).Color = Dt_MkCl{5};

xlim([TOfD(1),TOfD(end)])
xticks([0,8,12,20])
xline(8,"LineStyle",":")
xline(20,"LineStyle",":")
xlabel('Hour Of Day')

yline(0)
ylabel('log(BPW)')
YL = max(abs([YL, ylim]));
ylim (YL*[-1, 1])

legend(S1,[H,P],{'mean','CI'})
legend(S2,p,{'Cortical', 'Sub-cortical'})

subplot(3,1,2)
ylim (YL*[-1, 1])

%% Genrate data set for Part 2

%file = [cd,'/data/MDL_Taper.mat'];
%save(file,'MDL_df')
