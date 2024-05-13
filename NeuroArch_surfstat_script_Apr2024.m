% script for vertex-wise cortical volume analysis (SurfStat) in "Examining the neurostructural architecture of intelligence: The Lothian Birth Cohort 1936 Study"
% Author: Danielle Page

% READ DATA 
avsurf = SurfStatReadSurf( { '/home/Danielle/surf/lh.pial' '/home/Danielle/surf/rh.pial' });
%load avsurf.mat
mask = SurfStatROI([0;24;5],11, avsurf) == 0 & SurfStatROI([0;-29;10],18, avsurf) == 0 & SurfStatROI([0;-12;2],29, avsurf) == 0 & SurfStatROI([0;-1;15],16, avsurf) == 0 & SurfStatROI([0;11;9],16, avsurf) == 0 & SurfStatROI([0;12;-1],12, avsurf) == 0;
[lbc36no sex agedays speed cryst memory visuo g bi_speed bi_cryst bi_memory bi_visuo bi_g icv VolLeft VolRight SALeft SARight ThkLeft ThkRight] = textread('matlabDF_2_clean_edit' , ' %s %s %f %f %f %f %f %f %f %f %f %f %f %f %s %s %s %s %s %s'); 


% read volumetrics if not already
if ~exist('Vol', 'var')
    Vol = SurfStatReadData ([VolLeft  VolRight ]);
end

% CREATE TERMS 
Lbc36no = term(lbc36no);
Sex = term(sex);
Agedays = term(agedays);
Speed = term(speed);
Cryst = term(cryst);
Memory = term(memory);
Visuo = term(visuo);
G = term(g);
Bi_speed = term(bi_speed);
Bi_cryst = term(bi_cryst);
Bi_memory = term(bi_memory);
Bi_visuo = term(bi_visuo);
Bi_g = term(bi_g);
Icv = term(icv);

% INDIVIDUALLY MODELLED DOMAINS & G %
% SET UP MODELS M1-M6
M1 = 1 + Agedays + Sex + Icv + Speed;
M2 = 1 + Agedays + Sex + Icv + Cryst;
M3 = 1 + Agedays + Sex + Icv + Memory;
M4 = 1 + Agedays + Sex + Icv + Visuo;
M5 = 1 + Agedays + Sex + Icv + G;
M6 = 1 + Agedays + Sex + Icv + Speed + Cryst + Memory + Visuo + G;

% bounding box size for surface plots
bbox = [0 0 8 6];

cmap = winter_autumn(420);

% SINGLE EFFECTS % SIMPLE VOLUME MODELS AGE, SEX, ICV CORRECTED
%PROCESSING SPEED%
slm1 = SurfStatLinMod(Vol,  M1,  avsurf);
slm1a = SurfStatT(slm1,  speed);
SurfStatView(slm1a.t.*mask,  avsurf,  'Processing speed on Volume t-values (df = 617)');  % gives t maps
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Speed_Vol_T.png');
SurfStatView(SurfStatP(slm1a, mask),  avsurf,  'Processing speed on Volume p-values (df = 617)');  % gives p maps & RFT clusters
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Speed_Vol_P.png');
SurfStatView(SurfStatQ(slm1a, mask),  avsurf,  'Processing speed on Volume q-values (df = 617)');  % gives FDR q maps
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Speed_Vol_Q.png');



%CRYSTALLISED%
slm2 = SurfStatLinMod(Vol,  M2,  avsurf);
slm2a = SurfStatT(slm2,  cryst);
SurfStatView(slm2a.t.*mask,  avsurf,  'Crystallised ability on Volume t-values (df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Cryst_Vol_T.png');
SurfStatView(SurfStatP(slm2a, mask),  avsurf,  'Crystallised ability on Volume p-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Cryst_Vol_P.png');
SurfStatView(SurfStatQ(slm2a, mask),  avsurf,  'Crystallised ability on Volume q-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Cryst_Vol_Q.png');

%MEMORY%
slm3 = SurfStatLinMod(Vol,  M3,  avsurf);
slm3a = SurfStatT(slm3,  memory);
SurfStatView(slm3a.t.*mask,  avsurf,  'Memory on Volume t-values (df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Memory_Vol_T.png');
SurfStatView(SurfStatP(slm3a, mask),  avsurf,  'Memory on Volume p-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Memory_Vol_P.png');
SurfStatView(SurfStatQ(slm3a, mask),  avsurf,  'Memory on Volume q-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Memory_Vol_Q.png');

%VISUOSPATIAL%
slm4 = SurfStatLinMod(Vol,  M4,  avsurf);
slm4a = SurfStatT(slm4,  visuo);
SurfStatView(slm4a.t.*mask,  avsurf,  'Visuospatial ability on Volume t-values (df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Visuo_Vol_T.png');
SurfStatView(SurfStatP(slm4a, mask),  avsurf,  'Visuospatial ability on Volume p-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Visuo_Vol_P.png');
SurfStatView(SurfStatQ(slm4a, mask),  avsurf,  'Visuospatial ability on Volume q-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Visuo_Vol_Q.png');

%G%
slm5 = SurfStatLinMod(Vol,  M5,  avsurf);
slm5a = SurfStatT(slm5,  g);
SurfStatView(slm5a.t.*mask,  avsurf,  'g on Volume t-values (df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/g_Vol_T.png');
SurfStatView(SurfStatP(slm5a, mask),  avsurf,  'g on Volume p-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/g_Vol_P.png');
SurfStatView(SurfStatQ(slm5a, mask),  avsurf,  'g on Volume q-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/g_Vol_Q.png');

%OVERLAP & UNIQUE%
qval = SurfStatQ(slm1a, mask);
SpeedQM = 1*double(qval.Q < 0.05);
qval = SurfStatQ(slm2a, mask);
CrystQM = 2*double(qval.Q < 0.05);
qval = SurfStatQ(slm3a, mask);
MemoryQM = 3*double(qval.Q < 0.05);
qval = SurfStatQ(slm4a, mask);
VisuoQM = 4*double(qval.Q < 0.05);
qval = SurfStatQ(slm5a, mask);
GQM = 5*double(qval.Q < 0.05);

SurfStatView(((SpeedQM>0) + (CrystQM>0) + (MemoryQM>0) + (VisuoQM>0) + (GQM>0)).*mask, avsurf);
SurfStatColormap('hot')
sum(((SpeedQM>0) + (CrystQM>0) + (MemoryQM>0) + (VisuoQM>0) + (GQM>0))>0)
%209990/327684 vertices total

%% To get the overlapping areas:
SurfStatView((SpeedQM>0) + (CrystQM>0) + (MemoryQM>0) + (VisuoQM>0) + (GQM>0), avsurf);
SurfStatColormap('hot');
SurfStatColLim([0, 5]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/overlap.png');

%% To get ALL non-overlapping areas:
Speed_test = ((SpeedQM) & ~((CrystQM>0) | (MemoryQM>0) | (VisuoQM>0) | (GQM>0)))*1;
Cryst_test = ((CrystQM) & ~((SpeedQM>0) | (MemoryQM>0) | (VisuoQM>0) | (GQM>0)))*2;
Memory_test = ((MemoryQM) & ~((SpeedQM>0) | (CrystQM>0) | (VisuoQM >0) | (GQM>0)))*3;
Visuo_test = ((VisuoQM) & ~((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (GQM>0)))*4;
G_test = ((GQM) & ~((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (VisuoQM >0)))*5;
TOTAL_MASK = Speed_test + Cryst_test + Memory_test + Visuo_test + G_test;
SurfStatView(TOTAL_MASK, avsurf); 
SurfStatColormap(cmap)
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/unique.png'); 
sum(TOTAL_MASK >0)
%53790/327684 vertices

%% To get DOMAIN non-overlapping areas:
Speed_test = ((SpeedQM) & ~((CrystQM>0) | (MemoryQM>0) | (VisuoQM>0) | (GQM>0)))*1;
Cryst_test = ((CrystQM) & ~((SpeedQM>0) | (MemoryQM>0) | (VisuoQM>0) | (GQM>0)))*2;
Memory_test = ((MemoryQM) & ~((SpeedQM>0) | (CrystQM>0) | (VisuoQM >0) | (GQM>0)))*3;
Visuo_test = ((VisuoQM) & ~((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (GQM>0)))*4;
TOTAL_MASK = Speed_test + Cryst_test + Memory_test + Visuo_test;
SurfStatView(TOTAL_MASK, avsurf); 
SurfStatColormap(cmap)
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/unique_domains.png'); 

%% To get G ONLY 
G_test = ((GQM) & ~((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (VisuoQM >0)))*1
TOTAL_MASK = G_test
SurfStatView(TOTAL_MASK, avsurf)
SurfStatColormap(cmap)
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/unique_g.png');


%domain vertices (how many vertices for each domain)
sum(SpeedQM>0) 
%speed = 152085
sum(CrystQM>0) 
%cryst = 41274
sum(MemoryQM>0)
%memory = 39937
sum(VisuoQM>0) 
%visuo = 183014
sum(GQM>0) 
%g  = 177816

%domain vertices unique (how many vertices unique from g)
sum((SpeedQM) & ~(GQM>0)) 
% speed = 14068
sum((CrystQM) & ~(GQM>0)) 
% cryst = 162
sum((MemoryQM) & ~(GQM>0)) 
% memory = 642
sum((VisuoQM) & ~(GQM>0)) 
% visuo = 29794
sum((GQM) & ~((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (VisuoQM>0)))
% g = 5565


%% Map g vs non g domains vs overlap:
Speed_test = ((SpeedQM) & ~((CrystQM>0) | (MemoryQM>0) | (VisuoQM>0) | (GQM>0)))*1;
Cryst_test = ((CrystQM) & ~((SpeedQM>0) | (MemoryQM>0) | (VisuoQM>0) | (GQM>0)))*1;
Memory_test = ((MemoryQM) & ~((SpeedQM>0) | (CrystQM>0) | (VisuoQM >0) | (GQM>0)))*1;
Visuo_test = ((VisuoQM) & ~((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (GQM>0)))*1;
G_test = ((GQM) & ~((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (VisuoQM >0)))*2;
Overlap_test = ((GQM) & ((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (VisuoQM >0)))*3;
sum((GQM) & ((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (VisuoQM >0)))
TOTAL_MASK = Speed_test + Cryst_test + Memory_test + Visuo_test + G_test + Overlap_test;
SurfStatView(TOTAL_MASK, avsurf); 
SurfStatColormap(cmap)
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/domain_g_overlap.png'); 



%calculate overlap between domains only
sum((SpeedQM>0) & (CrystQM>0)) 
sum((SpeedQM>0) & (CrystQM>0) & ~(GQM>0)) 
% speed+cryst = 32769
% speed+cryst-g = 0
sum((SpeedQM>0) & (MemoryQM>0)) 
sum((SpeedQM>0) & (MemoryQM>0) & ~(GQM>0)) 
% speed+mem = 32275
% speed+mem-g = 0
sum((SpeedQM>0) & (VisuoQM>0)) 
sum((SpeedQM>0) & (VisuoQM>0) & ~(GQM>0)) 
% speed+visuo = 130027
% speed+visuo-g = 5104
sum((CrystQM>0) & (MemoryQM>0)) 
sum((CrystQM>0) & (MemoryQM>0) & ~(GQM>0)) 
% cryst+mem = 16834
% cryst+mem-g = 0
sum((CrystQM>0) & (VisuoQM>0)) 
sum((CrystQM>0) & (VisuoQM>0) & ~(GQM>0)) 
% cryst+visuo = 39231
% cryst+visuo-g = 0
sum((MemoryQM>0) & (VisuoQM>0)) 
sum((MemoryQM>0) & (VisuoQM>0) & ~(GQM>0)) 
% mem+visuo = 36550
% mem+visuo-g = 57 

%plot speed+visuo-g and mem+visuo-g
Speed_Visuo_test = ((SpeedQM>0) & (VisuoQM>0) & ~(GQM>0))*1;
Memory_Visuo_test = ((MemoryQM>0) & (VisuoQM>0) & ~(GQM>0)) *2;
TOTAL_MASK = Speed_Visuo_test + Memory_Visuo_test;
SurfStatView(TOTAL_MASK, avsurf); 
SurfStatColormap(cmap)
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/non_g_domains.png'); 


% SCATTERPLOT
tval1a = (slm1a.t.*mask)
tval2a = (slm2a.t.*mask)
tval3a = (slm3a.t.*mask)
tval4a = (slm4a.t.*mask)
tval5a = (slm5a.t.*mask)

% scatter plot with marker transparency and equal x,y axis
scatter(tval1a, tval2a, '.', 'MarkerEdgeColor', [0, 0.5, 1], 'MarkerEdgeAlpha', 0.25)
ax = axis; 
ax([1,3]) = min(ax([1,3]))
ax([2,4]) = max(ax([2,4]))
axis(ax); 
axis square;
xlabel('T-value 1');
ylabel('T-value 2');


% SIMULTANEOUS MAIN EFFECTS % SIMPLE VOLUME MODELS AGE, SEX, ICV CORRECTED
slm6 = SurfStatLinMod (Vol,  M6,  avsurf);

%PROCESSING SPEED%
slm6a = SurfStatT(slm6,  speed);
SurfStatView(slm6a.t.*mask,  avsurf,  'Processing speed on Volume t-values(df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Speed_Vol_T.png');
SurfStatView(SurfStatP(slm6a, mask),  avsurf,  'Processing speed on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Speed_Vol_P.png');
SurfStatView(SurfStatQ(slm6a, mask),  avsurf,  'Processing speed on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Speed_Vol_Q.png');

%CRYSTALLISED ABILITY%
slm6b = SurfStatT(slm6,  cryst);
SurfStatView(slm6b.t.*mask,  avsurf,  'Crystallised ability on Volume t-values(df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Cryst_Vol_T.png')
SurfStatView(SurfStatP(slm6b, mask),  avsurf,  'Crystallised ability on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Cryst_Vol_P.png')
SurfStatView(SurfStatQ(slm6b, mask),  avsurf,  'Crystallised ability on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Cryst_Vol_Q.png')

%MEMORY%
slm6c = SurfStatT(slm6,  memory);
SurfStatView(slm6c.t.*mask,  avsurf,  'Memory on Volume t-values(df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Memory_Vol_T.png')
SurfStatView(SurfStatP(slm6c, mask),  avsurf,  'Memory on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Memory_Vol_P.png')
SurfStatView(SurfStatQ(slm6c, mask),  avsurf,  'Memory on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Memory_Vol_Q.png')

%VISUOSPATIAL ABILITY%
slm6d = SurfStatT(slm6,  visuo);
SurfStatView(slm6d.t.*mask,  avsurf,  'Visuospatial ability on Volume t-values(df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Visuo_Vol_T.png')
SurfStatView(SurfStatP(slm6d, mask),  avsurf,  'Visuospatial ability on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Visuo_Vol_P.png')
SurfStatView(SurfStatQ(slm6d, mask),  avsurf,  'Visuospatial ability on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/Visuo_Vol_Q.png')

%G%
slm6e = SurfStatT(slm6,  g);
SurfStatView(slm6e.t.*mask,  avsurf,  'g on Volume t-values(df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/G_Vol_T.png')
SurfStatView(SurfStatP(slm6e, mask),  avsurf,  'g on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/G_Vol_P.png')
SurfStatView(SurfStatQ(slm6e, mask),  avsurf,  'g on Volume t-values(df = 617)');  
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/G_Vol_Q.png')

%OVERLAP & UNIQUE%
qval = SurfStatQ(slm6a, mask);
SpeedQM = 1*double(qval.Q < 0.05);
qval = SurfStatQ(slm6b, mask);
CrystQM = 2*double(qval.Q < 0.05);
qval = SurfStatQ(slm6c, mask);
MemoryQM = 3*double(qval.Q < 0.05);
qval = SurfStatQ(slm6d, mask);
VisuoQM = 4*double(qval.Q < 0.05);
qval = SurfStatQ(slm6e, mask);
GQM = 5*double(qval.Q < 0.05);

SurfStatView(((SpeedQM>0) + (CrystQM>0) + (MemoryQM>0) + (VisuoQM>0) + (GQM>0)).*mask, avsurf);
SurfStatColormap('hot')
sum(((SpeedQM>0) + (CrystQM>0) + (MemoryQM>0) + (VisuoQM>0) + (GQM>0))>0)
%0/327684 vertices

%% To get the overlapping areas:
SurfStatView((SpeedQM>0) + (CrystQM>0) + (MemoryQM>0) + (VisuoQM>0) + (GQM>0), avsurf);
SurfStatColormap('hot');
set(gcf,'PaperPosition',bbox);
SurfStatColLim([0, 5]);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/overlap.png');

%% To get the non-overlapping areas:
Speed_test = ((SpeedQM) & ~((CrystQM>0) | (MemoryQM>0) | (VisuoQM>0) | (GQM>0)))*1;
Cryst_test = ((CrystQM) & ~((SpeedQM>0) | (MemoryQM>0) | (VisuoQM>0) | (GQM>0)))*2;
Memory_test = ((MemoryQM) & ~((SpeedQM>0) | (CrystQM>0) | (VisuoQM >0) | (GQM>0)))*3;
Visuo_test = ((VisuoQM) & ~((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (GQM>0)))*4;
G_test = ((GQM) & ~((SpeedQM>0) | (CrystQM>0) | (MemoryQM>0) | (VisuoQM >0)))*5;
TOTAL_MASK = Speed_test + Cryst_test + Memory_test + Visuo_test + G_test;
SurfStatView(TOTAL_MASK, avsurf); 
SurfStatColormap(cmap)
set(gcf,'PaperPosition',bbox);
SurfStatColLim([0, 5]);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Simultaneous/unique.png'); 
sum(TOTAL_MASK >0)
%0/327684 vertices



% BIFACTOR MODELLING OF DOMAINS & G%%

% SET UP MODELS M7-M12
M7 = 1 + Agedays + Sex + Icv + Bi_speed;
M8 = 1 + Agedays + Sex + Icv + Bi_cryst;
M9 = 1 + Agedays + Sex + Icv + Bi_memory;
M10 = 1 + Agedays + Sex + Icv + Bi_visuo;
M11 = 1 + Agedays + Sex + Icv + Bi_g;
M12 = 1 + Agedays + Sex + Icv + Bi_speed + Bi_cryst + Bi_memory + Bi_visuo + Bi_g;

% SINGLE EFFECTS % SIMPLE VOLUME MODELS AGE, SEX, ICV CORRECTED
%PROCESSING SPEED%
slm7 = SurfStatLinMod(Vol,  M7,  avsurf);
slm7a = SurfStatT(slm7,  bi_speed);
SurfStatView(slm7a.t.*mask,  avsurf,  'Processing speed on Volume t-values (df = 617)');  % gives t maps
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_Speed_Vol_T.png');
SurfStatView(SurfStatP(slm7a, mask),  avsurf,  'Processing speed on Volume p-values (df = 617)');  % gives p maps & RFT clusters
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_Speed_Vol_P.png');
SurfStatView(SurfStatQ(slm7a, mask),  avsurf,  'Processing speed on Volume q-values (df = 617)');  % gives FDR q maps
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_Speed_Vol_Q.png');

%CRYSTALLISED%
slm8 = SurfStatLinMod(Vol,  M8,  avsurf);
slm8a = SurfStatT(slm8,  bi_cryst);
SurfStatView(slm8a.t.*mask,  avsurf,  'Crystallised ability on Volume t-values (df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_Cryst_Vol_T.png');
SurfStatView(SurfStatP(slm8a, mask),  avsurf,  'Crystallised ability on Volume p-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_Cryst_Vol_P.png');
SurfStatView(SurfStatQ(slm8a, mask),  avsurf,  'Crystallised ability on Volume q-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_Cryst_Vol_Q.png');

%MEMORY%
slm9 = SurfStatLinMod(Vol,  M9,  avsurf);
slm9a = SurfStatT(slm9,  bi_memory);
SurfStatView(slm9a.t.*mask,  avsurf,  'Memory on Volume t-values (df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_memory_Vol_T.png');
SurfStatView(SurfStatP(slm9a, mask),  avsurf,  'Memory on Volume p-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_memory_Vol_P.png');
SurfStatView(SurfStatQ(slm9a, mask),  avsurf,  'Memory on Volume q-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_memory_Vol_Q.png');

%VISUOSPATIAL%
slm10 = SurfStatLinMod(Vol,  M10,  avsurf);
slm10a = SurfStatT(slm10,  bi_visuo);
SurfStatView(slm10a.t.*mask,  avsurf,  'Visuospatial ability on Volume t-values (df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_Visuo_Vol_T.png');
SurfStatView(SurfStatP(slm10a, mask),  avsurf,  'Visuospatial ability on Volume p-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_Visuo_Vol_P.png');
SurfStatView(SurfStatQ(slm10a, mask),  avsurf,  'Visuospatial ability on Volume q-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_Visuo_Vol_Q.png');

%G%
slm11 = SurfStatLinMod(Vol,  M11,  avsurf);
slm11a = SurfStatT(slm11,  bi_g);
SurfStatView(slm11a.t.*mask,  avsurf,  'g on Volume t-values (df = 617)');  
SurfStatColormap(cmap);
SurfStatColLim([-7.1269, 7.1269]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_g_Vol_T.png');
SurfStatView(SurfStatP(slm11a, mask),  avsurf,  'g on Volume p-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_g_Vol_P.png');
SurfStatView(SurfStatQ(slm11a, mask),  avsurf,  'g on Volume q-values (df = 617)');  
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_g_Vol_Q.png');

%OVERLAP & UNIQUE%
qval = SurfStatQ(slm7a, mask);
Bi_speedQM = 1*double(qval.Q < 0.05);
qval = SurfStatQ(slm8a, mask);
Bi_crystQM = 2*double(qval.Q < 0.05);
qval = SurfStatQ(slm9a, mask);
Bi_memoryQM = 3*double(qval.Q < 0.05);
qval = SurfStatQ(slm10a, mask);
Bi_visuoQM = 4*double(qval.Q < 0.05);
qval = SurfStatQ(slm11a, mask);
Bi_gQM = 5*double(qval.Q < 0.05);

SurfStatView(((Bi_speedQM>0) + (Bi_crystQM>0) + (Bi_memoryQM>0) + (Bi_visuoQM>0) + (Bi_gQM>0)).*mask, avsurf);
SurfStatColormap('hot')
sum(((Bi_speedQM>0) + (Bi_crystQM>0) + (Bi_memoryQM>0) + (Bi_visuoQM>0) + (Bi_gQM>0))>0)
%180987/327684 vertices

%% To get the overlapping areas:
SurfStatView((Bi_speedQM>0) + (Bi_crystQM>0) + (Bi_memoryQM>0) + (Bi_visuoQM>0) + (Bi_gQM>0), avsurf);
SurfStatColormap('hot');
SurfStatColLim([0, 5]);
set(gcf,'PaperPosition',bbox);
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_overlap.png');

%% To get ALL non-overlapping areas:
Bi_speed_test = ((Bi_speedQM) & ~((Bi_crystQM>0) | (Bi_memoryQM>0) | (Bi_visuoQM>0) | (Bi_gQM>0)))*1;
Bi_cryst_test = ((Bi_crystQM) & ~((Bi_speedQM>0) | (Bi_memoryQM>0) | (Bi_visuoQM>0) | (Bi_gQM>0)))*2;
Bi_memory_test = ((Bi_memoryQM) & ~((Bi_speedQM>0) | (Bi_crystQM>0) | (Bi_visuoQM >0) | (Bi_gQM>0)))*3;
Bi_visuo_test = ((Bi_visuoQM) & ~((Bi_speedQM>0) | (Bi_crystQM>0) | (Bi_memoryQM>0) | (Bi_gQM>0)))*4;
Bi_g_test = ((Bi_gQM) & ~((Bi_speedQM>0) | (Bi_crystQM>0) | (Bi_memoryQM>0) | (Bi_visuoQM >0)))*5;
TOTAL_MASK = Bi_speed_test + Bi_cryst_test + Bi_memory_test + Bi_visuo_test + Bi_g_test;
SurfStatView(TOTAL_MASK, avsurf); 
SurfStatColormap(cmap)
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_unique.png')
sum(TOTAL_MASK >0)
%161537/327684 vertices

%% To get the DOMAIN non-overlapping areas:
Bi_speed_test = ((Bi_speedQM) & ~((Bi_crystQM>0) | (Bi_memoryQM>0) | (Bi_visuoQM>0) | (Bi_gQM>0)))*1;
Bi_cryst_test = ((Bi_crystQM) & ~((Bi_speedQM>0) | (Bi_memoryQM>0) | (Bi_visuoQM>0) | (Bi_gQM>0)))*2;
Bi_memory_test = ((Bi_memoryQM) & ~((Bi_speedQM>0) | (Bi_crystQM>0) | (Bi_visuoQM >0) | (Bi_gQM>0)))*3;
Bi_visuo_test = ((Bi_visuoQM) & ~((Bi_speedQM>0) | (Bi_crystQM>0) | (Bi_memoryQM>0) | (Bi_gQM>0)))*4;
TOTAL_MASK = Bi_speed_test + Bi_cryst_test + Bi_memory_test + Bi_visuo_test;
SurfStatView(TOTAL_MASK, avsurf); 
SurfStatColormap(cmap)
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_unique_domains.png')

%% To get G ONLY
Bi_g_test = ((Bi_gQM) & ~((Bi_speedQM>0) | (Bi_crystQM>0) | (Bi_memoryQM>0) | (Bi_visuoQM >0)))*1;
TOTAL_MASK = Bi_g_test;
SurfStatView(TOTAL_MASK, avsurf); 
SurfStatColormap(cmap)
set(gcf,'PaperPosition',bbox)
saveas(gcf, '/CCACE_Shared/Colin/Danielle_Cortex-g/NeuralArchImages/Individual/Bi_unique_g.png')

%domain vertices (how many vertices for each domain)
sum(Bi_speedQM>0) 
%speed = 20082
sum(Bi_crystQM>0) 
%cryst = 0
sum((Bi_memoryQM>0)>0) 
%memory = 0
sum((Bi_visuoQM>0)>0) 
%visuo = 3659
sum(Bi_gQM>0) 
%g = 177343

%domain vertices unique (how many vertices unique from g)
sum((Bi_speedQM) & ~(Bi_gQM>0)) 
% speed = 3175
sum((Bi_crystQM) & ~(Bi_gQM>0)) 
% cryst = 0
sum((Bi_memoryQM) & ~(Bi_gQM>0)) 
% memory = 0
sum((Bi_visuoQM) & ~(Bi_gQM>0)) 
% visuo = 469
sum((Bi_gQM) & ~((Bi_speedQM>0) | (Bi_crystQM>0) | (Bi_memoryQM>0) | (Bi_visuoQM >0)))
% g = 157893



%%CORRELATIONS BETWEEN T MAPS FOR G + DOMAINS
qval1a = SurfStatQ(slm1a,mask);
qval1a = 1*double(qval1a.Q < 0.05);
qval2a = SurfStatQ(slm2a,mask);
qval2a = 1*double(qval2a.Q < 0.05);
qval3a = SurfStatQ(slm3a,mask);
qval3a = 1*double(qval3a.Q < 0.05);
qval4a = SurfStatQ(slm4a,mask);
qval4a = 1*double(qval4a.Q < 0.05);
qval5a = SurfStatQ(slm5a,mask);
qval5a = 1*double(qval5a.Q < 0.05);

%%PROCESSING SPEED
%for whole cortex%
[ r p ] = corrcoef(slm1a.t.*mask,slm2a.t.*mask, 'rows', 'pairwise') 
[ r p ] = corrcoef(slm1a.t.*mask,slm3a.t.*mask, 'rows', 'pairwise') 
[ r p ] = corrcoef(slm1a.t.*mask,slm4a.t.*mask, 'rows', 'pairwise') 
[ r p ] = corrcoef(slm1a.t.*mask,slm5a.t.*mask, 'rows', 'pairwise') 
%for significant q vertices only%
[ r p ] = corrcoef(slm1a.t.*qval1a,slm2a.t.*qval2a, 'rows', 'pairwise')
[ r p ] = corrcoef(slm1a.t.*qval1a,slm3a.t.*qval3a, 'rows', 'pairwise')
[ r p ] = corrcoef(slm1a.t.*qval1a,slm4a.t.*qval4a, 'rows', 'pairwise')
[ r p ] = corrcoef(slm1a.t.*qval1a,slm5a.t.*qval5a, 'rows', 'pairwise')


%%CRYSTALLISED ABILITY
%for whole cortex%
[ r p ] = corrcoef(slm2a.t.*mask,slm3a.t.*mask, 'rows', 'pairwise')
[ r p ] = corrcoef(slm2a.t.*mask,slm4a.t.*mask, 'rows', 'pairwise')
[ r p ] = corrcoef(slm2a.t.*mask,slm5a.t.*mask, 'rows', 'pairwise')
%for significant q vertices only%
[ r p ] = corrcoef(slm2a.t.*qval2a,slm3a.t.*qval3a, 'rows', 'pairwise')
[ r p ] = corrcoef(slm2a.t.*qval2a,slm4a.t.*qval4a, 'rows', 'pairwise')
[ r p ] = corrcoef(slm2a.t.*qval2a,slm5a.t.*qval5a, 'rows', 'pairwise')


%%MEMORY
%for whole cortex%
[ r p ] = corrcoef(slm3a.t.*mask,slm4a.t.*mask, 'rows', 'pairwise')
[ r p ] = corrcoef(slm3a.t.*mask,slm5a.t.*mask, 'rows', 'pairwise')
%for significant q vertices only%
[ r p ] = corrcoef(slm3a.t.*qval3a,slm4a.t.*qval4a, 'rows', 'pairwise')
[ r p ] = corrcoef(slm3a.t.*qval3a,slm5a.t.*qval5a, 'rows', 'pairwise')

%%VISUOSPATIAL
%for whole cortex%
[ r p ] = corrcoef(slm4a.t.*mask,slm5a.t.*mask, 'rows', 'pairwise')
%for significant q vertices only%
[ r p ] = corrcoef(slm4a.t.*qval4a,slm5a.t.*qval5a, 'rows', 'pairwise')


speedt = slm1a.t.*mask;
crystt = slm2a.t.*mask;
memoryt = slm3a.t.*mask;
visuot  = slm4a.t.*mask;
gt = slm5a.t.*mask;
%scatter(speedt, gt);
%scatter(crystt, gt);
%scatter(memoryt, gt);
%scatter(visuot, gt);


% write vertex-wise CSV files for each domain, containing t-stats, p-values, and q-values
P1a = SurfStatP(slm1a, mask);
Q1a = SurfStatQ(slm1a, mask);
tab1a = array2table(vertcat(slm1a.t.*mask, P1a.P, Q1a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab1a, 'tables/Individual_Speed_Vol.csv', 'WriteVariableNames', true);

P2a = SurfStatP(slm2a, mask);
Q2a = SurfStatQ(slm2a, mask);
tab2a = array2table(vertcat(slm2a.t.*mask, P2a.P, Q2a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab2a, 'tables/Individual_Cryst_Vol.csv', 'WriteVariableNames', true);

P3a = SurfStatP(slm3a, mask);
Q3a = SurfStatQ(slm3a, mask);
tab3a = array2table(vertcat(slm3a.t.*mask, P3a.P, Q3a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab3a, 'tables/Individual_Memory_Vol.csv', 'WriteVariableNames', true);

P4a = SurfStatP(slm4a, mask);
Q4a = SurfStatQ(slm4a, mask);
tab4a = array2table(vertcat(slm4a.t.*mask, P4a.P, Q4a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab4a, 'tables/Individual_Visuo_Vol.csv', 'WriteVariableNames', true);

P5a = SurfStatP(slm5a, mask);
Q5a = SurfStatQ(slm5a, mask);
tab5a = array2table(vertcat(slm5a.t.*mask, P5a.P, Q5a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab5a, 'tables/Individual_g_Vol.csv', 'WriteVariableNames', true);

P6a = SurfStatP(slm6a, mask);
Q6a = SurfStatQ(slm6a, mask);
tab6a = array2table(vertcat(slm6a.t.*mask, P6a.P, Q6a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab6a, 'tables/Simultaneous_Speed_Vol.csv', 'WriteVariableNames', true);

P6b = SurfStatP(slm6b, mask);
Q6b = SurfStatQ(slm6b, mask);
tab6b = array2table(vertcat(slm6b.t.*mask, P6b.P, Q6b.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab6b, 'tables/Simultaneous_Cryst_Vol.csv', 'WriteVariableNames', true);

P6c = SurfStatP(slm6c, mask);
Q6c = SurfStatQ(slm6c, mask);
tab6c = array2table(vertcat(slm6c.t.*mask, P6c.P, Q6c.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab6c, 'tables/Simultaneous_Memory_Vol.csv', 'WriteVariableNames', true);

P6d = SurfStatP(slm6d, mask);
Q6d = SurfStatQ(slm6d, mask);
tab6d = array2table(vertcat(slm6d.t.*mask, P6d.P, Q6d.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab6d, 'tables/Simultaneous_Visuo_Vol.csv', 'WriteVariableNames', true);

P6e = SurfStatP(slm6e, mask);
Q6e = SurfStatQ(slm6e, mask);
tab6e = array2table(vertcat(slm6e.t.*mask, P6e.P, Q6e.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab6e, 'tables/Simultaneous_g_Vol.csv', 'WriteVariableNames', true);

P7a = SurfStatP(slm7a, mask);
Q7a = SurfStatQ(slm7a, mask);
tab7a = array2table(vertcat(slm7a.t.*mask, P7a.P, Q7a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab7a, 'tables/Individual_Bi_Speed_Vol.csv', 'WriteVariableNames', true);

P8a = SurfStatP(slm8a, mask);
Q8a = SurfStatQ(slm8a, mask);
tab8a = array2table(vertcat(slm8a.t.*mask, P8a.P, Q8a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab8a, 'tables/Individual_Bi_Cryst_Vol.csv', 'WriteVariableNames', true);

P9a = SurfStatP(slm9a, mask);
Q9a = SurfStatQ(slm9a, mask);
tab9a = array2table(vertcat(slm9a.t.*mask, P9a.P, Q9a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab9a, 'tables/Individual_Bi_Memory_Vol.csv', 'WriteVariableNames', true);

P10a = SurfStatP(slm10a, mask);
Q10a = SurfStatQ(slm10a, mask);
tab10a = array2table(vertcat(slm10a.t.*mask, P10a.P, Q10a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab10a, 'tables/Individual_Bi_Visuo_Vol.csv', 'WriteVariableNames', true);

P11a = SurfStatP(slm11a, mask);
Q11a = SurfStatQ(slm11a, mask);
tab11a = array2table(vertcat(slm11a.t.*mask, P11a.P, Q11a.Q)', 'VariableNames', {'t', 'p', 'q'});
writetable(tab11a, 'tables/Individual_Bi_g_Vol.csv', 'WriteVariableNames', true);
