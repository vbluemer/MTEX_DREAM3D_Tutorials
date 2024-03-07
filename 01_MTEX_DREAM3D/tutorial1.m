
%% 
% In this tutorial, we are loading a ctf file, cropping the observed region, 
% plotting the ebsd, reconstructing grains and fitting some 
% functions to the grain properties

close all
show_plots = 0;
if not(show_plots)
    disp('Plotting is turned off')
end


%%
% First of all, add the MTEX folder to the path
mtex_pth = 'mtex-5.9.0';

s       = pathsep;
pathStr = [s, path, s];
onPath  = contains(pathStr, [s, mtex_pth, s], 'IgnoreCase', ispc);

if not(onPath)
    addpath(mtex_pth)
    startup_mtex
end
%% Specify Crystal and Specimen Symmetries

% this set of crystal symmetries is generated from the EBSD import wizard
% which is started by running import_wizard('EBSD') in the command line.
% Special care needs to be taken with the options regarding Sample and
% Acquisition coordinate systems. Euler angles are recorded w.r.t.
% intrinsic (updated) ZXZ axes of the Sample CS, spatial x,y coordinates
% are recorded in the Acquisition CS. In our data, both CS are aligned so
% a rotation of [0,0,0] between the two was specified in the wizard.
CS = {... 
  'notIndexed',...
  'notIndexed',...
  'notIndexed',...
  'notIndexed',...
  'notIndexed',...
  crystalSymmetry('6/mmm', [3 3 4.7], 'X||a*', 'Y||b', 'Z||c*', 'mineral', 'Ti-Hex', 'color', [0 0 0.55]),...
  crystalSymmetry('m-3m', [3.2 3.2 3.2], 'mineral', 'Titanium cubic', 'color', [0 0.39 0])};

% Plotting convention. This only affects how your data appears on screen
% (which way is up). The relation of Sample and Acquisition CS is
% unaffected 
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');
%% Specify File Names and import
disp('Load EBSD ...')

% EBSD is loaded from file
pname = 'ebsd_inputs/';
fname = 'CTF_File_68.ctf';
pth = [pname fname];

ebsd = EBSD.load(pth,CS,'interface','ctf');
%%
% A couple of things you can do with the newly created ebsd variable
if show_plots    
    ebsd(1:10)
    ebsd(1:10).orientations
    ebsd('Titanium cubic')
    ebsd_cubic=ebsd('Titanium cubic');
end
%% Plotting
if show_plots    
    figure
    plot(ebsd);
    figure
    plot(ebsd('Ti-Hex'),ebsd('Ti-Hex').orientations);
end
%% Focus on Ti-Hex phase from now on 
ebsd=ebsd('Ti-Hex');

disp('Cropping ...')

% Specify region to be cropped
% Origin x, Origin y, Length x, Length y
region = [5 2 30 50];

% Show region
if show_plots    
    plot(ebsd,ebsd.orientations);
    rectangle('position',region,'edgecolor','r','linewidth',2);
end
% Apply cropping
condition = inpolygon(ebsd,region);
ebsd = ebsd(condition);

% Plot new ebsd
if show_plots    
    plot(ebsd,ebsd.orientations)
end


%% Reconstruct grains
% So far the ebsd is only pixel data, so we cannot compute any statistics
% over the grains or convert to a DAMASK domain 

disp('Reconstructing grains ...')

% Grain misorientation threshold in degrees
alpha=5;
% Grains are reconstructed. New variable grains is created. ebsd gets new
% attribute "grainId", which records the association of each pixel to a
% grain 
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',alpha*degree);
plot(grains,grains.meanOrientation)

% Some examples of what you can do with the new variables
if show_plots    
    ebsd(1:10).grainId
    grains(1:10).aspectRatio
    grains(1:10).area
    grains(1:10).meanOrientation
    plot(grains(1:70),grains(1:70).meanOrientation);
end
%% Extract some statistics from the grains

disp('Extracting statistics ...')

% MTEX works with equivalent radius, Dream3D with equivalend diameter
eq_diameter = 2*grains.equivalentRadius;
% fit a lognormal distribution which Dream3D expects as input. It is
% characterized by mu and sigma
pd   = fitdist(eq_diameter,'Lognormal');

% The components can be accessed like so
pd.mu
pd.sigma

% Plot the distribution
f=figure;
p=plot(pd);

% Some graphical adjustments
p(1).LineWidth = 2
xlabel('Equivalent Sphere Diameter [μm]')
ylabel('Probability Distribution [-]')
f.Position = [100 100 900 600];
ylim([0,3.5])
xlim([0,4])
%% 
% Small grains are quite dominant by number (less so by area of course).
% Here is how to remove them if needed.
ebsd_backup = ebsd;

disp('Remove small grains ...')


% Pixels that belong to small grains are removed from the EBSD
ebsd(grains(grains.equivalentRadius<=0.2))   = [];

% Grains need to be re-calculated. Grains variable is overwritten
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd,'angle',alpha*degree);

% Plot new grains
if show_plots    
    plot(grains,grains.meanOrientation)
end


% Redo the computation of statistics
eq_diameter = 2*grains.equivalentRadius;
pd_new   = fitdist(eq_diameter,'Lognormal');

% Plot figure and graphical adjustments
f=figure;
p=plot(pd_new);
p(1).LineWidth = 2;
xlabel('Equivalent Sphere Diameter [μm]')
ylabel('Probability Distribution [-]')
ylim([0,3.5])
xlim([0,4])
f.Position = [100 100 900 600];

%% Aspect Ratio

% A similar computation can be done for aspect ratios. Dream3D expects a
% fitted beta-distribution in the dialig window which is characterized by
% parameters a and b  

% take the inverse to fit Dream3D syntax
aspect_ratios = 1./grains.aspectRatio;

%fit distribution
pd_aspect = fitdist(aspect_ratios,'beta');

% individual components can be accessed like so
pd_aspect.a
pd_aspect.b

% plot figure
f=figure;
p=plot(pd_aspect);
p(1).LineWidth = 2
xlabel('Aspect Ratio [-]')
ylabel('Probability Distribution [-]')
ylim([0,7])
xlim([0,1])
f.Position = [100 100 900 600];

%%
disp('Export in Damask syntax')

% Damask works with different data style than EBSDs, as in, data needs to
% be provided in a grid as supposed to a list as is the case in EBSDs
ebsd = ebsd_backup;
% ebsd is "gridified"
ebsd_grid = gridify(ebsd);
% empty space is filled by nearest neighbor
ebsd_grid = fill(ebsd_grid);

[grains_grid,ebsd_grid.grainId,ebsd_grid.mis2mean] = calcGrains(ebsd_grid,'angle',alpha*degree);
loop_condition = any(isnan(grains_grid.meanOrientation));

while loop_condition
    disp('NaN found in grain orientations, truncating EBSD by one layer')
    ebsd_grid=ebsd_grid([2:end-1],[2:end-1]);
    [grains_grid,ebsd_grid.grainId,ebsd_grid.mis2mean] = calcGrains(ebsd_grid,'angle',alpha*degree);
    loop_condition = any(isnan(grains_grid.meanOrientation));
end

disp('No more NaNs in grain orientations')
% grains are calculated over the updated ebsd_grid to make sure there are
% no NaNs in the data

% Three files are written:
% (1) grid_grain_IDs.txt contains the grid of all points that make up the
%     simulation domain. At each gridpoint, the ID of the grain that it 
%     belongs to, is saved. 
% (2) grain_euler_angles.txt contains the a list of the rotations of the
%     grains. The index in the list coresponds to the grain ID
% (3) dimensions.txt contains the physical length dimensions of the
%     simulation domain.


filename = 'matlab_outputs\grid_grain_IDs.txt';
% to fit Damask convention:
%   axes are transposed 
%   1 is subtracted from index
grain_Id = transpose(ebsd_grid.grainId-1);
writematrix(grain_Id,filename)

filename = 'matlab_outputs\grain_euler_angles.txt';
eulers = rad2deg([grains_grid.meanOrientation.phi1 grains_grid.meanOrientation.Phi grains_grid.meanOrientation.phi2]);
writematrix(round(eulers,2),filename,Delimiter=',')

filename = 'matlab_outputs\dimensions.txt';
writematrix(round([ebsd_grid.xmin ebsd_grid.xmax ebsd_grid.ymin ebsd_grid.ymax ebsd_grid.dx ebsd_grid.dy],4),filename,Delimiter=',')

disp('Finished successfully')
