clear
addpath('~/Documents/MATLAB/fieldtrip');
ft_defaults;
addpath(genpath('~/Documents/MATLAB/iso2mesh'));
addpath(genpath('~/Documents/MATLAB/NIRFASTer'));

load('standard_bem.mat')
elec = ft_read_sens('standard_1005.elc');
all_elec = elec;
% source = ft_read_headshape('cortex_5124.surf.gii');
source = ft_read_headshape('cortex_8196.surf.gii');
% source = ft_read_headshape('cortex_20484.surf.gii');
% load('fsaverage5/sourcemodel.mat')

% load('2022-03-11_001/2022-03-11_001_probeInfo.mat');
load('Pilot/2022-03-28/2022-03-28_001/2022-03-28_001_probeInfo.mat');
probes = probeInfo.probes;
labels = {'Fp1', 'AF7', 'AF3', 'F1', 'F3', 'F5', 'F7', 'FT7', 'FC5', 'FC3', 'FC1', 'C1', 'C3', 'C5', 'T7', 'TP7', ...
    'CP5', 'CP3', 'CP1', 'P1', 'P3', 'P5', 'P7', 'P9', 'PO7', 'PO3', 'O1', 'Iz', 'Oz', 'POz', 'Pz', 'CPz', ...
    'Fpz', 'Fp2', 'AF6', 'AF4', 'AFz', 'Fz', 'F2', 'F4', 'F6', 'F8', 'FT8', 'FC6', 'FC4', 'FC2', 'FCz', 'Cz', ...
    'C2', 'C4', 'C6', 'T8', 'TP8', 'CP6', 'CP4', 'CP2', 'P2', 'P4', 'P6', 'P8', 'P10', 'PO8', 'PO4', 'O2'}';

headmodel = [];
headmodel.vol = vol;
headmodel.grid.pos = source.pos;
headmodel.grid.inside = 1:size(source.pos,1);

chanpos = zeros(length(labels), size(elec.chanpos, 2));
elecpos = zeros(length(labels), size(elec.elecpos, 2));
chantype = cell(length(labels), size(elec.chantype, 2));
chanunit = cell(length(labels), size(elec.chanunit, 2));
label = cell(length(labels), 1);
for i = 1:length(labels)
    tmp = strcmpi(labels(i), elec.label);
    chanpos(i, :) = elec.chanpos(tmp, :);
    elecpos(i, :) = elec.elecpos(tmp, :);
    chantype(i, :) = elec.chantype(tmp, :);
    chanunit(i, :) = elec.chanunit(tmp, :);
    label(i) = elec.label(tmp);
end
if ~all(strcmpi(label, labels))
    error('Label mismatch!');
end
elec.chanpos = chanpos;
elec.elecpos = elecpos;
elec.chantype = chantype;
elec.chanunit = chanunit;
elec.label = label;
% adjust the electrode locations
cfg = [];
cfg.method = 'project';
cfg.headshape = vol.bnd(1);
cfg.elec = elec;
elec = ft_electroderealign(cfg);
% Have a look. Note that the model is de-faced
% figure,ft_plot_mesh(vol.bnd(1))
% hold on
% scatter3(elec.chanpos(:,1),elec.chanpos(:,2),elec.chanpos(:,3),'r','filled')
% camlight
headmodel.elec = elec;

lf = ft_prepare_leadfield(headmodel);

triangulated = triangulation(source.tri, source.pos);
normals = vertexNormal(triangulated);
fprintf('Applying orientation constraints to dipoles...\n')
L = zeros(length(lf.label), length(lf.leadfield));
for i = 1:size(L, 2)
    L(:,i) =  lf.leadfield{i} * normals(i, :)';
end
if any(isnan(L(:)))
	warning('NaN values found in the leadfield matrix were set to zero');
	L(isnan(L)) = 0;
end

% save('templatebrain/eegmodel8196','headmodel','L','source');
save('fsaverage5/eegmodel_fs','headmodel','L','source','source_inflated')
%% NIRS starts here

% Mesh
% [node,elem,face]=s2m([vol.bnd(1).pos;vol.bnd(2).pos;vol.bnd(3).pos],[vol.bnd(1).tri;vol.bnd(2).tri+length(vol.bnd(1).pos);vol.bnd(3).tri+length(vol.bnd(1).pos)+length(vol.bnd(2).pos)],1,5,'tetgen',[80,-42,2;74,-35,1,;0,0,0]);
[node2,elem2]=surfboolean(vol.bnd(1).pos,vol.bnd(1).tri,'all',vol.bnd(2).pos,vol.bnd(2).tri,'all',vol.bnd(3).pos,vol.bnd(3).tri);   % fix intersecting problems
[node,elem,face]=s2m(node2,elem2,1,10,'tetgen',[80,0,0;68,0,0,;0,0,0]);
% Have a look
% figure; hold on;    % this is iso2mesh function titled plotmesh, but renamed to not be confused with the nirfast plotmesh
% iso2mesh_plotmesh(node,elem(elem(:,5)==3,:),'y>0','FaceColor',[0.35 0.35 0.35],'EdgeAlpha',0.6) %%brain
% iso2mesh_plotmesh(node,elem(elem(:,5)==2,:),'y>0','FaceColor',[1 1 0.9],'EdgeAlpha',0.6) %%bone
% iso2mesh_plotmesh(node,elem(elem(:,5)==1,:),'y>0','FaceColor',[1 0.8 0.7],'EdgeAlpha',0.6)%% scalp/muscle

mesh = [];
mesh.ele = elem;
mesh.node = node;
mesh.nnpe = 4;
mesh.dim = 3;
% solidmesh2nirfast(mesh,'/home/jiaming/Dropbox (CMU Biophotonics Lab)/CMU Biophotonics/Users/Jiaming/Research/VisualNIRX/templatebrain/nirs_mesh8196_full','stnd');
% nirs_mesh = load_mesh('templatebrain/nirs_mesh8196_full');
solidmesh2nirfast(mesh,'/home/jiaming/Dropbox (CMU Biophotonics Lab)/CMU Biophotonics/Users/Jiaming/Research/VisualNIRX/fsaverage5/nirs_mesh_fs_full','stnd');
nirs_mesh = load_mesh('fsaverage5/nirs_mesh_fs_full');
%%
% Place optodes
src_label = probes.labels_s;
det_label = probes.labels_d;
idx_src = zeros(length(src_label),1);
idx_det = zeros(length(det_label),1);
for i=1:length(src_label)
    idx_src(i) = find(strcmpi(src_label{i}, all_elec.label));
end
for i=1:length(det_label)
    idx_det(i) = find(strcmpi(det_label{i}, all_elec.label));
end

src = [];
src.chanpos = all_elec.elecpos(idx_src,:);
src.chantype = all_elec.chantype(idx_src,:);
src.chanunit = all_elec.chanunit(idx_src,:);
src.elecpos = all_elec.elecpos(idx_src,:);
src.label = all_elec.label(idx_src,:);
src.type = all_elec.type;
src.unit = all_elec.unit;
cfg = [];
cfg.method = 'project';
cfg.headshape = vol.bnd(1);
cfg.elec = src;
src = ft_electroderealign(cfg);

det = [];
det.chanpos = all_elec.elecpos(idx_det,:);
det.chantype = all_elec.chantype(idx_det,:);
det.chanunit = all_elec.chanunit(idx_det,:);
det.elecpos = all_elec.elecpos(idx_det,:);
det.label = all_elec.label(idx_det,:);
det.type = all_elec.type;
det.unit = all_elec.unit;
cfg = [];
cfg.method = 'project';
cfg.headshape = vol.bnd(1);
cfg.elec = det;
det = ft_electroderealign(cfg);

nirs_mesh.source.coord = src.chanpos;
nirs_mesh.meas.coord = det.chanpos;
nirs_mesh.link = probes.index_c;
nirs_mesh.link = [nirs_mesh.link, ones(length(nirs_mesh.link), 1)];
nirs_mesh.source.num = (1:size(nirs_mesh.source.coord,1))';
nirs_mesh.source.fwhm = zeros(size(nirs_mesh.source.coord,1),1);
nirs_mesh.source.fixed =0;
nirs_mesh.source.distributed =0;
nirs_mesh.meas.num = (1:size(nirs_mesh.meas.coord,1))';
nirs_mesh.meas.fixed =0;
nirs_mesh.ri = 1.4 * ones(size(nirs_mesh.ri));  % refractive index
% save_mesh(nirs_mesh,'/home/jiaming/Dropbox (CMU Biophotonics Lab)/CMU Biophotonics/Users/Jiaming/Research/VisualNIRX/templatebrain/nirs_mesh8196_full');
% nirs_mesh = load_mesh('/home/jiaming/Dropbox (CMU Biophotonics Lab)/CMU Biophotonics/Users/Jiaming/Research/VisualNIRX/templatebrain/nirs_mesh8196_full');   % Need to load again to move optodes to correct positions
save_mesh(nirs_mesh,'/home/jiaming/Dropbox (CMU Biophotonics Lab)/CMU Biophotonics/Users/Jiaming/Research/VisualNIRX/fsaverage5/nirs_mesh_fs_full');
nirs_mesh = load_mesh('/home/jiaming/Dropbox (CMU Biophotonics Lab)/CMU Biophotonics/Users/Jiaming/Research/VisualNIRX/fsaverage5/nirs_mesh_fs_full');

% Get standard Jacobian first
% use 750nm numbers for 760 for now, data from Eggebrecht et al. 2012
fprintf('Calculating Jacobian for 760nm...\n')
nirs_mesh.mua(nirs_mesh.region == 1) = 0.0170;      % scalp; mm^-1
nirs_mesh.mua(nirs_mesh.region == 2) = 0.0116;      % skull
nirs_mesh.mua(nirs_mesh.region == 3) = 0.0180;      % brain; assume all gray matter

nirs_mesh.mus(nirs_mesh.region == 1) = 0.74;        % mm^-1; mus_prime
nirs_mesh.mus(nirs_mesh.region == 2) = 0.94;
nirs_mesh.mus(nirs_mesh.region == 3) = 0.8359;
J_760 = jacobian_FD(nirs_mesh, 0);

% 850nm, data from Eggevrecht et al. 2012
fprintf('Calculating Jacobian for 850nm...\n')
nirs_mesh.mua(nirs_mesh.region == 1) = 0.0190;      % mm^-1
nirs_mesh.mua(nirs_mesh.region == 2) = 0.0139;
nirs_mesh.mua(nirs_mesh.region == 3) = 0.0192;

nirs_mesh.mus(nirs_mesh.region == 1) = 0.64;        % mm^-1; mus_prime
nirs_mesh.mus(nirs_mesh.region == 2) = 0.84;
nirs_mesh.mus(nirs_mesh.region == 3) = 0.6726;
J_850 = jacobian_FD(nirs_mesh, 0);

% Now get the spectral Jacobian; origional is "del OD/del mua", convert to "del OD/del HbO (HbD)"
% [HbO_lambda1 HbD_lambda1]
% [HbO_lambda2 HbD_lambda2]
excoeff = importdata('excoef.txt');
excoeff = excoeff.data;
epsilon = [excoeff(find(excoeff(:,1)==760), 2:3); excoeff(find(excoeff(:,1)==850), 2:3)];
jacobian_full = [J_760.complete*epsilon(1,1), J_760.complete*epsilon(1,2); J_850.complete*epsilon(2,1), J_850.complete*epsilon(2,2)];

% jacobian_nirs = zeros(size(jacobian_full,1), size(source.pos,1)*2);
% for i=1:size(jacobian_nirs,1)
%     tmp1 = griddata(nirs_mesh.nodes(:,1),nirs_mesh.nodes(:,2),nirs_mesh.nodes(:,3),jacobian_full(i,1:length(node)),source.pos(:,1),source.pos(:,2),source.pos(:,3));
%     tmp2 = griddata(nirs_mesh.nodes(:,1),nirs_mesh.nodes(:,2),nirs_mesh.nodes(:,3),jacobian_full(i,length(node)+1:end),source.pos(:,1),source.pos(:,2),source.pos(:,3));
%     jacobian_nirs(i,:) = [tmp1;tmp2]';
% end
save('templatebrain/nirsmodel8196_full.mat', 'nirs_mesh', 'J_760', 'J_850', 'jacobian_full', '-v7.3');
% save('fsaverage5/nirsmodel_fs_full.mat', 'nirs_mesh', 'J_760', 'J_850', 'jacobian_full', '-v7.3');
