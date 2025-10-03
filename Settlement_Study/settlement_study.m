% moon_settlement_MC.m
% Monte Carlo settlement study for 20k people on the Moon
% Focuses ONLY on settlement distribution and clustering (no comms)
% Author: ChatGPT (template)

clear; close all; rng('shuffle');

%% Parameters
R_moon_km = 1737.4;         % mean lunar radius (km)
N_people = 20000;           % population per realization
numRuns = 10;               % number of Monte Carlo runs

% Settlement "centers" [lat_deg, lon_deg, weight, sigma_deg]
components = [
    90,    0,   0.30, 4;    % North pole
   -90,    0,   0.30, 4;    % South pole
    32,  -15,   0.15, 8;    % Mare Imbrium region
     0,  -60,   0.10, 10;   % Oceanus Procellarum
   -43,  -11,   0.15, 8     % Tycho crater region
];
components(:,3) = components(:,3) / sum(components(:,3));  % normalize weights

% DBSCAN clustering params (distance in km)
dbscan_eps_km = 50;   % cluster radius threshold
dbscan_minpts = 10;

%% Precompute component vectors
nComp = size(components,1);
comp_mu = zeros(3,nComp);
comp_sigma = components(:,4);
for k=1:nComp
    lat = deg2rad(components(k,1));
    lon = deg2rad(components(k,2));
    x = cos(lat)*cos(lon);
    y = cos(lat)*sin(lon);
    z = sin(lat);
    comp_mu(:,k) = [x; y; z];
end

%% Storage for results
results = struct('polarNorth',[],'polarSouth',[],'nClusters',[],'meanNN_km',[]);

for run=1:numRuns
    % Assign each person to a component
    % Weighted random sampling without toolbox
    edges = [0; cumsum(components(:,3))];
    r = rand(N_people,1);
    [~, idx] = histc(r, edges);
    
    % Sample positions
    pts = zeros(3,N_people);
    for k=1:nComp
        mask = (idx==k);
        nk = sum(mask);
        if nk==0, continue; end
        sigma_rad = deg2rad(comp_sigma(k));
        noise = randn(3,nk) * sigma_rad;
        samples = comp_mu(:,k) * ones(1,nk) + noise;
        samples = samples ./ vecnorm(samples); % normalize to sphere
        pts(:,mask) = samples;
    end
    
    % Convert to lunar coordinates
    xs = pts(1,:)*R_moon_km;
    ys = pts(2,:)*R_moon_km;
    zs = pts(3,:)*R_moon_km;
    lat_deg = rad2deg(asin(pts(3,:)));
    lon_deg = rad2deg(atan2(pts(2,:), pts(1,:)));
    
    % Polar populations
    north_count = sum(lat_deg >= 80);
    south_count = sum(lat_deg <= -80);
    
    % Clustering using DBSCAN (Euclidean = chord distance)
    coords_km = [xs' ys' zs'];
    labels = dbscan(coords_km, dbscan_eps_km, dbscan_minpts);
    nclusters = numel(unique(labels(labels>0)));
    
    % Mean nearest-neighbor distance (subset for speed)
    sampleIdx = randperm(N_people, min(2000,N_people));
    D = pdist2(coords_km(sampleIdx,:), coords_km(sampleIdx,:));
    D(1:size(D,1)+1:end) = inf;
    nn = min(D,[],2);
    arc_nn = 2*R_moon_km .* asin(min(1, nn./(2*R_moon_km)));
    meanNN_km = mean(arc_nn);
    
    % Store
    results(run).polarNorth = north_count;
    results(run).polarSouth = south_count;
    results(run).nClusters = nclusters;
    results(run).meanNN_km = meanNN_km;
    
    fprintf('Run %2d: North=%d South=%d Clusters=%d meanNN=%.1f km\n', ...
        run, north_count, south_count, nclusters, meanNN_km);
end

%% Summary
allNorth = [results.polarNorth];
allSouth = [results.polarSouth];
allClusters = [results.nClusters];
allNN = [results.meanNN_km];

fprintf('\nSummary over %d runs:\n', numRuns);
fprintf('  Avg North pole pop: %.1f ± %.1f\n', mean(allNorth), std(allNorth));
fprintf('  Avg South pole pop: %.1f ± %.1f\n', mean(allSouth), std(allSouth));
fprintf('  Avg clusters:       %.1f ± %.1f\n', mean(allClusters), std(allClusters));
fprintf('  Mean NN distance:   %.1f ± %.1f km\n', mean(allNN), std(allNN));

%% Visualization (last run)
figure('Name','Moon Settlements - 3D Scatter');
[xs_s, ys_s, zs_s] = sphere(120);
surf(R_moon_km*xs_s, R_moon_km*ys_s, R_moon_km*zs_s, ...
    'FaceAlpha',0.3,'EdgeColor','none'); colormap(gray); hold on;
scatter3(xs, ys, zs, 8, lat_deg, 'filled');
colorbar; title('Settlement points colored by latitude (deg)');
axis equal; lighting phong; camlight;

% Density heatmap
nLat = 180; nLon = 360;
lat_edges = linspace(-90,90,nLat);
lon_edges = linspace(-180,180,nLon);
H = histcounts2(lat_deg, lon_deg, lat_edges, lon_edges);
latc = (lat_edges(1:end-1)+lat_edges(2:end))/2;
lonc = (lon_edges(1:end-1)+lon_edges(2:end))/2;
[LonC,LatC] = meshgrid(lonc, latc);
Xg = R_moon_km * cosd(LatC).*cosd(LonC);
Yg = R_moon_km * cosd(LatC).*sind(LonC);
Zg = R_moon_km * sind(LatC);

figure('Name','Settlement Density Heatmap');
surf(Xg,Yg,Zg,H,'EdgeColor','none');
colormap(jet); colorbar; title('Population density heatmap');
axis equal; lighting phong; camlight;
view(30,30);
