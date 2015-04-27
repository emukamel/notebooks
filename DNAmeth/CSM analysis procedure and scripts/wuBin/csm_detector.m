% this script can analyze simple format of inputs directly
% Last modified on May 11, 2014
clear;
% add path
datapath = pwd;
addpath(datapath);
codepath = '/home/mingansun/03.CSM_model/wuBin'
addpath(genpath(codepath));

% read file
%filename = [datapath, '/testData.500K.pval.distance.txt'];
filename = 'methylPattern.txt'
fid = fopen(filename);
temp = textscan(fid, '%s');
fclose(fid);
seg_tmp = temp{1};
numseg = length(seg_tmp);
%numseg = 10000;
unit = round(numseg / 100);
depth = NaN(1, numseg);
% numunipolar = zeros(1, numseg);
pval_vec = NaN(1, numseg);
dist_vec = NaN(1, numseg);

% set parameters
delta = 0.5;
tau = 0.5;
nperm = 10000;

% get p-values
rtime = tic;
for i = 1 : numseg
	tmp1 = regexp(seg_tmp(i), '[;:]', 'split');
	tmp2 = tmp1{1};
	tmp3 = tmp2(1 : (end - 1));
	l = length(tmp3);
	tmp4 = reshape(tmp3, [2, l / 2]);
	ix_valid = cellfun(@isempty, strfind(tmp4(1, :), '?'));
	datai = tmp4(:, ix_valid);
    X_tmp = reshape(str2num([datai{1, :}]'), 4, []);
    X_rep = cellfun(@str2num, datai(2, :));
    X = NaN(4, sum(X_rep));
    e = cumsum(X_rep);
    b = [1, e(1 : end - 1) + 1];
    for j = 1 : length(X_rep)
        X(:, b(j) : e(j)) = repmat(X_tmp(:, j), [1, X_rep(j)]);
    end
% 	if (sum(sum(X == 0)) == 0 || sum(sum(X == 1)) == 0)
% 		numunipolar(i) = 1;
% 	end
	depth(i) = size(X, 2);
	[pval, gh, d] = CSMdetect_2step(X, [], delta, tau, nperm);
    pval_vec(i) = pval;
    dist_vec(i) = d;
    if mod(i, unit) == 0
        disp([num2str(floor(i / unit)), '% finished!']);
    end
end
toc(rtime);


% output p-values
fid = fopen('pval.d5t5.txt', 'w');
fprintf(fid, '%f\n', pval_vec);
fclose(fid);

dist = fopen('dist.d5t5.txt', 'w');
fprintf(dist, '%f\n', dist_vec);
fclose(dist);

%ix_csm = find(pval_vec <= 0.05);
%ix_ncsm = setdiff(1 : numseg, ix_csm);
%dist_csm = dist_vec(ix_csm);
%dist_ncsm = dist_vec(ix_ncsm);
%[f1, x1] = ksdensity(dist_csm, 'function', 'pdf');
%[f2, x2] = ksdensity(dist_ncsm, 'function', 'pdf');
%figure;
%plot(x1, f1);
%hold on;
%plot(x2, f2, 'r');
