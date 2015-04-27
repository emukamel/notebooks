clear;
datapath = 'G:/CSM/Jan22/Data';
addpath(datapath);
codepath = 'G:/CSM/Jan22/Code';
addpath(genpath(codepath));
filename = 'infile.txt';
fid = fopen(filename);
temp = textscan(fid, '%*s %d;%d;%d;%d; %d %d %d %d %f %f %s');
fclose(fid);
seg_pos = [temp{1}, temp{2}, temp{3}, temp{4}];
is_CSM_old = temp{8}';
seg_tmp = temp{11};
numseg = length(seg_tmp);
pval_vec = NaN(1, numseg);
cell_prop = NaN(1, numseg);

delta = 0.4;
nperm = 1000;

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
	[pval, gh] = CSMdetect_2step(X, [], delta, nperm);
    pval_vec(i) = pval;
    gh_uniq = unique(gh);
    cell_prop(i) = sum(gh == gh_uniq(1)) / length(gh);
end
toc(rtime);

fid = fopen('pval.txt', 'w');
fprintf(fid, '%f', pval_vec);
fclose(fid);
is_CSM_new = (pval_vec < 0.05);
