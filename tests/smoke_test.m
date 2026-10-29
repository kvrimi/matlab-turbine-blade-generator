% 从仓库根目录运行：run('tests/smoke_test.m')
root_dir = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(root_dir,'blade_profile.m'));
assert(numel(x_mid) == 398);
assert(height(T) == 500);
assert(max([dist_Ca_Cp,dist_Cp_Cb,dist_Cb_Cq,dist_Cq_Ca]) < 1e-9);
assert(norm(xyz(1,:)-xyz(end,:)) < 1e-9);
assert(all(isfinite(xyz(:))));
assert(all(xyz(:,3) == 0));
for j = 1:numel(names)
    data = readmatrix(fullfile(out_dir,names{j}));
    assert(size(data,2) == 3);
    assert(max(abs(data(:)-curves{j}(:))) < 1e-8);
end
fprintf('PASS: default geometry, joins, closure, and exported XYZ files.\n');
