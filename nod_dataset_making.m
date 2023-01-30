clear

skip = 3;   %スキップ幅(default 1)
dT  = 30;

path = uigetdir(pwd);
list = dir(fullfile(path, "*.csv"));
for num = 1:length(list)

%[file,path] = uigetfile('*.csv');   %choose movie
file = list(num).name;
csv_file = fullfile(path,file);
disp(csv_file);
in_file = csv_file;
confidences_threshold = 0;
% Where to store the output
output_dir = pwd;
% Most of the features will be in the csv file in the output directory with
% the same name as the input file
[~,name,~] = fileparts(in_file);
output_csv = sprintf('%s/%s.csv', output_dir, name);

% First read in teatures
tab = readtable(output_csv);
% particular f(output_csv);
column_names = tab.Properties.VariableNames;

% Read all of the data
all_params  = dlmread(output_csv, ',', 1, 0);

% This indicates which frames were succesfully tracked

% Find which column contains success of tracking data and timestamp data
valid_ind = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'success'));
time_stamp_ind = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'timestamp'));
%confidences
confidences_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'confidence'));
for j = 1:height(all_params)
    if(all_params(j, confidences_inds) >= confidences_threshold)
        all_params(j, confidences_inds) = true;
    else
        all_params(j, confidences_inds) = false;
    end

end

% Extract tracking success data and only read those frame
valid_cells = logical(all_params(:,valid_ind));
confidences_cells = logical(all_params(:, confidences_inds));

for i = 1:height(all_params)
    if (valid_cells(i, :) == true) && (confidences_cells(i, :) == true)
       all_params(i, valid_ind) = true;
    else
        all_params(i, valid_ind) = false;
    end
end

%confidences = all_params(valid_frames, confidences_inds);
valid_frames = logical(all_params(:,valid_ind));


% Get the timestamp data
time_stamps = all_params(valid_frames, time_stamp_ind);

% Finding which header line starts with p_ (basically model params)
shape_inds = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'p_'));

pose_inds = cellfun(@(x) ~isempty(x) && x == 1, strfind(column_names, 'pose_'));
pose_params = all_params(valid_frames, pose_inds);

% Output rigid (first 6) and non-rigid shape parameters
shape_params  = all_params(valid_frames, shape_inds);


% Demonstrate 3D landmarks
landmark_inds_x = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'X_'));
landmark_inds_y = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'Y_'));
landmark_inds_z = cellfun(@(x) ~isempty(x) && x==1, strfind(column_names, 'Z_'));

Xs = all_params(valid_frames, landmark_inds_x);
Ys = all_params(valid_frames, landmark_inds_y);
Zs = all_params(valid_frames, landmark_inds_z);


pose_tx = pose_params(:, 1);
pose_ty = pose_params(:, 2);
pose_tz = pose_params(:, 3);
pose_rx = pose_params(:, 4);
pose_ry = pose_params(:, 5);
pose_rz = pose_params(:, 6);

fpsper2 = 7;
coeff15flame = ones(1, fpsper2) / fpsper2;
pose_tx = filter(coeff15flame, 1, pose_tx);
pose_ty = filter(coeff15flame, 1, pose_ty);
pose_tz = filter(coeff15flame, 1, pose_tz);
pose_rx = filter(coeff15flame, 1, pose_rx);
pose_ry = filter(coeff15flame, 1, pose_ry);
pose_rz = filter(coeff15flame, 1, pose_rz);

time_stamps = time_stamps - time_stamps((fpsper2+1)/2, 1);

fps_h = (fpsper2+1)/2;
pose_tx = pose_tx(fps_h:end-fpsper2, 1);
pose_ty = (pose_ty(fps_h:end-fpsper2, 1));
pose_tz = (pose_tz(fps_h:end-fpsper2, 1));
pose_rx = (pose_rx(fps_h:end-fpsper2, 1));
pose_ry = (pose_ry(fps_h:end-fpsper2, 1));
pose_rz = (pose_rz(fps_h:end-fpsper2, 1));
time_stamps = time_stamps(fps_h:end-fpsper2, 1);

%length(time_stamps) >> length(time_stamps) - (fpsper2 - 1)

data_length = length(pose_rx);

nodtiming = readmatrix(strcat(name, "_nodtiming.xlsx"));

[row, ~] = size(nodtiming);
nod_bool = zeros(length(time_stamps), 1);
for u = 1:length(time_stamps)
    for s = 1:row
        if nodtiming(s, 1) < time_stamps(u, 1) && time_stamps(u, 1) < nodtiming(s, 2)
            %nod_bool(u-14, 1) = nodtiming(s, 1);
            nod_bool(u, 1) = 1;
        end
    end
end

countI = 0;
AngleNum = 2*dT/skip+1;
fm = (zeros(data_length-2*dT,(dT/skip+1)*3+2));
fm_A = zeros(data_length-2*dT,(2*dT/skip+1)*3);

w = (gausswin(AngleNum));
prev_d = 0;
% calculate relative eul angle
alpha_12 = zeros(1, 2*dT/skip+1);
beta_12 = zeros(1, 2*dT/skip+1);
gamma_12 = zeros(1, 2*dT/skip+1);

% gaussian window
W_alpha_12 = zeros(1, 2*dT/skip+1);
W_beta_12 = zeros(1, 2*dT/skip+1);
W_gamma_12 = zeros(1, 2*dT/skip+1);
for i = 1+dT:data_length-dT
    if(rem(i, 100) == 0)
        disp(i);
    end
    countI = countI + 1;
    %% calcurate relative rotate matrix
    % X_t >> X_cへの変換のため、角度はすべてマイナスにする
    alpha_1 = -pose_rx(i);
    beta_1 = -pose_ry(i);
    gamma_1 = -pose_rz(i);

    eul_1 = [gamma_1 beta_1 alpha_1];
    rot_1 = eul2rotm(eul_1);

    t_1 = [pose_tx(i); pose_ty(i); pose_tz(i)]; %縦ベクトル
    countJ = 0;
    for j = -dT:skip:dT
    % eul angle >> rotate matrix
    % X_t >> X_c (3*3)

    countJ = countJ + 1;

    alpha_2 = -pose_rx(i+j);
    beta_2 = -pose_ry(i+j);
    gamma_2 = -pose_rz(i+j);

    eul_2 = [gamma_2 beta_2 alpha_2];
    rot_2 = eul2rotm(eul_2);

    %calculate relative Rot and eul angle
    Rot_12 = rot_1\rot_2;
    eul_12 = rotm2eul(Rot_12);
    
    % calculate relative eul angle
    alpha_12(1, countJ) = eul_12(1);
    beta_12(1, countJ) = eul_12(2);
    gamma_12(1, countJ) = eul_12(3);

    % gaussian window
    W_alpha_12(1, countJ) = alpha_12(1, countJ) * w(countJ, 1);
    W_beta_12(1, countJ) = beta_12(1, countJ) * w(countJ, 1);
    W_gamma_12(1, countJ) = gamma_12(1, countJ) * w(countJ, 1);
    
    %% calculate Trans mat and D
    t_2 = [pose_tx(i+j); pose_ty(i+j); pose_tz(i+j)];
    Trans_12 =rot_1\(t_2 - t_1);
    Rot_12 = Rot_12;
    
    [OP,flag,relres,iter,resvec] = lsqr((1.0001 * eye(3) - Rot_12), Trans_12);
    [V, lambda] = eig(Rot_12);
    
    % find the position of eigenvalue 1
    [row, col] = find(abs(lambda - 1) < 0.0001);
     if ~isscalar(col)
        %d(1, s) = 0 ;
        %disp("以上あり")
        d(countJ, 1) = prev_d;
        continue;
     end

     u = V(:, col);% eigenvector of eigenvalue 1
     theta = acos(dot(u, OP) / (norm(u) * norm(OP)));
     d(countJ, 1) = norm(OP) * sin(theta);
     prev_d = d(countJ, 1);
     
    
    end

    %% calculate f_Rot
    Angle_12 = cat(2, alpha_12, beta_12, gamma_12);
    W_Angle_12 = cat(2, W_alpha_12, W_beta_12, W_gamma_12);

    FFT_alpha_12 = fft(W_alpha_12);
    FFT_beta_12 = fft(W_beta_12);
    FFT_gamma_12 = fft(W_gamma_12);

    FFT_angle_12 = cat(2, FFT_alpha_12, FFT_beta_12, FFT_gamma_12);

    norm_alpha_12 = arrayfun(@(x) norm(x), FFT_alpha_12);
    norm_beta_12 = arrayfun(@(x) norm(x), FFT_beta_12);
    norm_gamma_12 = arrayfun(@(x) norm(x), FFT_gamma_12);

    f_Rot = cat(2, norm_alpha_12(1:dT/skip+1), norm_beta_12(1:dT/skip+1), norm_gamma_12(1:dT/skip+1));
    
    %% calculate f_Axis(M)
    f_Axis = [mean(d) max(d)]/1000;
    f = cat(2, f_Rot, f_Axis);
    fm(countI, :) = f;
    fm_A(countI, :) = Angle_12;

end
fm = gather(fm);

%% dataset making
fm_time = time_stamps(dT+1:end-dT, 1);
nod_bool2 = nod_bool(dT+1:end-dT, 1);

max_d = max(fm(:, end));
fm(:, end-1:end) = fm(:, end-1:end)/max_d;

fmW = cat(2, fm_time, nod_bool2, fm);
fm_A2 = cat(2, fm_time, nod_bool2, fm_A);

dataset_name = strcat(name, "_dataset5", "_skip_", num2str(skip), ".csv");
angle_name = strcat(name, "_RelativeAngle", "_skip_", num2str(skip), ".csv");
writematrix(fmW, dataset_name);
writematrix(fm_A2, angle_name);

%writematrix(fm, fullfile(savefolder, strcat(name, "eur_angle.csv")));
DirName_D = strcat("datasetZYX","_skip_", num2str(skip));
DirName_A = strcat("AngleZYX","_skip_", num2str(skip));

mkdir(output_dir,DirName_D);
mkdir(output_dir, DirName_A);

movefile(dataset_name, fullfile(output_dir,DirName_D));
movefile(angle_name, fullfile(output_dir,DirName_A));

end

