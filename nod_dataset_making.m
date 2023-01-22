clear

skip = 3;   %スキップ幅(default 1)
dT  = 30;

[file,path] = uigetfile('*.csv');   %choose movie

f = 4.0;%焦点距離 https://support.logi.com/hc/ja/articles/360053005413-C505-StreamCam-Technical-Specifications
csv_file = fullfile(path,file);
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
pose_ty = pose_ty(fps_h:end-fpsper2, 1);
pose_tz = pose_tz(fps_h:end-fpsper2, 1);
pose_rx = pose_rx(fps_h:end-fpsper2, 1);
pose_ry = pose_ry(fps_h:end-fpsper2, 1);
pose_rz = pose_rz(fps_h:end-fpsper2, 1);
time_stamps = time_stamps(fps_h:end-fpsper2, 1);

%length(time_stamps) >> length(time_stamps) - (fpsper2 - 1)

frame_sum = sum(logical(time_stamps));

nodtiming = readmatrix(strcat(name, "_nodtiming.xlsx"));

[row, ~] = size(nodtiming);
nod_bool = zeros(length(time_stamps), 1);
for u = 1:length(time_stamps)
    for s = 1:row
        if nodtiming(s, 1) < time_stamps(u, 1) && time_stamps(u, 1) < nodtiming(s, 2)
            %nod_bool(u-14, 1) = nodtiming(s, 1);
            nod_bool(u, 1) = 1;
            continue;
        end
    end
end



Rt1t2 = zeros(3,3,length(time_stamps)-skip);     %relative rotate(3*3*num of using data)
Tt1t2 = zeros(3,1,length(time_stamps)-skip);     %relative translate(3*3*num of using data)
eul_3 = zeros(length(time_stamps)-skip, 3);

fm = zeros(1,(dT+1)*3+2);
count = 1;
for i = 1:length(time_stamps)-skip
    % calcurate relative rotate and translate matrix
    % X_t >> X_cへの変換のため、角度はすべてマイナスにする
    alpha_1 = -pose_rx(i);
    alpha_2 = -pose_rx(i + skip);

    beta_1 = -pose_ry(i);
    beta_2 = -pose_ry(i + skip);

    gamma_1 = -pose_rz(i);
    gamma_2 = -pose_rz(i + skip);

    eul_1 = [gamma_1 beta_1 alpha_1];
    eul_2 = [gamma_2 beta_2 alpha_2];
    % eul angle >> rotate matrix
    % X_t >> X_c (3*3)
    rot_1 = eul2rotm(eul_1);
    rot_2 = eul2rotm(eul_2);

    T_01 = [pose_tx(i); pose_ty(i); pose_tz(i)];
    T_02 = [pose_tx(i+skip); pose_ty(i+skip); pose_tz(i+skip)];

    I = [1 0 0; 0 -1 0; 0 0 -1];

    T_w = [0; 0; 1000];
    %X_f >>( X_c )>> X_w 
    R_1 = I*rot_1;%3*3
    T_1 = I*T_01 + T_w;%3*1

    R_2 = I*rot_2;
    T_2 = I*T_02 + T_w;
    % relative rotate and translate matrix
    Rt1t2(:,:,count) = R_1\R_2;
    Tt1t2(:,:,count) =  R_1\(T_2 - T_1);
    
    %% calcurate eul angle from relative matrix(for Rotation axis features)
    eul_3(count, :) = rotm2eul(Rt1t2(1:3, 1:3, count));
    
    %% calcurate OP and d(for Head rotattion frequency features)
    OP = lsqr((1.0001 * eye(3) - Rt1t2(:,:,count)), Tt1t2(:,:,count));
    
    [V, lambda] = eig(Rt1t2(:,:,count));
    % find the position of eigenvalue 1
    [row, col] = find(abs(lambda - 1) < 0.0001);
    if ~isscalar(col)
        %d(1, s) = 0 ;
        disp("以上あり")
        d(count, 1) = prev_d;
        count = count + 1;

        if i > 2*dT + skip  %
        index_c = i - (dT+skip);  %center position
        angle_z = zeros(1, dT*2+1);
        angle_y = zeros(1, dT*2+1);
        angle_x = zeros(1, dT*2+1);
        d_list = zeros(1, dT*2+1);
        w = gausswin(2*dT+1);

        for j = -dT:dT
            index_n = index_c + j; % now position
            angle = eul_3(index_n, :);
            angle_z(1, j+dT+1) = w(j+dT+1)*angle(1, 1);
            angle_y(1, j+dT+1) = w(j+dT+1)*angle(1, 2);
            angle_x(1, j+dT+1) = w(j+dT+1)*angle(1, 3);

            d_list(1, j+dT+1) = d(index_n, 1);
        end

        fft_angle_z = fft(angle_z);
        fft_angle_y = fft(angle_y);
        fft_angle_x = fft(angle_x);

        norm_angle_z = arrayfun(@(x) norm(x), fft_angle_z);
        norm_angle_y = arrayfun(@(x) norm(x), fft_angle_y);
        norm_angle_x = arrayfun(@(x) norm(x), fft_angle_x);
%{
        fm_rot = cat(2, norm_angle_x(1, dT+1:2*dT+1), ...
            norm_angle_y(1, dT+1:2*dT+1), norm_angle_z(1, dT+1:2*dT+1));
%}      
        % fft-> [X(0), X(df), X(2df), ... , X(dT*df), X(-(dT-1)*df), ..., X(-df)]
        fm_rot = cat(2, norm_angle_x(1, 1:dT+1), ...
            norm_angle_y(1,1:dT+1), norm_angle_z(1, 1:dT+1));

        fm_axis = cat(2, mean(d_list), max(d_list));
        fm = cat(1, fm, cat(2, fm_rot, fm_axis));
    end

        continue
    end

    u = V(:, col);% eigenvector of eigenvalue 1
    theta = acos(dot(u, OP) / (norm(u) * norm(OP)));

    d(count, 1) = norm(OP) * sin(theta);
    prev_d = d(count, 1);
    count = count + 1;
    %% process of making 2 features
    if i > 2*dT + skip  %
        index_c = i - (dT+skip);  %center position
        angle_z = zeros(1, dT*2+1);
        angle_y = zeros(1, dT*2+1);
        angle_x = zeros(1, dT*2+1);
        d_list = zeros(1, dT*2+1);
        w = gausswin(2*dT+1);

        for j = -dT:dT
            index_n = index_c + j; % now position
            angle = eul_3(index_n, :);
            angle_z(1, j+dT+1) = w(j+dT+1)*angle(1, 1);
            angle_y(1, j+dT+1) = w(j+dT+1)*angle(1, 2);
            angle_x(1, j+dT+1) = w(j+dT+1)*angle(1, 3);

            d_list(1, j+dT+1) = d(index_n, 1);
        end

        fft_angle_z = fft(angle_z);
        fft_angle_y = fft(angle_y);
        fft_angle_x = fft(angle_x);

        norm_angle_z = arrayfun(@(x) norm(x), fft_angle_z);
        norm_angle_y = arrayfun(@(x) norm(x), fft_angle_y);
        norm_angle_x = arrayfun(@(x) norm(x), fft_angle_x);
%{
        fm_rot = cat(2, norm_angle_x(1, dT+1:2*dT+1), ...
            norm_angle_y(1, dT+1:2*dT+1), norm_angle_z(1, dT+1:2*dT+1));
%}      
        % fft-> [X(0), X(df), X(2df), ... , X(dT*df), X(-(dT-1)*df), ..., X(-df)]
        fm_rot = cat(2, norm_angle_x(1, 1:dT+1), ...
            norm_angle_y(1,1:dT+1), norm_angle_z(1, 1:dT+1));

        fm_axis = cat(2, mean(d_list), max(d_list));
        fm = cat(1, fm, cat(2, fm_rot, fm_axis));
    end
end

%% dataset making
time_stamps(1:dT, :) = [];
time_stamps(end+1-(dT+2*skip):end, :) = [];

nod_bool(1:dT, :) = [];
nod_bool(end+1-(dT+2*skip):end, :) = [];

fm(1, :) = [];

max_d = max(fm(:, end));
fm(:, end-1:end) = fm(:, end-1:end)/max_d;

fm = cat(2, time_stamps, nod_bool, fm);

dataset_name = strcat(name, "_dataset3.csv");
writematrix(fm, dataset_name);

mkdir(output_dir, "dataset")
movefile(dataset_name, strcat(output_dir, "/dataset/"))




