clear
[csvfile,  path] = uigetfile("*.csv");
dT = 30;

fm = readmatrix(fullfile(path, csvfile));
time = fm(:, 1);
fm_Y = fm(:, 2);
fm_X = fm(:, 3:end);
path = uigetdir(pwd);

list = dir(fullfile(path, "*.mat"));

fm_Y_d = circshift(fm_Y, 1);
span_v = fm_Y - fm_Y_d;
index_down_v = find(span_v == -1);
index_up_v = find(span_v == 1);
span_v2 = cat(2, index_up_v, index_down_v);

% SVM保存先のパスに移動
pre_d = pwd;
cd(path)

 nodAcc_frame = zeros(length(list), 1);
 nodAcc_timing = zeros(length(list), 1);
 countK = zeros(length(list), 1);
 

for k = 1:length(list)
    filename = list(k).name;
    filepath = fullfile(path, filename);
    load(filepath);

    Y_pre = predict(SVMModel, fm_X);
    count_h0 = 0;
    count_h1 = 0;
    Y_pre2 = Y_pre;

    for i = dT+1:length(Y_pre)-dT
        sum = 0;
        
        if(Y_pre(i,1) == 1)
    
            for j = -dT:dT
                sum = sum + Y_pre(i+j, 1);
            end
            if sum<10
                Y_pre2(i, 1) = 0;
                count_h0 = count_h0 + 1;
                continue
    
            end
    
        else
            for j = -dT:dT
                sum = sum + Y_pre(i+j, 1);
            end
            if sum>50
                Y_pre2(i, 1) = 1;
                count_h1 = count_h1 + 1;
                continue
    
            end
    
        end
        
    end

    count1 = 0;
    count2 = 0;
    count3 = 0;
    count4 = 0;
    
    for i = 1: length(fm_Y)
        if (fm_Y(i, 1) == 0 && (fm_Y(i, 1) == Y_pre2(i, 1)))
            count1 = count1 + 1;
           
            
        elseif fm_Y(i, 1) == 0 && (fm_Y(i, 1) ~= Y_pre2(i, 1))
            count2 = count2 + 1;
            
        elseif fm_Y(i, 1) == 1 && (fm_Y(i, 1) ~= Y_pre2(i, 1))
            count3 = count3 + 1;
            
        elseif (fm_Y(i, 1) == 1 && (fm_Y(i, 1) == Y_pre2(i, 1)))
            count4 = count4 + 1;
            
        end
    end
    nodAcc_frame(k, 1) = (count1 + count4) / (count1 + count2 + count3 + count4);

    Y_pred2_d = circshift(Y_pre2, 1);
    span_n = Y_pre2 - Y_pred2_d;
    index_down_n = find(span_n == -1);
    index_up_n = find(span_n == 1);

    span_n2 = cat(2, index_up_n, index_down_n);
    nodtiming_n_c = round((span_n2(:, 1) + span_n2(:, 2))/2);
    count_nod = 0;

    for n = 1:length(nodtiming_n_c)
        for m = 1:length(span_v2)
            if(span_v2(m, 1) <= nodtiming_n_c(n, 1) && nodtiming_n_c(n, 1) <= span_v2(m, 2))
                count_nod = count_nod + 1;
                continue;
            end

        end
    end
    countK(k, 1) = k;

    nodAcc_timing(k, 1) = count_nod/length(span_v2);

    filearray(k, 1) = string(filename);


    disp(k);


end
alldata = num2cell(cat(2, countK, nodAcc_frame, nodAcc_timing));

A = cat(2, filearray, alldata);

writematrix(A, fullfile(path, "Accracy.csv"));