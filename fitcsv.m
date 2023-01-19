clear
flag = true;
countf = 0;

% dataset make
csvname = "";
while(flag)
    [file,  path] = uigetfile("*.csv");
    if(file == 0)
        break
    end

    % file name save
    csvname = strcat(csvname, file, "   ");
    fm = csvread(strcat(path, file), 0, 1);
    fm_Y = fm(:, 1);
    fm_X = fm(:, 2:end);
    if(countf == 0)
        X = fm_X;
        Y = fm_Y;
        countf = countf + 1;
        continue
    end

    X = cat(1, X, fm_X);
    Y = cat(1, Y, fm_Y);
    countf = countf + 1;

end
%{
extractX = X(:, 15:16);
extractx = cat(2, extractX, X(:, 31:32), X(:, 47:48));
%}


Y_index1 = find(Y == 1);
len_Y1 = length(Y_index1);

Y_index0 = find(Y == 0);
len_Y0 = length(Y_index0);
diff_lenY = len_Y0 - len_Y1;
learning_weight = 1; % default 1, 

for i = 0.01:0.01:1
    frame = len_Y1 + round(diff_lenY * (1-learning_weight));
    Y_index00 = randsample(Y_index0, frame);
    
    Y_index = cat(1, Y_index1, Y_index00);
    Y_index = sort(Y_index);
    
    Y_extracted = Y(Y_index, :);
    X_extracted = X(Y_index, :);
    
    % svm
    
    X_g =  gpuArray(X_extracted);
    Y_g =  gpuArray(Y_extracted);
    function_name = 'rbf';
    tic
    SVMModel = fitcsvm(X_g, Y_g,'Standardize',true,'KernelFunction',function_name,...
        'KernelScale','auto');
    time = toc;
    disp("SVM end");
    %CVSVMModel = crossval(SVMModel);
    %disp("Closs validation end");
    %classLoss = kfoldLoss(CVSVMModel);
    d = datestr(now, 'yyyy-mm-dd_HH-MM-ss');
    
    matfilename = strcat(d, "_", function_name, "_.mat");
    save(fullfile( "./SVMfile3", matfilename), "SVMModel");
    
    % 
    
    S1 = sprintf("datafile :%s\n", csvname);
    S2 = sprintf("original learning_rate(Label0 : Label1)   :%d :   %d\n", len_Y0, len_Y1);
    S3 = sprintf("learning_wight:   %f\n", learning_weight);
    S4 = sprintf("learning_rate((Label0 : Label1))  %d  :   %d\n", frame, len_Y1);
    S5 = sprintf("learning time:    %f [s]\n", time);

    rformat='_Result.txt';
    result_file=strcat(d,rformat);
    writematrix(S1,result_file,'WriteMode','append');
    writematrix(S2,result_file,'WriteMode','append');
    writematrix(S3,result_file,'WriteMode','append');
    writematrix(S4,result_file,'WriteMode','append');
    writematrix(S5,result_file,'WriteMode','append');
    mkdir("./SVMfile3");

    movefile(result_file, "./SVMfile3");

    learning_weight = 1 - i;
end