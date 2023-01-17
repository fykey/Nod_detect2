clear

[file,  path] = uigetfile("*.mat");
load(file)

[file,  path] = uigetfile("*.csv");

fm = csvread(strcat(path, file), 0, 0);
time = fm(:, 1);
fm_Y = fm(:, 2);
fm_X = fm(:, 3:end);
extractX = fm_X(:, 15:16);
%fm_X = cat(2, extractX, fm_X(:, 31:32), fm_X(:, 47:48));

Y_pre = predict(SVMModel, fm_X);
count_h0 = 0;
count_h1 = 0;
Y_pre2 = Y_pre;

for i = 16:length(Y_pre)-15
    sum = 0;
    
    if(Y_pre(i, 1) == 1)

        for j = -15:15
            sum = sum + Y_pre(i+j, 1);
        end
        if sum<10
            Y_pre2(i, 1) = 0;
            count_h0 = count_h0 + 1;
            continue

        end

    else
        for j = -15:15
            sum = sum + Y_pre(i+j, 1);
        end
        if sum>14
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

[~, name, ~] = fileparts(file);
data = cat(2, time, fm_Y, Y_pre2);
csvwrite(strcat(name, "_ectracted.csv"), data);