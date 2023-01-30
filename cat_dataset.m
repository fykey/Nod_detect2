clear
path = uigetdir(pwd);
list = dir(fullfile(path, "*.csv"));
length_pre = 0;
for num = 1:length(list)
    name = list(num).name;
    filepath = fullfile(path, name);
    data = readmatrix(filepath);
    data_length = length(data);
    if(num == 1)
        Alldata(1:(length_pre+data_length), :) = data(:, 2:end);
    else
        Alldata(length_pre+1:(length_pre+data_length), :) = data(:, 2:end);
    end
    length_pre = data_length;

end

fm_Y = Alldata(:, 1);
fm_X = Alldata(:, 2:end);

Y_index1 = find(fm_Y == 1);
len_Y1 = length(Y_index1);

Y_index0 = find(fm_Y == 0);
len_Y0 = length(Y_index0);

Y_index00 = randsample(Y_index0,len_Y1);

Y_index = cat(1, Y_index1, Y_index00);
Y_index = sort(Y_index);

Y_extracted = fm_Y(Y_index, :);
X_extracted = fm_X(Y_index, :);

dataset = cat(2, Y_extracted, X_extracted);


csvwrite(fullfile(path, "dataset.csv"), Alldata);