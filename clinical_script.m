% Clinical Batch Script

close all; clear; clc;
id = [12 14 15 16 17 18 19 21 22 23];

for i = 1 : length(id)
    id(i)
    imgname = sprintf('IRB_%03d',id(i));
    Exp_Clinical;
end