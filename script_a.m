%%
clc;clear all
addpath('~/expfys2/termo/matningar/Data/');

files_t=dir('~/expfys2/termo/matningar/Data/Al*t.lvm');
files_v=dir('~/expfys2/termo/matningar/Data/Al*v.lvm');


for i = 1:length(files_v)
    files_t(i).name
    files_v(i).name
    time=load(files_t(i).name, '-ascii');
    volt=load(files_v(i).name, '-ascii');
    
    
end
