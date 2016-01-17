%%
clc;clf;clear all

f=@(x) sin(2*x)-.5*cos(x);
t_start=0;
t_end=4*pi;

test_movie(f, t_start, t_end )

