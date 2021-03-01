%% root_locus examples
%%
sys = tf(1,[1 3 2 0]);
root_locus(sys);
%%
sys = tf([1 4] ,[1 2 0]);
root_locus(sys);
%%
sys = tf([1 0] ,[1 5 4 20]);
root_locus(sys);
%%
sys = tf([1 0 1] ,[1 1 0]);
root_locus(sys);