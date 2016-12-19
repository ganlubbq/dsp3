% MEX all *.c files in the current folder

for file = dir('*.c')'; eval(['mex -v -g ' file.name]) ; end