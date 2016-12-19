% MEX all *.c files in the current folder

for file = dir('*.c')'; eval(['mex ' file.name]) ; end