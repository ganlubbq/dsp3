clear
fprintf('============================================================\n');
currentMachine  = computer('arch');
currentOS       = getenv('OS');
currentTime     = datestr(now);
currentPath     = pwd;

fprintf('Running on a %s machinen with %s! Current path is %s\\ \n',currentMachine,currentOS,currentPath);
fprintf('Setting paths...\n');

addpath([currentPath '/']);
addpath([currentPath '/c']);
addpath([currentPath '/lib']);
addpath([currentPath '/test']);

fprintf('Checking file changes...\n'); tic;

listing = dir([currentPath '/']);
for ii = 3:length(listing)
    if listing(ii).isdir
        listing(ii).name = [pwd '/' listing(ii).name '/'];
    else
        listing(ii).name = [pwd '/' listing(ii).name];
    end
end

fileNum = 0;
folderNum = 0;
while length(listing)>2
    if listing(3).isdir
        sublisting = dir(listing(3).name);
        for ii = 3:length(sublisting)
            if sublisting(ii).isdir
                sublisting(ii).name = [listing(3).name sublisting(ii).name '/'];
            else
                sublisting(ii).name = [listing(3).name sublisting(ii).name];
            end
            listing = [listing; sublisting(ii)];
        end
        listing = [listing(1:2); listing(4:end)];
        folderNum = folderNum + 1;
    else
        fprintf('%s\t%s\n',datestr(listing(3).datenum),listing(3).name);
        listing = [listing(1:2); listing(4:end)];
        fileNum = fileNum + 1;
    end
end; et = toc;

pause(0.1);

fprintf('Totally %d files and %d folders are checked in %.2f seconds.\n',fileNum,folderNum,et);
fprintf('Time now is %s.\n',currentTime);
fprintf('============================================================\n');


