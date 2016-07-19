base_estimates_directory = '/home/marius/Documents/Database/DSD100/outputs/measures/'; 
songs_list = dir(fullfile(base_estimates_directory,'*.mat') ); 

methods_folders = dir(base_estimates_directory);
isub = [methods_folders(:).isdir]; %# returns logical vector
methods_folders = {methods_folders(isub).name}';
methods_folders(ismember(methods_folders,{'.','..'})) = [];

result = struct;
subsets_names = {'Test','Dev'};
for imethod_name = 1:length(methods_folders )
    for id=1:100
        maxDir=50; %50 songs for Dev and Test
        imethod_name = floor((id-1)/(maxDir*2))+1;
        rest = mod((id-1),maxDir*2)+1;
        subset_index = floor((rest-1)/maxDir)+1;
        song_index = mod(rest-1,maxDir)+1;
        
        method_name = methods_folders{imethod_name};
        estimates_folder = fullfile(base_estimates_directory,method_name,subsets_names{subset_index}); 

        song_name = songs_list(song_index+2).name;  
        load(fullfile(estimates_folder,song_name),'results');
        result.(lower(subsets_names{subset_index}))(song_index).results =results;
    end 
    %save the results for the method 
    result_file = fullfile(base_estimates_directory,[method_name,'.mat']);
    save(result_file,'result');
end
