% ----- Chạy 1 file -----
EEG = emotiv_edf_to_set('VEP-EDF/Apple/A1/sub3_A1.edf');

% ----- Hoặc chạy hàng loạt trong một thư mục -----
%{
files = dir('VEP-EDF/Apple/A1/*.edf');
for i = 1:numel(files)
    try
        emotiv_edf_to_set(fullfile(files(i).folder, files(i).name));
    catch ME
        fprintf(2, 'Error on %s: %s\n', files(i).name, ME.message);
    end
end
%}
