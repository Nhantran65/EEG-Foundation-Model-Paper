function EEG = emotiv_edf_to_set(edfFile, outDir)
% Read Emotiv EDF, keep 14 EEG channels, build events from Marker*,
% and save EEGLAB .set with proper EEG.event.
%
% Example (run từ Command Window hoặc script khác):
%   EEG = emotiv_edf_to_set('VEP-EDF/Apple/A1/sub3_A1.edf');

if nargin < 2 || isempty(outDir), outDir = fileparts(edfFile); end

% --- deps
eeglab nogui;
addpath(fullfile(fileparts(which('eeglab.m')),'plugins','biosig'));

% --- read whole file with biosig
[dat, hdr] = sload(edfFile);             % dat: N x Ch
labels = string(hdr.Label);

% --- whitelist 14 EEG channels (Emotiv EPOC layout)
EEG_CH = ["AF3","F7","F3","FC5","T7","P7","O1","O2","P8","T8","FC6","F4","F8","AF4"];

% Lấy đúng theo thứ tự EEG_CH (đảm bảo không sai thứ tự)
ixEEG = arrayfun(@(nm) find(labels==nm,1), EEG_CH);
if any(isnan(ixEEG) | ixEEG==0)
    missing = EEG_CH(isnan(ixEEG) | ixEEG==0);
    error('Không tìm thấy các kênh EEG sau trong EDF: %s', strjoin(missing, ', '));
end

ixIdx  = find(labels=="MarkerIndex",     1);
ixType = find(labels=="MarkerType",      1);
ixVal  = find(labels=="MarkerValueInt",  1);

assert(~isempty(ixIdx) && ~isempty(ixType) && ~isempty(ixVal), ...
    'Thiếu 1 trong 3 kênh MarkerIndex/MarkerType/MarkerValueInt.');

% --- build RAW EEG matrix (channels x samples) cho EEGLAB
X = dat(:, ixEEG).';                     % (14 x N)
fs = isscalar(hdr.SampleRate) * hdr.SampleRate + ...
     ~isscalar(hdr.SampleRate) * mode(hdr.SampleRate(:));

% --- import vào EEGLAB (truyền biến X, không phải chuỗi 'X')
EEG = pop_importdata('dataformat','array','data', X, ...
    'srate', fs, 'nbchan', numel(ixEEG), 'pnts', size(X,2), 'xmin', 0);

% Đặt tên gọn
[~, base] = fileparts(edfFile);
EEG.setname = base;

% --- gán label kênh theo đúng thứ tự
for k = 1:numel(EEG_CH)
    EEG.chanlocs(k).labels = EEG_CH(k);
end

% --- trích event từ 3 kênh marker
eventsTbl = extract_emotiv_markers_from_arrays( ...
                 dat(:,ixIdx), dat(:,ixType), dat(:,ixVal), fs);

% In số lượng sự kiện & 5 event đầu
fprintf('Detected %d events in file: %s\n', height(eventsTbl), edfFile);
if height(eventsTbl) > 0
    disp(eventsTbl(1:min(5,height(eventsTbl)), :));
end

% --- đẩy vào EEG.event (type để là value; có thể đổi theo nhu cầu)
EEG.event = struct('type',{},'latency',{},'urevent',{});
for i = 1:height(eventsTbl)
    EEG.event(i).type    = num2str(eventsTbl.value(i));  % '1','2',...
    EEG.event(i).latency = double(eventsTbl.pos(i));     % sample index (1-based)
    EEG.event(i).urevent = i;
end
EEG = eeg_checkset(EEG,'eventconsistency');

% --- save .set
outfile = fullfile(outDir, [base '.set']);
EEG = pop_saveset(EEG, 'filename', [base '.set'], 'filepath', outDir);
fprintf('Saved: %s  (events: %d)\n', outfile, numel(EEG.event));
end


% ===================== Helper ======================
function events = extract_emotiv_markers_from_arrays(mIdx, mTyp, mVal, fs)
% Input: 3 vector cột (N x 1) của MarkerIndex/Type/ValueInt và Fs
% Output: table {pos, t_sec, index, type, value}

% ép integer (Emotiv hay lưu float)
mIdx = int64(round(mIdx));
mTyp = int64(round(mTyp));
mVal = int64(round(mVal));

% 1) phát hiện lần thay đổi trạng thái (so với mẫu trước)
chg = [true; any(diff(double([mIdx mTyp mVal])) ~= 0, 2)];
pos = find(chg);                         % sample indices (1-based)

idx_at = mIdx(pos); 
typ_at = mTyp(pos); 
val_at = mVal(pos);

% 2) lọc: bỏ trạng thái rỗng 0-0-0
if ~isempty(pos)
    valid = ~(idx_at==0 & typ_at==0 & val_at==0);
    pos    = pos(valid);
    idx_at = idx_at(valid);
    typ_at = typ_at(valid);
    val_at = val_at(valid);
end

% Nếu sau lọc không còn gì -> trả bảng rỗng an toàn
if isempty(pos)
    events = table('Size',[0 5], ...
        'VariableTypes', {'double','double','int64','int64','int64'}, ...
        'VariableNames', {'pos','t_sec','index','type','value'});
    return;
end

% 3) chống spam: chỉ nhận "rising" khi MarkerIndex tăng
% (đảm bảo rising có cùng kích thước với pos)
if numel(idx_at) >= 2
    rising = [true; diff(idx_at) > 0];
else
    rising = true(size(idx_at));  % chỉ 1 event -> giữ lại
end

% Áp rising
pos    = pos(rising);
idx_at = idx_at(rising);
typ_at = typ_at(rising);
val_at = val_at(rising);

% Nếu tiếp tục rỗng -> trả bảng rỗng
if isempty(pos)
    events = table('Size',[0 5], ...
        'VariableTypes', {'double','double','int64','int64','int64'}, ...
        'VariableNames', {'pos','t_sec','index','type','value'});
    return;
end

% 4) tính thời gian giây
t_sec = (pos-1) ./ fs;

% 5) xuất bảng
events = table(pos, t_sec, idx_at, typ_at, val_at, ...
    'VariableNames', {'pos','t_sec','index','type','value'});
end
