function FPRAAT = loadpraat(praatfln)
% Input: path to text file containing formant
% estimation values yielded by PRAAT
% Output: matrix having a structure similar to that created by Michel
% Scheffers' FOREST algorithm for FORmant ESTimation:
% one line per frame, containing 11 columns. The first four formants 
% are in columns 4, 6, 8, and 10 respectively: 
%  TIME  (left empty here) RMS    GAIN(3rd column: left empty here) F1   B1     F2   B2     F3   B3     F4   B4
% The 12th and 13th column correspond to the 5th formant.

fid = fopen(praatfln, 'rt');
% excluding first 5 lines of header; recovering value in 6th line,
% which indicates how many frames there are
for k = 1:6
   tline = fgetl(fid);
end

nbframes = str2double(tline);
% excluding the next three lines (metadata)
for k = 1:3
    tline = fgetl(fid);
end

% loop for each frame
for i = 1:nbframes
    % first line: RMS amplitude. Needs to be converted from the
    % PRAAT notation: e.g. 3.6290753875823685e-06.
    tline = fgetl(fid);
    FPRAAT(i,2) = str2double(tline);

    % second line: number of formants detected
    tline = fgetl(fid);
    nfor = str2double(tline);

    for k = 1:nfor
        % For each formant, 2 lines: its frequency; its bandwidth.
            % Formant frequency: 
            FPRAAT(i,2 + 2 * k) = str2double(fgetl(fid));
            % Formant bandwidth:
            FPRAAT(i,3 + 2 * k) = str2double(fgetl(fid));
    end
end

% closing the file
fclose(fid)