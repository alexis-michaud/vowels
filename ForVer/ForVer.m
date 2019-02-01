% This is a script for applying formant detection software 
% to a set of items and checking the results visually, one by one. 
% In its present form, the script takes as input the output of <peakdet>,
% a script for analysis of electroglottographic signals
% available from the following repository:
% https://github.com/alexis-michaud/egg 
%
% Under this workflow, the time codes (beginning
% and end) of intervals to be analysed are retrieved from the results file
% produced by <peakdet>. The <ForVer> script can easily be modified 
% to load the beginning and endpoints from another source, such as a text file.
%
% The script calls the standard BURG formant detection procedure
% implemented by PRAAT. 

% The script can also call the algorithm FOREST (software for FORmant ESTimation 
% developed by M. Scheffers et al., Kiel University). But preliminary
% results (from 2005) showed that FOREST raised some problems for /i/ and /u/. 
% So the current version leaves aside this method and just calls Praat.

% The values yielded by this analysis are plotted onto a broad-band
% spectrogram, so that the user can check the results visually.  
% (The raw results are also stored in the output matrices.)
%
% The script could be adapted to base the filtering on
% canonical values; on this method and its limits, see e.g. Gendrot &
% Decker 2004, who tested excluding values that differ from canonical values 
% by more than 200 Hz (same threshold for all formants). 
% Reference: 
% "Analyses formantiques automatiques de voyelles
% orales : évidence de la réduction vocalique en langues française et
% allemande", MIDL: Identification des langues et des variétés dialectales
% par les humains et par les machines (Paris, 2004), pp. 7-12.   
%
% Note that the parameters passed on to the formant-detection algorithms
% are suitable for male speakers and not for female speakers: the maximum
% formant value in PRAAT needs to be set to 5,500 Hz for female speakers.
% On the whole, FOREST seems to do less well for female voices.
%
% This script requires 
% - the Windows executable forest.exe created by Michel Scheffers
% - the software package PRAAT version 4.2 or above. 
% - the script <loadpraat.m>, to load the text file containing formant
% estimation values yielded by PRAAT
% - the script <dr.m> (for DRawing a figure), which itself requires
% - the scripts <myspecgram.m> (spectrogram calculation) and <plotspecgram>
% (self-explanatory: drawing the output of myspecgram.m)
% - the script <averformant.m> for making averages
% The user can check one to five formants. (Usually: I check three.)
% The script can easily be modified for correction of four
% formants, or of only one or two, as desired. This value,
% <nbcheck>, is set within the script (about line 51).
%
% The choices of settings & data are made within the script (about line 55).

% Clearing workspace
clear

% setting number of formants to be visually checked
nbcheck = 3;
disp(['Number of formants to be checked: ' num2str(nbcheck)])

% Important: PRAAT must be open when this script is run. If it's started
% automatically by the script, thus:
%     dos('C:\Program Files\MATLAB704\work\forver\praat.exe')
%     disp('PRAAT has been started.')
% then MATLAB does not accept any further commands. So you need to start it
% manually.
attent = input('Check that the PRAAT programme is running, then press ENTER.');

% Clearing figure
figure(1)
clf

% setting the highest frequency of spectrogramme at 6,000 as a default; this
% can be changed by the user as the programme is being run.
maxfreq = 6000;

% loading file that contains beginning and endpoint of relevant
% portions of the signal, and the number of the item: the results file
% created by running PEAKDET.
[textfilenew,textpathfilenew] = uigetfile('*.*','Please choose the .mat file that contains the results of EGG analysis');

% Reading the results file. 
load([textpathfilenew textfilenew])

if exist([textpathfilenew textfilenew(1:(length(textfilenew) - 4)) '.wav'])
    pathEGGnew = textpathfilenew;
    EGGfilenamenew = [textfilenew(1:(length(textfilenew) - 4)) '.wav'];
else
    [EGGfilenamenew,pathEGGnew] = uigetfile([textpathfilenew '*.*'], 'Please choose the recording to be downloaded');
end

% finding out the characteristics of the sound file
[Y,FS,NBITS] = wavread([pathEGGnew EGGfilenamenew],1);

siz = wavread([pathEGGnew EGGfilenamenew],'size');

% If there is a single channel: no difficulty.
% If the signal contains two or more channels: creating a new file containing only 
% the first channel (left channel of stereo file) by convention. If this
% file has already been created (e.g. by an earlier execution of this
% script), it is recognised.
if siz(2) > 1
    if exist([pathEGGnew 'MONO' EGGfilenamenew])
        EGGfilenamenew = ['MONO' EGGfilenamenew];    
    else
        % writing the audio into a mono file, for PRAAT to use (there is no
        % 'Open stereo long sound file' command in the version of the PRAAT
        % software used in this study, hence the need to demultiplex)
        [SIG,FS,NBITS] = wavread([pathEGGnew EGGfilenamenew]);
        SIG = SIG(:,1);
        wavwrite(SIG,FS,NBITS,[pathEGGnew,'MONO',EGGfilenamenew]);
        EGGfilenamenew = ['MONO' EGGfilenamenew];    
        % clearing SIG, a huge amount of data
        clear SIG; 
    end
end

% loading the sound file into PRAAT: individual vowels are then cut out of this file 
dos (['sendpraat praat "Open long sound file... ' pathEGGnew EGGfilenamenew '"']);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% loop for individual vowels
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iii = 1:length(data(1,1,:))
    % opening a figure for plotting data 
    figure(1)
    clf

    % clearing previous results
    clear F; clear Fcor; % clear FFOREST; clear FFORESTcor;
    
    % beginning and endpoint: first and last detected glottis-closure-instants,
    % i.e. beginning of interval (= LENG(iii,1) ) plus duration up to first GCI
    % (= data(1,1,iii) ).

    % Depending on the software that has been used to create the <data>
    % matrix, that matrix contains the times of beginning and end
    % either in absolute terms (relative to the total length of the
    % recording), IN MILLISECONDS, or in relative terms, IN SECONDS, 
    % in which case the values in
    % <LENG> must be added to obtain the absolute position.
    % To test for this: the time value of the last sheet in <data>
    % is evaluated. 
    MATI = data(1,1,length(data(1,1,:)));
    % number of items: 
    itnb = length(data(1,1,:));
    minT = data(1,1,iii);
    maxT = data(length(nonzeros(data(:,1,iii))),2,iii);
    if MATI/itnb > 500
        % The lower threshold is empirically set at 500 ms per item
        % (including the interval in-between two regions to be analysed);
        % this threshold is OK for the Naxi data with which I dealt, 
        % but may not work in every case.
        debtime = minT / 1000;
        endtime = maxT / 1000;
    else
        % If the interval is relative to the selected portion of
        % signal: adding the time value of the beginning of interval.
        % This is in the column before last of LENG.
        LENGli = length(LENG(1,:));
        if LENGli > 2
            LENGpt = LENGli - 1;
        elseif LENGli == 2
            LENGpt = 1;
        elseif LENGli == 1
            LENGpt = 1;
        end
        debtime = LENG(iii,LENGpt) / 1000 + minT;
        endtime = LENG(iii,LENGpt) / 1000 + maxT;
    end

    disp(['Currently treating item on line ' num2str(iii) '.'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% Using PRAAT to calculate formants
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Selecting sound file in PRAAT
    dos(['sendpraat praat "select LongSound ' EGGfilenamenew(1:length(EGGfilenamenew) - 4) '"']);

    % Putting together the command to PRAAT to make a sound excerpt. Time of beginning
    % and end: in seconds.
    commandepr = ['sendpraat praat "Extract part... ' num2str(debtime) ' ' num2str(endtime) ' no"'];
    % Running the 1st command
    dos(commandepr);
    % Calling the formant detection algorithm
    dos('sendpraat praat "To Formant (burg)... 0 5 5000 0.025 50"');


    % Sending the results to a text file
    if EGGfilenamenew(1:4) == 'MONO'
        praatfln = [pathEGGnew EGGfilenamenew(5:length(EGGfilenamenew) - 4) ...
             num2str(iii) 'PRAAT.Formant'];
    else
        praatfln = [pathEGGnew EGGfilenamenew(1:length(EGGfilenamenew) - 4) ...
             num2str(iii) 'PRAAT.Formant'];
    end
    dos(['sendpraat praat "Write to short text file... ' praatfln '"']);


    % Removing the Formants object
    dos(['sendpraat praat "select Formant ' EGGfilenamenew(1:length(EGGfilenamenew) - 4) '"']);

    dos(['sendpraat praat "Remove"']);


    % Removing the sound sample
    dos(['sendpraat praat "select Sound ' EGGfilenamenew(1:length(EGGfilenamenew) - 4) '"']);

    dos(['sendpraat praat "Remove"']);


    % loading results file
    F = loadpraat(praatfln);

    % adding time values to the formants matrix: the times have not
    % been written into the formant file created by PRAAT
    [PL,PC] = size(F);
    step = (endtime - debtime) / (PL - 1);
    for nn = 1:PL
        F(nn,1) = (nn - 1) * step; 
    end
    % making a copy of the results: a matrix into which
    % the filtered data will be placed 
    Fcor = F; 
    % This matrix has thirteen columns. The first four formants are in columns 4, 6,
    % 8, and 10 respectively. 

    % suppressing the text file created by PRAAT, to avoid clutter
    dos(['del ' praatfln]); 

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%% FORmant ESTimation 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % For 2019: use of the FOREST script is discontinued, as it has no
    % clear advantage in terms of accuracy. To turn this on again, some
    % 'housekeeping' is required in the code: bit-depth higher than 16-bit
    % requires preprocessing the audio file. Uncomment the lines below,
    % preprocess, and debug as necessary!
%
%     % Putting together the command to call FOREST. Time of beginning
%     % and end: in seconds.
%     commande = ['forest -b=',num2str(debtime),' -e=',num2str(endtime),' -oA ',...
%             pathEGGnew,EGGfilenamenew,' -of=',pathEGGnew,'form',num2str(iii),...
%             EGGfilenamenew(1:length(EGGfilenamenew) - 4),'.txt'];
%     % Running FOREST
%     dos(commande);


    %%%%%%%%%%%% Checking the results by looking at the spectrogram
    % recovering portion of signal
    [SIG,FS,NBITS] = wavread([pathEGGnew EGGfilenamenew],[round(debtime * FS) round(endtime * FS)]);

    % playing the sound
    sound(SIG,FS)


    %%% loading formant values. The procedure is somewhat roundabout: as the <load>
    %%% function cannot handle the header of the file, another file containing all
    %%% the lines except the first two (=except the header) is created and then
    %%% loaded by MatLab.

%     % indicating path to raw results
%     filename = [pathEGGnew,'form',num2str(iii),...
%             EGGfilenamenew(1:length(EGGfilenamenew) - 4),'.txt'];
%     % creating name of a "short" file, without header
%     fnameshort = [pathEGGnew,'formshort',num2str(iii),...
%             EGGfilenamenew(1:length(EGGfilenamenew) - 4),'.txt'];
%     % opening these files
%     fid = fopen(filename, 'rt');
%     fidshort = fopen(fnameshort, 'w');
%     % excluding first two lines (metadata/header)
%     for k = 0:1
%        tline = fgetl(fid);
%     end
%     % placing all the other lines, one by one, into the short file
%     while feof(fid) == 0
%        tline = '';
%        tline = fgetl(fid);
%        fprintf(fidshort,'%s',tline);
%        fprintf(fidshort,'\r');
%     end
%     % closing the files
%     fclose(fid);
%     fclose(fidshort);
% 
%     % loading the short file into a matrix
%     FFOREST = load(fnameshort); 
%     % RMS amplitude values are in the second column of the FOREST results.
% 
% 
%     %%%%%% deleting the text files, to avoid clutter
%     dos(['del ' filename]);      
%     dos(['del ' fnameshort]); 


    % returning an error message if the results are empty
    if isempty(F)
        error('No single formant value calculated for this item. Check that the interval of signal is of sufficient length.')
    end
    % Calculating the spectrogram of the sound
    [Timefreq, STFT, ret_n, n, offset, f_r, t_r] = myspecgram(SIG, FS);


    %%%%%%%%%%%%%%%%%%% Plotting the formants, and visual verification

    for ii = 1:nbcheck % At this stage: checking only the first 3 formants. 
        % The script can easily be modified to allow correction of four
        % formants, or of only one or two, as desired. This value,
        % <nbcheck>, is set at the top of the script.

        % clearing figure, and drawing
        figure(1)
        clf
        dr(Fcor, ii, maxfreq, FS, Timefreq, STFT, ret_n, n, offset, f_r, t_r)
        % All the formant values are shown in green stars; the values that
        % are being checked, or that have already been confirmed, are
        % in circles (yellow for F1, red for F2, green for F3).

        % creating a variable to serve for formant-confirmation loop 
        cornb = 6;
        while cornb ~= 0
            disp(['To validate the F',num2str(ii),' values displayed, type 0 (zero).'])
            disp(' ')
            disp('If some values need to be suppressed (i.e. if no formant value can be detected safely')
            disp('at some of the time points concerned), enter 1.')
            disp(' ')
            disp(['To indicate the approximate position of the formant F' num2str(ii) ', enter 2.'])
            disp('To replace all values by the values for the next formant (e.g. F3 for F2, ')
            disp('if the values indicated as "F3" actually correspond to F2), enter 3.')
            disp('To listen to the sound again, enter 4.')
            disp('To change the upper frequency threshold for the spectrogram')
            disp(['(which is presently at ',num2str(maxfreq),'), enter 5.'])
            disp('To exclude all values for this formant, enter 9.')
            disp(' ')

            % creating a condition to check that the 'CurrentPoint' is
            % chosen AFTER the above message has been displayed.
            azero = get(gca,'CurrentPoint');
            cornb = 6;
            while ~ismember(cornb,[0:5 9])
                cornb = input(' > ','s');
                % reducing input to one single character
                if length(cornb) > 0
                    cornb = cornb(1);
                else
                    cornb = 6;
                end
                if ismember(cornb,['0' '1' '2' '3' '4' '5' '9'])
                     cornb = str2num(cornb);
                else
                    cornb = 6;
                end
            end
            if cornb == 1
              % Need for a second loop, as this requires several
              % modifications on the part of the user
              while cornb == 1  
                disp('Now click on the data points on the figure that need to be suppressed (=set at 0).')
                disp('When you are finished, click to the left of the plot within the figure window, ')
                disp('i.e. at a point whose abscissa is less than that of the first data point.')
                v = axis;
                aaa(1,1) = v(1) + 1;
                % Tant que l'utilisateur ne clique pas à gauche du
                % cadre (pour signaler qu'il a fini) :
                while aaa(1,1) > v(1)
                    pause(0.1)
                    aaa = get(gca,'CurrentPoint');
                    % Tant que l'utilisateur n'a pas fourni une valeur
                    % différente de la valeur de départ : on attend
                    % qu'il clique
                    while aaa == azero
                        pause(0.1)
                        aaa = get(gca,'CurrentPoint');
                    end
                    % getting abscissa of point concerned
                    absci = aaa(1,1);
                    if aaa(1,1) > v(1)
                        % computing difference between this point and the time points at
                        % which formant values were computed
                        dis = Fcor(:,1) - Fcor(1,1) - absci;
                        % finding minimum of dis
                        absdis = abs(dis);
                        [Ydis,Idis] = min(absdis);
                        % assigning value zero to the corresponding point
                        Fcor(Idis,2 + 2 * ii) = 0;
                        % re-drawing
                        figure(1)
                        clf
                        dr(Fcor, ii, maxfreq, FS, Timefreq, STFT, ret_n, n, offset, f_r, t_r)
                    end
                end
                % on réinitialise la valeur <azero> pour l'item suivant
                azero = get(gca,'CurrentPoint');
                disp('  ')
                disp('  ')
                disp('  ')
                disp('  ')
                disp('  ')

                cornb = 6;
                while ~ismember(cornb,[0:5])
                    cornb = input('If it''s OK now, type 0. To suppress other values, type 1 > ','s');
                    % reducing input to one single character
                    cornb = cornb(1);
                    if ismember(cornb,['0' '1' '2' '3' '4' '5' '9'])
                         cornb = str2num(cornb);
                    else
                        cornb = 6;
                    end
                end
              end
            % If user wants to indicate approximative value manually:
            elseif cornb == 2
                  % Need for a second loop, as this requires several
                  % modifications on the part of the user
                  while cornb == 2  
                        disp(['Now click in the area where you detect visually formant ' num2str(ii) '.'])
                        disp('When you are finished, click to the left of the plot within the figure window, ')
                        disp('i.e. at a point whose abscissa is less than that of the first data point.')
                        v = axis;
                        aaa(1,1) = v(1) + 1;
                        % While the user does not click to the left of the
                        % figure (to indicate that verification is over): 
                        while aaa(1,1) > v(1)
                            pause(0.1)
                            aaa = get(gca,'CurrentPoint');
                            % While the value is not different from the
                            % initial value: waiting for the user to click
                            while aaa == azero
                                pause(0.1)
                                aaa = get(gca,'CurrentPoint');
                            end
                            % getting abscissa of point concerned
                            absci = aaa(1,1);
                            if aaa(1,1) > v(1)
                                % computing difference between this point and the time points at
                                % which formant values were computed
                                dis = Fcor(:,1) - Fcor(1,1) - absci;
                                % finding minimum of dis
                                absdis = abs(dis);
                                [Ydis,Idis] = min(absdis);
                                % looking for closest resonance for the corresponding point
                                Fapprox = aaa(1,2);
                                % Formants 1, 2, 3 and 4 are in columns 4, 6,
                                % 8, 10 of Fcor, respectively.
                                candidates = Fcor(Idis,[4 6 8 10]);
                                [YclosestF,IclosestF] = min(abs(candidates - Fapprox));
                                % if the value closest to that
                                % suggested by the user is in fact higher than 
                                % the original index (e.g. resonance labelled
                                % as "F3" must be substituted 
                                % for F2), then the F4 value must be
                                % substituted for F3, and F5 for F4; and F5
                                % must be set to zero.
                                if  IclosestF > ii
                                    % Loop: depends on number of
                                    % detected formants. In some cases,
                                    % PRAAT only yielded four. This is
                                    % why, in the loop, the ceiling for
                                    % <dec> is not 4 in all cases, but
                                    % is calculated as: 
                                    % round(length(Fcor(1,:))/2)
                                    for dec = (ii + 1):(round(length(Fcor(1,:))/2) - 3)
                                        Fcor(Idis,2 + 2 * dec) = Fcor(Idis,4 + 2 * dec);
                                    end
                                    % setting F5 to zero
                                    Fcor(Idis,2 + 2 * 5) = 0;
                                end
    %                             % If the value is lower than the original index
    %                             % (e.g. F2 value substituted for F3), then all
    %                             % the other formants must be pushed up one
    %                             % step. Starting from F4. But technically:
    %                             % would require looking for the 'lost formant'
    %                             % among the list of candidates proposed by
    %                             % PRAAT: one cannot select 3 formants from a
    %                             % list of 2. This is not developed in the
    %                             % present version of the script.
    %                             if  IclosestF < ii
    %                                 for dec = (ii):4
    %                                     Fcor(Idis,14 - 2 * dec) = Fcor(Idis,12 - 2 * dec);
    %                                 end
    %                                 % setting F5 to zero
    %                                 Fcor(Idis,2 + 2 * 5) = 0;
    %                             end

                                % putting the closest value into the corrected
                                % formant matrix
                                Fcor(Idis,2 + 2 * ii) = candidates(IclosestF);

                                % re-drawing
                                figure(1)
                                clf
                                dr(Fcor, ii, maxfreq, FS, Timefreq, STFT, ret_n, n, offset, f_r, t_r)
                            end
                        end
                        % on réinitialise la valeur <azero> pour l'item suivant
                        azero = get(gca,'CurrentPoint');
                        disp('  ')
                        disp('  ')
                        disp('  ')
                        disp('  ')
                        disp('  ')
                        cornb = 6;
                        while ~ismember(cornb,[0:5 9])
                            cornb = input('If it''s OK now, type 0. To correct other values, type 2 > ','s');
                            % reducing input to one single character
                            cornb = cornb(1);
                            if ismember(cornb,['0' '1' '2' '3' '4' '5' '9'])
                                cornb = str2num(cornb);
                            else
                                cornb = 6;
                            end
                        end
                  end
            elseif cornb == 3
                % number of column where to do replacement:
                COL = 2 + 2 * ii;
                for LI = 1:length(Fcor(:,1))
                        % formant frequency:
                        Fcor(LI,COL) = Fcor(LI,COL + 2);
                        % bandwidth:
                        Fcor(LI,COL + 1) = Fcor(LI,COL + 3);
                        % Also reporting the change onto the higher
                        % formants: all undergo downstepping
                        if COL + 5 < 12
                            % formant frequency:
                            Fcor(LI,COL + 2) = Fcor(LI,COL + 4);
                            % bandwidth:
                            Fcor(LI,COL + 3) = Fcor(LI,COL + 5);
                        end
                        if COL + 7 < 12
                            % formant frequency:
                            Fcor(LI,COL + 4) = Fcor(LI,COL + 6);
                            % bandwidth:
                            Fcor(LI,COL + 5) = Fcor(LI,COL + 7);
                        end
                end
                % re-plotting
                figure(1)
                clf
                dr(Fcor, ii, maxfreq, FS, Timefreq, STFT, ret_n, n, offset, f_r, t_r)
            elseif cornb == 4
                sound(SIG,FS);
            elseif cornb == 5
                % The value must be in-between reasonable limits
                maxfreq = 0;
                while or(maxfreq < 200,maxfreq > 20000)
                    maxfreq = input('Enter the new value, in Hz: > ');
                end
            elseif cornb == 9
                % number of column where to do replacement:
                COL = 2 + 2 * ii;
                for LI = 1:length(Fcor(:,1))
                        % formant frequency:
                        Fcor(LI,COL) = 0;
                end
                % re-plotting
                figure(1)
                clf
                dr(Fcor, ii, maxfreq, FS, Timefreq, STFT, ret_n, n, offset, f_r, t_r)
            end
        end   


%         % making corrections in the results given by FOREST, on the basis of
%         % the user-corrected PRAAT values. Acceptable differences: 150
%         % Hz for F1, 300 Hz for F2, 450 Hz for F3. If difference is
%         % greater: set all values at 0 in FOREST results.
%         FFORESTcor = FFOREST;
%         DIF = mean(nonzeros(FFOREST(:,2 * ii + 2)) - mean(nonzeros(Fcor(:,2 * ii + 2))));
%         if abs(DIF) > 150 * ii
%             FFORESTcor(:,2 * ii + 2) = 0;
%         end
    % end of loop for formants 1 to 3
    end

    %%%%%%%% placing data in matrix
    % checking that there is no doubling of the last line (this occasional
    % problem results in a bug that I have not identified, which causes
    % the last line to be written twice into the <datafile> matrix)
    ld = length(Fcor(:,1));
    if ld > 1
        if Fcor(ld,:) == Fcor(ld - 1,:)
            Fcor = Fcor(1:ld - 1,:);
        end
    end
    % calculating the number of measurement points (= nb of lines)
    period_nb = size(Fcor,1);
    % calculating the number of columns
    nbcol = length(Fcor(1,:));
    % assigning values in data matrix
    for q = 1:length(F(1,:))
        for r = 1:period_nb
            dataF(r,q,iii) = F(r,q);
            dataFcor(r,q,iii) = Fcor(r,q);
        end
    end

%     % same for values calculated by FOREST: the raw values (in
%     % dataFPRAAT), and the filtered values (in dataFPRAATcor).
%     % calculating the number of measurement points (= nb of lines)
%     period_nb = size(FFOREST,1);
%     % calculating the number of columns
%     nbcol = length(FFOREST(1,:));
%     % assigning values in data matrix
%     for q = 1:nbcol
%         for r = 1:period_nb
%             dataFFOREST(r,q,iii) = FFOREST(r,q);
%             dataFFORESTcor(r,q,iii) = FFORESTcor(r,q);
%         end
%     end
%% end of loop for individual vowels
end

% Calculating averages from the data
disp('Averaging raw PRAAT data.')
dataresamp = averformant (dataF,100);
disp('Averaging corrected PRAAT data.')
dataresampcor = averformant (dataFcor,100);
% disp('Averaging raw ForEst data.')
% dataresampF = averformant (dataFFOREST,100);
% disp('Averaging filtered ForEst data.')
% dataresampFcor = averformant (dataFFORESTcor,100);

% Sorting into individual matrices for each vowel needs to be done further
% down the line, on the basis of this script's output. 

% removing LongSound object from PRAAT
dos(['sendpraat praat "select LongSound ' EGGfilenamenew(1:length(EGGfilenamenew) - 4) '"']);
dos(['sendpraat praat "Remove"']);

% clearing temporary variables: signal portion, results sheet etc.
clear datafile; clear F; clear Fcor; % clear FFOREST; clear FFORESTcor; 

%% saving final results, adding a '_withF' suffix (for 'with formant estimation')
save([textpathfilenew textfilenew(1:(length(textfilenew) - 4)) '_withF.mat'])

% Farewell messages
disp('Total number of formant data points corrected or set at zero by user:')
disp(length(nonzeros(dataF - dataFcor)))
disp('i.e. proportion of data points changed by user (in percent):')
disp(100 * (length(nonzeros(dataF - dataFcor)))/(3 * length(nonzeros(dataF(:,4,:)))))
