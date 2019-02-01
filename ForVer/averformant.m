function [dataresamp] = averformant (data,samplenumber)
% A version of the function AVER adapted to operate on formant estimations.

processed_items = 0; 

% loop processing the syllables of the set.

for n = 1:length(data(1,1,:))
    n
    % re-initializing variables
    TR = []; datafile = [];
    
    % placing the n-th page of data into a matrix
	datafile = data(:,:,n);
	
    % condition to exclude empty pages (corresponding to syllables with
	% initial consonants, for which formant measurements were not made
	% separately)
    if ~isempty(nonzeros(datafile))
        % calculating the number of points at which measurements were made.
        % Because all pages of the 3-dimensional array have the same number of lines, 
        % there are some empty lines in most data files. So the file must be trimmed.
        ZT = datafile(:,5);
        % Note: the column investigated is 5 and not 1, as column 1 may
        % begin with zero on the first line: time zero, although there is
        % actually a measurement for this portion of signal.
        NZS = size(nonzeros(ZT)); % in previous version: used "length" and not "size"
        meas_nb = NZS (1);
        % storing this value in vector for later calculation
        PER (n) = meas_nb;

        for q = 1:11
          for r = 1:meas_nb
              TR(r,q) = datafile(r,q);
          end
        end
        datafile = [];
        datafile = TR; % re-assigning the trimmed data into (re-initialized) datafile

    % If there is a single data point in file: this value is assigned to the entire
    % resampled curve, but the user is made specifically aware of
    % this issue so a principled decision can be made, in view of the
    % linguistic data at issue and implications that this has for the
    % investigation. 
        if length(nonzeros(datafile(:,5))) < 2
            disp(['Only one single measurement within results file for item ',num2str(n),'.'])
            disp(['This value is assigned to all of the ',num2str(samplenumber),' slots.'])
            pause
            for np = 1:samplenumber
                for COLU = 1:length(datafile(1,:))
                    dataresamp(np,COLU) = datafile(1,COLU);
                end
            end
        else
            % calculating the length between detected closures at onset of first and
            % last periods
            LEN(n) = datafile(meas_nb,1) - datafile(1,1);
	
	
                                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                         %%%      interpolation       %%%
                                         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
            % calculating vector of time points at which samples at to be made for the
            % item
                   step = LEN(n)  /  (samplenumber - 1);
            for p = 1:(samplenumber-1)
                  timepts(p) = (datafile(1,1)+(p-1)*step);
            end
            timepts(samplenumber) = datafile(meas_nb,1);
            dataresamp (:,1,n) = timepts;
	
	
            for zz = 2:11
                % For formant values, there are empty slots = values that are set at
                % zero by convention, for want of a detectable formant value.
                if ~isempty (nonzeros(datafile(:,zz)) )
                    % finding out how many continuous sections there are
                    SE = 1;
                    % initializing variables used to count closings and openings of
                    % syllable portions for which there are continuous 
                    % measurements. Default values: beginning at first line; closing: left
                    % empty, treated later.
                    nbOPE = 0; 
                    if datafile(1,zz) ~= 0 
                        OPE(1) = 1;
                        nbOPE = nbOPE + 1;
                    end
                    nbCLO = 0; CLO = []; 
                    for p = 2 : meas_nb
                              if and((datafile(p - 1,zz) ~= 0),(datafile(p,zz) == 0))
                                  nbCLO = nbCLO + 1;
                                  CLO(nbCLO) = p - 1;
                              elseif and((datafile(p - 1,zz) == 0),(datafile(p,zz) ~= 0))
                                  nbOPE = nbOPE + 1;
                                  OPE(nbOPE) = p;
                              end
                    end
	
                    % if the last value is present, the last period is also a 
                    % "closing" point.
                    if datafile(meas_nb,zz) ~= 0
                        CLO(nbCLO + 1) = meas_nb;
                        nbCLO = nbCLO + 1;
                    end
	
                    % loop for continuous portions
                    for r = 1:nbOPE
                        % algorithm: calculate duration; find closest timepoints for beginning
                        % and end; interpolate; integrate into matrix.
	
                        % Retrieving beginning and end:
                        BEG = datafile(OPE(r),1);
                        EN = datafile(CLO(r),1);
	
                        % Find closest timepoints:
                        [valBEG,indexBEG] = min(abs(timepts - BEG));
                        [valEN,indexEN] = min(abs(timepts - EN));
	
                         % Case in which several values in succession are present: 
                         if CLO(r) - OPE(r) > 0
                            % interpolation 
                            re = interp1(datafile(:,1),datafile(:,zz),timepts(indexBEG:indexEN));
                            % placing the results in corresponding portion of matrix
                            for t = indexBEG : indexEN
                                dataresamp(t,zz,n) = re(t - (indexBEG - 1));
                            end
                         % Case of isolated value 
                         elseif CLO(r) == OPE(r)
                             % calculating the Number of Values (NV) that must be assigned
                             % in the reinterpolated curve: e.g. if there are 33 periods in
                             % the original measures, each period roughly corresponds to 3
                             % points in the reinterpolated curve; an isolated value in
                             % the original measures can therefore be assigned to 3
                             % successive points in the reinterpolated curve
                             NV = round (100 / meas_nb);
                             % calculation of PE: the point in 100-pt curve that corresponds most
                             % closely to the time position of the measured value. PE is
                             % a value in %.
                             PE = (datafile(OPE(r)) - datafile(1,1)) / LEN(n);
                             % transformation into an index out of 100:
                             INDEX = round (100 * PE);
                             % calculation of index where to begin the assignment of values:
                             LC = INDEX - round(NV / 2);
                             % assignment: 
                             for u = LC : (LC + NV - 1)
                                 % adding a condition in order not to go beyond the range of
                                 % existing indices, i.e. [1:100]
                                 if ismember(u,[1:100])
                                     dataresamp(u,zz,n) = datafile(OPE(r),zz);
                                 end
                             end
                         end
	
                    % end of loop for portions    
                    end
                % end of condition on non-emptiness of column
                end
            % end of loop for columns
            end
        % end of condition on non-uniqueness of measurement
        end
    % end of condition on non-emptiness of data sheet
    end
% end of loop for syllables
processed_items = processed_items + 1;
end