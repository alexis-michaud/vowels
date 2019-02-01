function dr(F, ii, maxfreq, FS, Timefreq, STFT, ret_n, n, offset, f_r, t_r)
% Drawing results for <forver.m>.
% F: matrix of formant estimations.
% ii: index of formant being plotted.
% Other variables: data for spectrogram plotting.

	% 	%% The 4 lines of code below correspond to code used at an earlier stage; the
	% 	%% resulting spectrogram was not the broad-band spectrogram necessary for
	% 	%% speech analysis.
	% 	spectrogram(SIG,FS,145)
	% 	ax = axis;
	% 	% setting top line of spectrogram at 8,000 Hz
	% 	ax(4) = 8;
	% 	axis(ax);

% plotting the spectrogram
plotspecgram (maxfreq, n, FS, Timefreq, STFT, ret_n, offset, f_r, t_r)
hold on

for i = 1:ii
    %% plotting the values. They do not need to be divided by 1,000, as the
    %% scale is in Hz and not in KHz.
    if i == 1
        coul = 'y';
    elseif i == 2
        coul = 'r';
    elseif i == 3
        coul = 'g';
    end
    plot(F(:,1) - F(1,1), F(:,2 + 2 * i),'o','MarkerEdgeColor','k','MarkerFaceColor',coul,...
                'MarkerSize',10)
    % plotting extra points for first and last data points, which may have
    % been excluded by limits of the axes
	ax = axis;
	if F(1,1) < ax(1)
        plot(ax(1),F(1,2 + 2 * i),'o','MarkerEdgeColor','k','MarkerFaceColor',coul,...
                'MarkerSize',10);
	end
	if F(length(F(:,1)),1) > ax(2)
        plot(ax(2),F(length(F(:,1)),2 + 2 * i),'o','MarkerEdgeColor','k','MarkerFaceColor',coul,...
                'MarkerSize',10);
	end
    % special treatment for single-data points: plotting at centre
    if length(nonzeros(F(:,4))) == 1
        plot(((ax(2) - ax(1)) / 2),F(1,2 + 2 * i),'o','MarkerEdgeColor','k','MarkerFaceColor',coul,...
                'MarkerSize',10)
    end
end

% Plotting higher (unverified) formants in black
if ii < 4
    coul = 'b';
    for j = (ii+1):4
         plot(F(:,1) - F(1,1), F(:,2 + 2 * j),'o','MarkerEdgeColor','k','MarkerFaceColor',coul,...
                'MarkerSize',7)
    end
end