function fixedTimes = checkDIOTimes(nTrials, dioTimes, commentTimes, maxCommentDiff)


if nargin <= 3
    maxCommentDiff = 0.5;
end

nDIOTimes = length(dioTimes);
nComments = length(commentTimes);

if nDIOTimes ~= nTrials
    fprintf(1, 'WARNING: Number of trial onsets does not match the number of trials, trying to fix it.\n');
    % There are more trial onset times than there are trials!
    % First, check to see if comments is correct:
    if ~isempty(commentTimes)
        if nComments ~= nTrials
            error(['Something is wrong with this dataset, the number of trial ' ...
                'DIO times does not match the number of trials or the number of comments! ' ...
                'nDIO=%0.0f nComments=%0.0f nTrials=%0.0f\n'], length(dioTimes), length(commentTimes), nTrials);
        end
        
        if isempty(dioTimes)
            fprintf(2, ['WARNING: There are no DI events!']);
            fixedTimes = commentTimes;
            return;
        end
        
        % Old way, check all dio times to see if they are close enough to a
        % comment.
        %         for iComment = 1: length(commentTimes)
        %             commentTime = commentTimes(iComment - nRemoved);
        %             onsetTime   = dioTimes(iComment - nRemoved);
        %
        %             if abs(onsetTime - commentTime) > maxCommentDiff
        %                 % This onset is too far from a comment.
        %                 dioTimes(iComment - nRemoved) = [];
        %                 nRemoved = nRemoved + 1;
        %             end
        %         end
        
        % New way, for each comment find the closest dio event.
        
        fixedTimes = zeros(1, nTrials);
        for iTrial = 1: nTrials
            
            commentTime = commentTimes(iTrial);
            diffTimes   = dioTimes - commentTime;
            [closestDIODiff, closestDIOIndex] = min(abs(diffTimes));
            
            if closestDIODiff > maxCommentDiff
                % Some problem here, just go with comment times.
                %fprintf(2, 'WARNING: Closest comment to DIO event was more than %0.2f secs away! Using comment times.\n', maxCommentDiff);
                %fixedTimes = commentTimes;
                %return;
                fixedTimes(iTrial) = commentTime;
            else
                try
                fixedTimes(iTrial) = dioTimes(closestDIOIndex);
                dioTimes(closestDIOIndex) = [];
                catch
                    fixedTimes(iTrial) = commentTime;
                end
            end
            
            
        end
        fprintf(1, 'Fixed hopefully.\n');
    else % hack for two array set up with no comments.
        % find the most common time between dioTimes and use that...

        fixedTimes = zeros(1, nTrials);
        old_dio = dioTimes(1);
        fixedTimes(1) = old_dio; % assume first is right...
        for iTrial = 2: nTrials
            diffTimes   = abs(dioTimes - (old_dio+maxCommentDiff)); % should be close to next stim??
            [closestDIODiff, closestDIOIndex] = min(abs(diffTimes));
            if dioTimes(closestDIOIndex) == old_dio 
                fixedTimes(iTrial) = old_dio + maxCommentDiff;
            else
            fixedTimes(iTrial) = dioTimes(closestDIOIndex);
            end
            old_dio = fixedTimes(iTrial);
        end
    end
else
    fixedTimes = dioTimes;
end

end