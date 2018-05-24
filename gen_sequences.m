function [sqc] = gen_sequences(subj, NSEQ)
%% making sequences
% (1) define cond, length of seq, agency cond, target colors
% (2) make sequences for first condition
% (3) mirror them for second condition
% (4) shuffle with constraints on repeating conditions and symbol pairs:
fprintf('>>> running gen_sequences\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO:
%%% //what if nseq_condition < symbolpairs ?\\
%%% do we control for nb of time any symbol is a target?
%%% add constraint for not three times the same sequence condition?? and not several times the same symbol??
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %% arg check: 
% if nargin < 1
%     subj = 99;
%     NSEQ = 32;
%     warning('Cle: Will take default subj (99) and nseq = (32)');
% elseif nargin == 1
%     NSEQ = 32;
%     warning('///Will take default nseq = (32)');
% end
% 
% if mod(NSEQ, 64) ~= 0
%     error('Cle: cannot balance probeside with %d trials, choose multiple of 64 (64 128 192)\n', NSEQ);
% end

% %32 if not controlling for good response side (32, 64, 96, 128, 192)
% %64 if controlling for motor response (64, 128, 192)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% defining conditions:

sqc = struct();

%%% target/open condition = 1s or 2s
nseq_condition  = NSEQ/2;
condition       = [repelem(1,nseq_condition) repelem(2,nseq_condition)]; %size=1*64 (FULL)

%%% lengths of sequences = 8,12,16,20 samples per sequence (MUST BE EVEN)
nlenghts        = 4;
nseq_length     = nseq_condition/nlenghts; % nb of seqs per length condition (4different lengths)
nsamples        = [repelem(8,nseq_length) repelem(12,nseq_length) repelem(16,nseq_length) repelem(20,nseq_length)]; %size=HALF

%%% agency on color = for each length : half counterfactual(1), half always same color(2)
nagt_length   = nseq_length/2;
agencyoncolor = repmat([repelem(1, nagt_length) repelem(2, nagt_length)], [1,nlenghts]); %size=HALF

%%% color target = for each agency condition, half is blue(1)/half is red(2)
ncolors_agency = nagt_length/2;
targetcolor    = repmat([repelem(1, ncolors_agency) repelem(2, ncolors_agency)], [1,nlenghts*2]); %size=HALF

%%% probeside is not controlled yet...
if mod(NSEQ, 64) == 0
    nprobeside = ncolors_agency/2;
    probeside = repmat([repelem(1, nprobeside) repelem(2, nprobeside)], [1,nlenghts*2*2]);
elseif mod(NSEQ, 64) ~= 0
    probeside = randi(2,[1,nseq_condition]);
end

%%% get a pool of symbol pairs with least amount of repeats, shuffled
symbolpairs = nchoosek(1:8,2); % 28 unique possible symbol pairs
npairs      = size(symbolpairs,1); % get pool to make sure not too many repeats
if nseq_condition <= npairs %if we need fewer pairs, take random pairs, no need to control > we'll only use to pilot loosely
    pool = symbolpairs(randperm(npairs, nseq_condition),:);
elseif nseq_condition > npairs % if we need more pairs, try to have least nb of repeats
    padsize = randperm(npairs, mod(nseq_condition,npairs)); % choose unique pairs to add up at the end
    pool = [repmat(symbolpairs, [1, floor(nseq_condition/npairs)]); symbolpairs(padsize,:)]; % pool as many full sets of all pairs + padding
end
pool = pool(randperm(size(pool,1)),:);
%tabulate([pool(:,1);pool(:,2)])

%%% draw col values: no control of difficulty of seq in the different conditions
colorvalues     = fx_getSeqColorValues(nsamples); %function to draw 0-1 color values
cs = load('color_space.mat');
cs.vvec = linspace(0,1,size(cs.xvec,2));

%% build sequences info in first condition (target):
for iseq = 1:nseq_condition
    
    sqc(iseq).subj = subj;
    
    sqc(iseq).condition = condition(iseq); 
    sqc(iseq).mirroring = iseq;
    
    sqc(iseq).nsamples = nsamples(iseq);
    
    sqc(iseq).agtcondition = agencyoncolor(iseq);
    
    sqc(iseq).targetcolor = targetcolor(iseq);
    
    sqc(iseq).symbolpair = pool(iseq,:);
    targetsymbol = randi(2,1); % randomly choose who's target
    sqc(iseq).targetsymbol = pool(iseq,targetsymbol);
    sqc(iseq).othersymbol = pool(iseq,3-targetsymbol);
    
    % assign colors to each side depending on counterf and targetcol:
    if agencyoncolor(iseq) == 1
        sqc(iseq).colorslin(1,:) = colorvalues(iseq).colors';
        sqc(iseq).colorslin(2,:) = (-1) * colorvalues(iseq).colors';
    elseif agencyoncolor(iseq) == 2
        if targetcolor(iseq) == 1 % both are blue
            sqc(iseq).colorslin(1,:) = colorvalues(iseq).colors';
            sqc(iseq).colorslin(2,:) = sqc(iseq).colorslin(1,:);
        elseif targetcolor(iseq) == 2 % both red
            sqc(iseq).colorslin(1,:) = (-1) * colorvalues(iseq).colors';
            sqc(iseq).colorslin(2,:) = sqc(iseq).colorslin(1,:);
        end
    end
    
    %compute rgb triplets of drawn colors:
    rgbcol = nan((nsamples(iseq)),3,2); %rows = samples, cols = rgb, 3D = blue/red
    for col = 1:nsamples(iseq)        
        rgbcol(col,:,1) = fx_get_rgb(sqc(iseq).colorslin(1,col), cs);
        rgbcol(col,:,2) = fx_get_rgb(sqc(iseq).colorslin(2,col), cs);
    end
    sqc(iseq).rgb = rgbcol(:,:,:);
        
    % randomise symbol side, nrepeated side = 3
    nrepeated = 3; %
    symside = zeros(1,nsamples(iseq)); % init asmany samples:
    while abs(diff(hist(symside,1:2)))>1 % while not equal number of 1s and 2s
        symside(1:nrepeated) = randi(2,1,nrepeated); %fill the first three
        for i = (nrepeated+1):nsamples(iseq) % loop on following ones:
            ipos = randi(2,1); % get one random
            if sum(symside((i-nrepeated):(i-1))==ipos)>=nrepeated %if last 3 were already == thisside
                symside(i) = 3-ipos;%take other side
            else
                symside(i) = ipos;%else take this side
            end
        end
    end
    
    sqc(iseq).symbolsides = symside;
    sqc(iseq).probeside = probeside(iseq);
    sqc(iseq).bluebuttonside = randi(2,1);
    
end%end for seq in condition #1

%% duplicate for 2nd condition (open condition = exact same sequences)

for iseq = nseq_condition+1 : NSEQ
    iduplicat                = iseq - nseq_condition;

    sqc(iseq).subj           = sqc(iduplicat).subj;

    sqc(iseq).condition     = condition(iseq);
    sqc(iseq).mirroring      = iduplicat;
    
    sqc(iseq).nsamples       = sqc(iduplicat).nsamples;
    sqc(iseq).agtcondition   = sqc(iduplicat).agtcondition;
    sqc(iseq).targetcolor    = sqc(iduplicat).targetcolor;

    
    sqc(iseq).symbolpair     = sqc(iduplicat).symbolpair;
    sqc(iseq).targetsymbol   = sqc(iduplicat).targetsymbol;
    sqc(iseq).othersymbol    = sqc(iduplicat).othersymbol;
    
    sqc(iseq).colorslin      = sqc(iduplicat).colorslin;
    sqc(iseq).rgb            = sqc(iduplicat).rgb;
    
    sqc(iseq).symbolsides    = sqc(iduplicat).symbolsides;
    sqc(iseq).probeside      = sqc(iduplicat).probeside;
    sqc(iseq).bluebuttonside = sqc(iduplicat).bluebuttonside;

end%function

%% random permute: condition not three of the same sequence??

shuffle = randperm(NSEQ);
sqc = sqc(shuffle);

%%
% %randomise order, controlling for not twice same symbol in a row:
% solution = false;
% retryloopingoverpairs=0
% while solution == false % if doesn't find ordering this time, retry at pair 1
%     
%     retryloopingoverpairs = retryloopingoverpairs + 1
%     seq pairs = NaN(nseq_condition,2);
%     currentpool = pool;
%     idx = randi(length(pool),1);
%     seqpairs(1,:) = currentpool(idx,:); % draw a first pair to start with 
%     currentpool(idx,:) = []; % and remove it from remaining pairs
%     
%     for tr = 1:nseq_condition
%         currentpair = seqpairs(tr,:); %current pair
%         templist = currentpool(sum(ismember(currentpool, currentpair),2) == 0,:); %templist of all pairs not containing the two symbols drawn
%         if size(templist,1) == 0; break; end% if no pair is different, break for loop, start over at begining
%         seqpairs(tr+1,:) = templist(randi(size(templist,1),1),:); %draw next pair from this temp list
%         removefrompool = find(sum(ismember(currentpool,seqpairs(tr+1,:)),2) == 2,1); %in case there are several occurrences of the same pairs: remove only first
%         currentpool(removefrompool,:) = []; %remove the drawn new pair from the list of remaining pairs to place
%     end% for each seq, get different pair
%     
%     if tr == nseq_condition && sum(ismember(seqpairs(1,:),seqpairs(nseq_condition,:)),2) == 0; solution = true; end
% end % while no solution, retry looping over trials
% 
end%function

function [rgb] = fx_get_rgb(v, cs)
    % get R/G/B guns for desired position on color space
    v = (v+1)/2;
    v = min(max(v,0),1);
    rgb = nan(1,3);
    for i = 1:3
        rgb(i) = interp1(cs.vvec,cs.xvec(i,:),v,'pchip');
    end
end