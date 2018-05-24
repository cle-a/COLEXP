function [colorvalues] = fx_getSeqColorValues(nsampleslist)
% returns structure of the color draws for each sequence
% (1) get nsim sequence simulations
% (2) for each specific seq length (=nsamples in the sequence):
%       define how many seq of this length we want to keep
%       draw xmany bundles of this many sequences
%       see which bundle is closest to mean difficulty of all the draws
%       keep bundle with least distance to mean
% (3) return colovalues structure = 
% seqid     = id of each sequence (1:length(nsampleslist)
% nsamples  = nb of samples in the seq (8, 12, 16, 20)
% seqnb     = nb of this seq amongst the nseq of this length (1:4?)
% colors    = colorvalues (sz = nsamples)
% // NB: RETURNS 0-1 COLORS => TRANSFORM TO ISOLUMINANT RGB AFTERWARDS \\

fprintf('>>> running fx_getSeqColorValues\n')

%% parameters
tic;
fprintf('Drawing color values...\n')

%get final nb of color draws to get:
%nsampleslist = [repelem(8,8) repelem(12,8) repelem(16,8) repelem(20,8)]; %1*32
seqlengths   = unique(nsampleslist); %
nseqlength   = numel(seqlengths);

nsim         = 10e4;  % number of seq simulations
ndraws       = 10e5;  % number of draws of nseq sequences to get mean

% model parameters for nsim simulations:
modeltype = 'targseek'; % target-seeking or information-seeking sampling of forms
alpha     = 0.4; % leak / learning rate from 0 to 1
zeta      = 1; % counterfactual strength (from 0=no update of unchosen form to 1=fully counterfactual update of unchosen form)
sigma     = 0.3; % leave sigma to 0.3 for now

% get colors distributions (Valentins):
PERR  = 1/3;   % probability of misclassification of the form based on a single sample/color
x = -1:0.002:+1;
x = 0.5*(x(1:end-1)+x(2:end));
b = fzero(@(p)trapz(x(x > 0),1./(1+exp(x(x > 0)*p))/trapz(x,1./(1+exp(x*p))))-PERR,[0,1000]);
p1 = 1./(1+exp(-x*b)); p1 = p1/trapz(x,p1);
%p2 = 1./(1+exp(+x*b)); p2 = p2/trapz(x,p2);

% init
colorvalues = struct(); % color values of drawn blocks
stored = 0;

%% for each seqlength draw 10e5 then draw sets of xmany sequences:
for iseqlenght = 1:nseqlength %for each possible seq length (=nsample)
    fprintf('draw colors for seqs of length #  %d/%d...\n', iseqlenght, nseqlength)

    %1// get the matrix of 10e5 simulations of 1 bloc
    nsmp = seqlengths(iseqlenght);%length = nb of sampling decisions in seq
    nseq = numel(find(nsampleslist == nsmp));%how many seq/blocks of this length
    
    vs = zeros(nsmp+1,nsim,2); % internal color values associated with each form at each moment in time
    rs = zeros(nsmp,nsim); % sampling decision (which form to look at)
    cv = zeros(nsmp,nsim); % color values presented at each sample
    
    for ismp = 1:nsmp% for each sample in sequence:
        %%% <a> choose which symbol to sample from
        if ismp == 1 % if first one, choose randomly
            rs(ismp,:) = 1 + (rand(1,nsim) < 0.5);
        else
            switch modeltype %else depends on models' goals:
                case 'targseek'% symbol with highest value = chooses the positive color
                    [~,rs(ismp,:)] = max(vs(ismp,:,:),[],3);%dim3 = 2symbols
                case 'infoseek'% symbol with least information
                    [~,rs(ismp,:)] = min(abs(vs(ismp,:,:)),[],3);
            end%SWITCH
        end%IF: decide how to sample/choose which form to look at
        
        %%% <b> update color values associated with each
        s = itrnd([1,nsim],x,p1,'exact');% sampling from p1: gets a color
        cv(ismp,:) = s;
        s = cat(3,s,-s);% adds value of p2 (is the opposite of p1 = symetry)
        vold = vs(ismp,:,:);% values from the last draw/stimulus
        % computing the new values (vnew):
        vnew = nan(1,nsim,2);
        % get linear indices of color values for chosen/unchosen symbols:
        ich = sub2ind(size(vold),ones(1,nsim),1:nsim,rs(ismp,:));
        iun = sub2ind(size(vold),ones(1,nsim),1:nsim,3-rs(ismp,:));
        % update expected color values assigned to each symbol:
        vnew(ich) = vold(ich) +    (alpha  +   sigma*randn(1,nsim)).*(+s(ich)-vold(ich));
        vnew(iun) = vold(iun) +    (alpha  +   sigma*randn(1,nsim)).*(-s(ich)*zeta-vold(iun));
        % update leanrt colors:
        vs(ismp+1,:,:) = vnew; % update *internal/expected* color values for next trial
    end
    
    %2// draw groups of x many & store color values
    iseq        = randi(nsim,[nseq ndraws]); %table of nseq x ndraws sequences idx
    idx         = 1:numel(iseq); % get indices of seq from the 10e5
    dseq        = reshape(vs(end, iseq(idx), 1) - vs(end,  iseq(idx), 2) , [nseq ndraws]); %table of nseq x ndraws sequences diff
    cl          = cv(:,iseq(idx));
    cl = reshape(cl, [nsmp, nseq, ndraws]);
    
    %3// sort, retreive blocks with least distance to mean => TODO store colvals
    sseq            = sort(dseq, 1); %sort diff by difficulty
    avg             = mean(sseq, 2); %average difficuly for each diff rank = rows
    distseq         = (sseq(:,:) - avg).^2; %dist of each diff to the mean
    avgd            = mean(distseq, 1);% avg difficuly of each seq = cols
    ibest           = cl(:,:,(avgd == min(avgd)));%
        
    %4// store color values in struct colorvalues
    for seqval = 1:nseq
        stored = stored +1;
        colorvalues(stored).seqid = stored;
        colorvalues(stored).nsamples = nsmp;
        colorvalues(stored).seqnb = seqval;
        colorvalues(stored).colors = ibest(:,seqval);
        
    end%storing colvals into struct
    
end%for each possible seq length, do the simulation and sampling

fprintf('Done drawing color values!\n')
toc;

end%function def

 