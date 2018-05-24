%% run experiment
% TODO add screen close textures?
% TODO add colors to instructions?
% Q: give feedback on cond target? => no we want sampling to be rewarded??
% Q: give feedback on their nb of points?
% Q: cond open but both colors same => feedback always right??!

% Clear workspace & screens:
sca;
close all;
clearvars;
fprintf('Ready player one\n')

% Valentin's toolbox:
addpath ./Toolboxes/IO


%% Init subj: get nb, create folder, init resp matrix:

% query subj number in command line:
subj = input('>>>>>>>>>>>>>>>>>>>>>  Input participant nb: ');
if ~isscalar(subj) || mod(subj,1) ~= 0
    error('Invalid subject number!');
end

% make folder:
subjdir = sprintf('../Data/S%02d',subj);
if ~exist(subjdir,'dir')
    mkdir(subjdir);
end

% create header:
header = struct();
header.subj = subj;
header.start = datestr(now,'yyyymmdd-HHMM');

%create some form of matrix for the results
sampling = struct();
finalchoice = struct();

%% load / make sequences
% Generate sequences information
fprintf('Generating sequences for participant # %i...\n', subj)
sqc = gen_sequences(subj, 64);
%load('/Users/clemence/Dropbox/LNC/AGEXP/experiment/Data/S99/output99.mat', 'sqc');

%% init PTB, set visual parameters & make textures:

% Set visual parameters
ppd               = 40; % number of screen pixels per degree of visual angle
shape_sz          = 6*ppd; % ppd pixels per degree of visual angle = 40
shape_offset      = 4.5*ppd; % shape offset
respbutton_offset = 10*ppd; % response button offset
rep_sz            = 4*ppd;
fb_sz             = 6*ppd; % feedback square size
black             = [0,0,0];
white             = [1,1,1];
red               = [1,0,0];
blue              = [0,0,1];
green             = [0,1,0];
load('./color_space.mat', 'lbgd');

% Default PTB setup  + skipping tests on my laptop:
PsychDefaultSetup(2);
Screen('Preference', 'SkipSyncTests', 2);% //!\\ force skipping tests, PTB wont work
screens      = Screen('Screens');
screenNumber = min(screens); %get external screen
[window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
Screen('Flip', window);
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');%smooth blending/antialiasing
[screenXpixels, screenYpixels] = Screen('WindowSize', window);
[xCenter, yCenter] = RectCenter(windowRect);
Screen('TextSize', window, 50);

% Get frame duration and priority levels
ifi = Screen('GetFlipInterval', window);
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
isiSecs   = 0.5;
isiFrames = round(isiSecs / ifi);

% Load symbols & feedback shapes & make textures
% symbols to sample
shape_tex = zeros(2,8); % shape_tex(1,i)=black (contour) shape_tex(2,i)=grey (inside)
for i = 1:8
    outline = double(imread(sprintf('./img/shape%dc.png',i)))/255;
    outline = imresize(outline,shape_sz/size(outline,1));
    shape_tex(1,i) = Screen('MakeTexture',window,cat(3,ones(size(outline)),outline),[],[],2);
    inside = double(imread(sprintf('./img/shape%d.png',i)))/255;
    inside = imresize(inside,shape_sz/size(inside,1));
    shape_tex(2,i) = Screen('MakeTexture', window, cat(3,ones(size(inside)),inside),[],[],2);
end
positions(1,:)  = CenterRectOnPoint(Screen('Rect',shape_tex(1)),xCenter-shape_offset,yCenter);
positions(2,:)  = CenterRectOnPoint(Screen('Rect',shape_tex(1)),xCenter+shape_offset,yCenter);
centerpos     = CenterRectOnPoint(Screen('Rect',shape_tex(1)),xCenter,yCenter);

% resp buttons
outrep = double(imread('./img/shape7c.png'))/255;
outrep = imresize(outrep,rep_sz/size(outrep,1));
inrep = double(imread('./img/shape7.png'))/255;
inrep = imresize(inrep,rep_sz/size(inrep,1));
rep_tex(1) = Screen('MakeTexture',window,cat(3,ones(size(outrep)),outrep),[],[],2);
rep_tex(2) = Screen('MakeTexture',window,cat(3,ones(size(inrep)),inrep),[],[],2);
rep_pos(1,:)  = CenterRectOnPoint(Screen('Rect',rep_tex(1)),xCenter-respbutton_offset,screenYpixels*0.80);
rep_pos(2,:)  = CenterRectOnPoint(Screen('Rect',rep_tex(1)),xCenter+respbutton_offset,screenYpixels*0.80);

%feedbacks squares
fbsquare         = double(imread('./img/fbsquare.png'))/255;
fbsquare         = imresize(fbsquare,fb_sz/size(fbsquare,1));
fbsquare         = Screen('MakeTexture',window,cat(3,ones(size(fbsquare)),fbsquare),[],[],2);
fbpos(1,:)       = CenterRectOnPoint(Screen('Rect',fbsquare),xCenter-respbutton_offset,screenYpixels*0.80);
fbpos(2,:)       = CenterRectOnPoint(Screen('Rect',fbsquare),xCenter+respbutton_offset,screenYpixels*0.80);

%feedback labels
Screen('TextSize', window, 50);
Screen('FillRect', window, lbgd);
Screen('Flip',window);
[~, ~, textBounds] = DrawFormattedText(window, 'erreur !', 'center', 'center', red);
Screen('FillRect', window, lbgd);% text size, in black, recovered in gray
fblabelrec = ones(ceil((textBounds(4) - textBounds(2)) * 1.1),...
    ceil((textBounds(3) - textBounds(1)) * 1.1)).* lbgd; % rect for Texture of the sz of txt
fbcorr = Screen('MakeTexture', window, fblabelrec);
Screen('TextSize', fbcorr, 50);
DrawFormattedText(fbcorr, ' bravo !', 'center', 'center', green);
fberr = Screen('MakeTexture', window, fblabelrec);
Screen('TextSize', fberr, 50);
DrawFormattedText(fberr, ' erreur !', 'center', 'center', red);
fblabelpos(1,:)  = CenterRectOnPoint(Screen('Rect',fbcorr),xCenter-respbutton_offset, screenYpixels * 0.60);
fblabelpos(2,:)  = CenterRectOnPoint(Screen('Rect',fberr),xCenter+respbutton_offset, screenYpixels * 0.60);

% Keyboard keys
KbName('UnifyKeyNames'); %across OSs
spaceKey    = KbName('space');
escapeKey   = KbName('ESCAPE');
leftKey     = KbName('q');
rightKey    = KbName('p'); 
continueKey    = KbName('t'); 
respKey     = [leftKey,rightKey]; % left=1 right=2
RestrictKeysForKbCheck([spaceKey, escapeKey, leftKey, rightKey, continueKey]); %reset by clearall []
scanlist    = zeros(1,256); % speedup wait for key by restricting nb of keys to check
scanlist(respKey) = 1;

% set loop parameters (for testing)
progressON = true; % adding progress at each seq: "xx many percent of this bloc"
scoreON = false; % adding score at each sequence
forcepauseON = true; % forcing 30'' breaks
labelsON = true;
testingON = false;

if testingON
    nseq = 5; % nseq = numel(sqc);
    nsamples = 5;%sqc(iseq).nsamples;
    nseqperblock = 16; %actually blocs can be either 8 or 16 seq long...
    forcepauseON = false; % forcing 30'' breaks
else
    nseq = numel(sqc);
    nseqperblock = 16;
    %test kb:
    fprintf('Press [ESPACE] to check keyboard responsiveness...: ');
    if WaitKeyPress([],30) == 0
        fprintf('\n\n');
        error('No key press detected after 30 seconds.');
    else;fprintf('Good.\n\n');
    end
    % hide cursor and stop spilling key presses into MATLAB windows:
    HideCursor;
    FlushEvents;
    ListenChar(2);
end %setting testing parameters

%% loop over sequences and trials: TRY:
try
% welcome screen:
% Screen('FillRect', window, lbgd);
% Screen('Flip',window);
line1 = 'Début de l''expérience :\n\n\n';
line2 = '    - Il y aura 4 blocs de 16 séquences de symboles,\n';
line3 = '    - Des pauses sont prévues entre les bloc.\n    - etc...';
DrawFormattedText(window, [line1 line2 line3], screenYpixels * 0.25, screenYpixels * 0.25, black);
spaceline = 'Appuyez sur [espace] pour commencer.';
DrawFormattedText(window, spaceline, 'center', screenYpixels * 0.90, black);
Screen('Flip', window);
WaitKeyPress(spaceKey);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% loop over sequences
currenttrial = 1;
fprintf('Starting looping over trials\n');
currentseqinblock = 0;
for iseq = 1:nseq
%%%
    if mod(iseq,5) == 0; fprintf('iseq # %d/%d...\n', iseq, nseq); end

    % Get this sequence's data
    condition    = sqc(iseq).condition;
    targetsymbol = sqc(iseq).targetsymbol;
    othersymbol  = sqc(iseq).othersymbol;
    targetcolor  = sqc(iseq).targetcolor;
    if      targetcolor == 1; colorstr = 'bleu';  
    elseif  targetcolor == 2; colorstr = 'rouge';   
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Draw sequence instructions
    sidea  = randi(2,1); % get random sides of symbols for the instructions
    sideb  = 3 - sidea;
    if condition == 1
        line = sprintf('Voici les symboles, essayez de tirer du %s !', colorstr);
    elseif condition == 2
        line = 'Voici les symboles, devinez leur couleur !';
    end
    DrawFormattedText(window, line, 'center', screenYpixels * 0.30, black);
    Screen('DrawTexture', window, shape_tex(1,targetsymbol), [], positions(sidea,:), 0, [], [], [0,0,0]);
    Screen('DrawTexture', window, shape_tex(2,targetsymbol), [], positions(sidea,:), 0, [], [], lbgd);
    Screen('DrawTexture', window, shape_tex(1,othersymbol), [], positions(sideb,:), 0, [], [], [0,0,0]);
    Screen('DrawTexture', window, shape_tex(2,othersymbol), [], positions(sideb,:), 0, [], [], lbgd);
    
    if progressON
        currentseqinblock = currentseqinblock+1; %this lags by 1 element 0=>1
        spaceline = 'Appuyez sur [espace] pour lancer la séquence.';
        progressline = sprintf('Séquence %.f/%.f dans ce bloc, continuez !\n', currentseqinblock, nseqperblock);
        nextline = [progressline spaceline];
    elseif progressON == false
        spaceline = 'Appuyez sur [espace] pour lancer la séquence.';
        nextline = spaceline;
    end
    DrawFormattedText(window, nextline, 'center', screenYpixels * 0.90, black); 
    Screen('DrawingFinished', window);
    Screen('Flip', window);
    WaitSecs(0.5);
    WaitKeyPress(spaceKey);
    
    % Draw wait/blank frames for 0.5 secs
    Screen('DrawDots', window, [xCenter yCenter], 15 , black, [0 0], 2);
    Screen('Flip', window);
    WaitSecs(0.8);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%loop over trials
    nsamples = sqc(iseq).nsamples;
    %countingresponses = NaN(1,nsamples);
    for isample = 1:nsamples
        
        % Init resp var
         awaitingresp = true;

        % Check if abort key is pressed
        if CheckKeyPress(escapeKey); sca; break
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Draw sampling alternatives:
        targetside  = sqc(iseq).symbolsides(isample);
        otherside   = 3 - sqc(iseq).symbolsides(isample);
        Screen('DrawDots', window, [xCenter yCenter], 15 , black, [0 0], 2);
        Screen('DrawTexture', window, shape_tex(1,targetsymbol), [], positions(targetside,:), 0, [], [], [0,0,0]);
        Screen('DrawTexture', window, shape_tex(2,targetsymbol), [], positions(targetside,:), 0, [], [], lbgd);
        Screen('DrawTexture', window, shape_tex(1,othersymbol), [], positions(otherside,:) , 0, [], [], [0,0,0]);
        Screen('DrawTexture', window, shape_tex(2,othersymbol), [], positions(otherside,:) , 0, [], [], lbgd);
        Screen('DrawingFinished', window);
        Screen('Flip', window);
        vbl = GetSecs;
        
        % Get sampling choice:
        while awaitingresp == 1 
            [keyIsDown,secs,k] = KbCheck([],scanlist);
            if keyIsDown & find(k==1) == respKey(1) | find(k==1) == respKey(2)
                response = find(respKey == (find(k==1)));
                rt  = secs - vbl;
                awaitingresp = false;
            end
        end%while awaiting for resp
        
        % Draw color feedback for 0.5 secs
        samplecolor = sqc(iseq).rgb(isample,:,response);
        if response == targetside
            Screen('DrawDots', window, [xCenter yCenter], 15 , black, [0 0], 2);
            Screen('DrawTexture', window, shape_tex(1,targetsymbol), [], positions(targetside,:), 0, [], [], [0,0,0]);
            Screen('DrawTexture', window, shape_tex(2,targetsymbol), [], positions(targetside,:), 0, [], [], samplecolor);
        elseif response == otherside
            Screen('DrawDots', window, [xCenter yCenter], 15 , black, [0 0], 2);
            Screen('DrawTexture', window, shape_tex(1,othersymbol), [], positions(otherside,:) , 0, [], [], [0,0,0]);
            Screen('DrawTexture', window, shape_tex(2,othersymbol), [], positions(otherside,:) , 0, [], [], samplecolor);
        end
        Screen('Flip', window);
        WaitSecs(0.5);

        % Draw wait/blank frames for 0.5 secs
        Screen('DrawDots', window, [xCenter yCenter], 15 , black, [0 0], 2);
        Screen('Flip', window);
        WaitSecs(0.5);
        
        % Store results:
%         if response == targetside; countingresponses(isample) = 1;
%         elseif response == otherside; countingresponses(isample) = 2;
%         end
        sampling(currenttrial).subj = subj;
        sampling(currenttrial).condition = condition;
        sampling(currenttrial).mirroring = sqc(iseq).mirroring;
        sampling(currenttrial).sequencenb = iseq;
        sampling(currenttrial).samplenb = isample;
        sampling(currenttrial).nsamples = sqc(iseq).nsamples;
        sampling(currenttrial).trialnb = currenttrial;
        sampling(currenttrial).symbolpair = sqc(iseq).symbolpair;
        sampling(currenttrial).targetsym = targetsymbol;
        sampling(currenttrial).targetcolor = targetcolor;
        sampling(currenttrial).targetside = targetside;
        sampling(currenttrial).samplingchoice = response;
        sampling(currenttrial).samplingcolor = samplecolor;
        sampling(currenttrial).samplingrt = rt;
        
        currenttrial = currenttrial+1;
    end%FOR each sample
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Draw final probe while no answer yet:
    
    % Get sides of color buttons for final probe :
    bluebuttonside = sqc(iseq).bluebuttonside;
    redbuttonside  = 3-bluebuttonside;
    if targetcolor == 1%blue
        goodanswer = bluebuttonside;
        wronganswer = redbuttonside;
    elseif targetcolor == 2%red
        goodanswer = redbuttonside;
        wronganswer = bluebuttonside;
    end
    
    % Depending on condition, draw final probe:
    switch condition
        case 1
            line = sprintf('Dans cette séquence,\navez-vous eu l''impression de tirer plutôt du : ');
            DrawFormattedText(window,  line, 'center', screenYpixels*0.20, black);
        case 2
            line = sprintf('Dans cette séquence,\nde quelle couleur était le symbole suivant : ');
            DrawFormattedText(window,  line, 'center', screenYpixels*0.20, black);
            Screen('DrawTexture', window, shape_tex(1,targetsymbol), [], centerpos, 0, [], [], [0,0,0]);
            Screen('DrawTexture', window, shape_tex(2,targetsymbol), [], centerpos, 0, [], [], lbgd);
    end
    Screen('DrawTexture', window, rep_tex(1), [], rep_pos(redbuttonside,:), 0, [], [], [0,0,0]);
    Screen('DrawTexture', window, rep_tex(2), [], rep_pos(redbuttonside,:), 0, [], [], red);
    Screen('DrawTexture', window, rep_tex(1), [], rep_pos(bluebuttonside,:), 0, [], [], [0,0,0]);
    Screen('DrawTexture', window, rep_tex(2), [], rep_pos(bluebuttonside,:), 0, [], [], blue);
    Screen('DrawingFinished', window);
    Screen('Flip', window);
    vbl = GetSecs;
    
    % Get final response & and wait 0.5secs:
    awaitingfinalresp = 1;
    while awaitingfinalresp
        [keyIsDown,secs,k] = KbCheck([],scanlist);
        if keyIsDown & find(k==1) == respKey(1) | find(k==1) == respKey(2)
            respprobe = find(respKey == (find(k==1)));
            rtprobe  = secs - vbl;
            awaitingfinalresp = false;
        end
    end
    Screen('FillRect', window, lbgd);
    Screen('Flip', window);
    WaitSecs(0.3);
    
    % Feedback for final response IN CONDITION OPEN ONLY:
    if condition == 2
        % redraw probe
        line = sprintf('Dans cette séquence,\nde quelle couleur était le symbole suivant : ');
        DrawFormattedText(window,  line, 'center', screenYpixels*0.20, black);
        Screen('DrawTexture', window, shape_tex(1,targetsymbol), [], centerpos, 0, [], [], [0,0,0]);
        Screen('DrawTexture', window, shape_tex(2,targetsymbol), [], centerpos, 0, [], [], lbgd);
        Screen('DrawTexture', window, rep_tex(1), [], rep_pos(redbuttonside,:), 0, [], [], [0,0,0]);
        Screen('DrawTexture', window, rep_tex(2), [], rep_pos(redbuttonside,:), 0, [], [], red);
        Screen('DrawTexture', window, rep_tex(1), [], rep_pos(bluebuttonside,:), 0, [], [], [0,0,0]);
        Screen('DrawTexture', window, rep_tex(2), [], rep_pos(bluebuttonside,:), 0, [], [], blue);
        
        % add feedback squares and labels
        switch respprobe 
            case goodanswer
                Screen('DrawTexture', window, fbsquare, [], fbpos(respprobe,:), [], [], [], green); %square
                if labelsON == 1; Screen('DrawTextures', window, fbcorr,[], fblabelpos(respprobe,:), [], []);
                end
            case wronganswer
                Screen('DrawTexture', window, fbsquare, [], fbpos(respprobe,:), [], [], [], red); %square
                if labelsON == 1; Screen('DrawTextures', window, fberr,[], fblabelpos(respprobe,:), [], []);
                end
        end
        Screen('DrawingFinished', window);
        Screen('Flip', window);
        WaitSecs(0.8);
        
    elseif condition == 1
        Screen('FillRect', window, lbgd);
        Screen('Flip', window);
        WaitSecs(0.8);
    end%IF cond1, show feedback
    
    % Store final response: cuisine pour remettre  info dans tous les samples: TODO?
    allthesesamples = (currenttrial -(nsamples)) : currenttrial-1;
    for i = currenttrial-1:-1:(currenttrial -(nsamples))
        sampling(i).finaltargetside = goodanswer;
        sampling(i).finalresponse = respprobe;
        sampling(i).finalrt = rtprobe;
    end
    finalchoice(iseq).subj = subj;
    finalchoice(iseq).condition = condition;
    finalchoice(iseq).mirroring = sqc(iseq).mirroring;
    finalchoice(iseq).sequencenb = iseq;
    finalchoice(iseq).nsamples = sqc(iseq).nsamples;
    finalchoice(iseq).symbolpair = sqc(iseq).symbolpair;
    finalchoice(iseq).targetsymbol = targetsymbol;
    finalchoice(iseq).targetcolor = targetcolor;
    finalchoice(iseq).goodanswerside  = goodanswer;
    finalchoice(iseq).finalresponse    = respprobe;
    finalchoice(iseq).finalrt          = rtprobe;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % If this is a block, insert pause & give score
    presSecs = 30:-1:0;
    if mod(iseq,nseqperblock) == 0
        fprintf('Participant is on their break!\nPress [space] and [t] to continue...\n')
        if forcepauseON
        for sec = 1:numel(presSecs)
            line1 = 'C''est la fin de ce bloc, bravo !\n';
            DrawFormattedText(window, line1, 'center', screenYpixels * 0.25, black);
            line2 = sprintf('PRENEZ UNE PAUSE !\nAu moins 30 secondes ; vous pouvez prendre plus.\n%d', presSecs(sec));
            DrawFormattedText(window, line2, 'center', screenYpixels * 0.50, black);
%             line3 = num2str(presSecs(sec));
%             DrawFormattedText(window, line3, 'center', screenYpixels * 0.50, black);
            Screen('Flip', window);
            WaitSecs(1);
        end
        end
        line1 = 'C''est la fin de ce bloc, bravo !\n';
        DrawFormattedText(window, line1, 'center', screenYpixels * 0.25, black);
        line2 = 'PRENEZ UNE PAUSE !\n\n';
        DrawFormattedText(window, line2, 'center', screenYpixels * 0.50, black);
        line3 = 'Quand vous êtes prêt-e à continuer,\nallez chercher l''expérimentateur.';
        DrawFormattedText(window, line3, 'center', screenYpixels * 0.75, black);
        Screen('DrawingFinished', window);
        Screen('Flip', window);
        WaitKeyPress(spaceKey);
        WaitKeyPress(continueKey);
        currentseqinblock = 0;
    end%if this is a bloc, insert pause
    
%%%
end%FOR each sequence of samples
%%%

% draw last screen
line = 'Merci d''avoir participé, l''expérience est finie !\n';
endline = 'TERMINER -- [espace] -- TERMINER';
DrawFormattedText(window, [line endline], 'center', yCenter, black);
Screen('DrawingFinished', window);
Screen('Flip', window);
WaitKeyPress(spaceKey);
sca;

% save results, sqc and header in single mat file in subj directory
header.end = datestr(now,'yyyymmdd-HHMM');
save([subjdir '/' 'output' num2str(subj)], 'sampling', 'finalchoice', 'sqc', 'header')

%closing
Priority(0);
FlushEvents;
ListenChar(0);
ShowCursor;
Screen('Preference', 'SkipSyncTests', 0);%restore synch testing
Screen('Preference','VisualDebuglevel',4);%restore warning verbosity
Screen('Preference','SuppressAllWarnings',0);%=1 == warning verbosity=0


%% catching results so far if crashes
catch ME
    % save sofar:
    header.end = datestr(now,'yyyymmdd-HHMM');
    header.aborted = 1;
    save([subjdir '/' 'output_sofar' num2str(subj)], 'sampling', 'finalchoice', 'sqc', 'header')

    % clean close Psychtoolbox
    Priority(0);
    Screen('CloseAll');
    FlushEvents;
    ListenChar(0);
    ShowCursor;
    
    % handle error
    rethrow(ME);
end% TRY/CATCH
