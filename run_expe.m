%% run experiment
% NB: native KB keys are not inhibited: don't press any key!

% iblock = block idx in nb of blocks
% iseq = seq idx in nb of seq in this block
% isqc = seq idx out of total nb of sequences
% isample = sample idx in nb of samples in this seq
% itotal = sample idx out of total nb of samples

%% Clear workspace & screens, addpaths: 
sca;
close all;
clearvars;

addpath ./Toolboxes/IO % Valentin's toolbox in Experiment's path
%addpath C:/Toolboxes/Psychtoolbox/ % PTB path specific to laptop
%addpath('/Applications/Psychtoolbox/PsychHardware/EyelinkToolbox/'); % PTB toolbox to interface w Eyelink
%addpath(genpath('/Applications/Eyelink/'));

%% Init subj: get nb, create folder, init resp matrix, load or generate sequences: 

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

% Generate sequences information
generate = input('>>>>>>>>>>>>>>>>>>>>>  Generate sequences? 1=yes anythingelse=no ');
if generate == 1
    fprintf('Generating sequences for participant # %i...\n', subj)
    sqc = gen_sequences(subj, 96);
    sqc = fx_makingblocks(sqc);
    sqc = fx_orderingblocks(sqc);
    save(sprintf('../Data/S%02d/sqc%02d.mat', subj, subj), 'sqc');
else
    fprintf('Loading existing sqc structure from ''/Data/S%02d/sqc%02d.mat''...\n', subj, subj)
    load(sprintf('../Data/S%02d/sqc%02d.mat', subj, subj),'sqc');
end

%% Running parameters: 
video.testingON = false;
video.synchflip = true;
video.eyeON = true;
video.mouseON = false;
video.eyelink.inputphys  = [530 295]; %measured => this gave the best calibrations in the end.
%video.eyelink.inputphys  = [400 300]; %old small screen
video.eyelink.inputscreendistance_mm = 660;

%% PTB parameters: 

%%% Set visual parameters
video.ppd               = 40; % number of screen pixels per degree of visual angle
video.shape_sz          = 6*video.ppd; % video.ppd pixels per degree of visual angle = 40
video.shape_offset      = 4.5*video.ppd; % shape offset
video.respbutton_offset = 10*video.ppd; % response button offset
video.rep_sz            = 5*video.ppd;
video.fb_sz             = 6*video.ppd; % feedback square size
black                   = [0,0,0];
white                   = [1,1,1];
red                     = [1,0,0];
green                   = [0,1,0];
load('./color_space.mat', 'lbgd');
video.lbgd              = lbgd;
orange                  = [1,0.5,0];
blue                    = [0,0,1];

%%% which screens?
screens = Screen('Screens');
video.screen = max(screens); %screen index on PC: 0=extended, 1=native, 2=ext

%%% set synch properties:
if video.synchflip && ispc
    % set screen synchronization properties -- workaround Valentin for PC
    % see 'help SyncTrouble',
    %     'help BeampositionQueries' or
    %     'help ConserveVRAMSettings' for more information
    Screen('Preference','VisualDebuglevel',3); % verbosity
    Screen('Preference','SyncTestSettings',[],[],0.2,10); % soften synchronization test requirements
    Screen('Preference','ConserveVRAM',bitor(4096,Screen('Preference','ConserveVRAM'))); % enforce beamposition workaround for missing VBL interval
    fprintf('Synching flips with softer requirements...\n');
elseif ~video.synchflip || ismac
    % skip synchronization tests altogether
    % //!\\ force skipping tests, PTB wont work = timing inaccurate
    Screen('Preference','SkipSyncTests',2); % assumes 60Hz etc..
    Screen('Preference','VisualDebuglevel',0);
    Screen('Preference','SuppressAllWarnings',1);
    fprintf('||| SYNCHFLIP OFF or running on OSX => TIMINGS WILL BE INACCURATE! |||\n\n')
end

%%% open main window:
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask','General','UseFastOffscreenWindows');
    PsychImaging('AddTask','General','NormalizedHighresColorRange');
    video.res = Screen('Resolution',video.screen);
    
%     % set res manually..?
%     res  = []; % screen resolution, [] = full screen
%     fps  = []; % screen refresh rate
%     if ~isempty(res) && ~isempty(fps) % set screen res and refresh rate
%     r = Screen('Resolutions',screennb);
%     i = find([r.width] == res(1) & [r.height] == res(2));
%     if isempty(i) || ~any([r(i).hz] == fps)
%         error('Cannot set screen to %d x %d at %d Hz.',res(1),res(2),fps);
%     end
%     Screen('Resolution',video.screen, res(1), res(2), fps);
%     end
    
    fprintf('before openwindow\n\n')
    [video.window, video.windowRect] = PsychImaging('OpenWindow',video.screen,video.lbgd);
    fprintf('passed OpenWindow!\n\n')

    [video.screenXpixels, video.screenYpixels] = Screen('WindowSize', video.window);
    [video.xCenter, video.yCenter] = RectCenter(video.windowRect);
    Screen('BlendFunction', video.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');%smooth blending/antialiasing
    video.blend = '''GL_SRC_ALPHA'', ''GL_ONE_MINUS_SRC_ALPHA''';
    video.colorrange = 1;
    Screen('ColorRange',video.window,video.colorrange);
    % 1 = to pass color values in OpenGL's native floating point color range of 0.0 to
    % 1.0: This has two advantages: First, your color values are independent of
    % display device depth, i.e. no need to rewrite your code when running it on
    % higher resolution hardware. Second, PTB can skip any color range remapping
    % operations - this can speed up drawing significantly in some cases.   ...??? 
    Screen('TextSize', video.window, 50);
    video.textsize = 50;
    video.ifi = Screen('GetFlipInterval', video.window);
    video.priority = Priority(MaxPriority(video.window));
    Screen('Flip', video.window);

%% Make textures: 

%%% symbols textures:
shape_tex = zeros(2,8); % shape_tex(1,i)=black (contour) shape_tex(2,i)=grey (inside)
for sh = 1:8
    outline = double(imread(sprintf('./img/shape%dc.png',sh)))/255;
    outline = imresize(outline,video.shape_sz/size(outline,1));
    shape_tex(1,sh) = Screen('MakeTexture',video.window,cat(3,ones(size(outline)),outline),[],[],2);
    inside = double(imread(sprintf('./img/shape%d.png',sh)))/255;
    inside = imresize(inside,video.shape_sz/size(inside,1));
    shape_tex(2,sh) = Screen('MakeTexture', video.window, cat(3,ones(size(inside)),inside),[],[],2);
end
positions(1,:)  = CenterRectOnPoint(Screen('Rect',shape_tex(1)),video.xCenter-video.shape_offset,video.yCenter);
positions(2,:)  = CenterRectOnPoint(Screen('Rect',shape_tex(1)),video.xCenter+video.shape_offset,video.yCenter);
centerpos     = CenterRectOnPoint(Screen('Rect',shape_tex(1)),video.xCenter,video.yCenter);

% resp buttons textures:
outrep = double(imread('./img/bulbout.png'))/255;
outrep = imresize(outrep,video.rep_sz/size(outrep,1));
inrep = double(imread('./img/bulbin.png'))/255;
inrep = imresize(inrep,video.rep_sz/size(inrep,1));
rep_tex(1,1) = Screen('MakeTexture',video.window,cat(3,ones(size(outrep)),outrep),[],[],2);
rep_tex(2,1) = Screen('MakeTexture',video.window,cat(3,ones(size(inrep)),inrep),[],[],2);
rep_tex(1,2) = Screen('MakeTexture',video.window,cat(3,ones(size(outrep)),flip(outrep,2)),[],[],2);
rep_tex(2,2) = Screen('MakeTexture',video.window,cat(3,ones(size(inrep)),flip(inrep,2)),[],[],2);
rep_pos(1,:)  = CenterRectOnPoint(Screen('Rect',rep_tex(1)),video.xCenter-video.respbutton_offset,video.screenYpixels*0.80);
rep_pos(2,:)  = CenterRectOnPoint(Screen('Rect',rep_tex(1)),video.xCenter+video.respbutton_offset,video.screenYpixels*0.80);

%feedbacks squares textures:
fbsquare         = double(imread('./img/fbsquare.png'))/255;
fbsquare         = imresize(fbsquare,video.fb_sz/size(fbsquare,1));
fbsquare         = Screen('MakeTexture',video.window,cat(3,ones(size(fbsquare)),fbsquare),[],[],2);
fbpos(1,:)       = CenterRectOnPoint(Screen('Rect',fbsquare),video.xCenter-video.respbutton_offset,video.screenYpixels*0.80);
fbpos(2,:)       = CenterRectOnPoint(Screen('Rect',fbsquare),video.xCenter+video.respbutton_offset,video.screenYpixels*0.80);

%% Init keys: 
clear PsychHID; % Force new enumeration of devices.
clear KbCheck;
KbName('UnifyKeyNames'); %across OSs
GetKeyboardIndices();

escapeKey   = KbName('ESCAPE'); %experimenter keys (native Kb)
continueKey = KbName('y');

spaceKey    = KbName('space'); %participant keys (numpad)
if ispc; leftKey = KbName('a');
elseif ismac; leftKey = KbName('q');
end
rightKey    = KbName('p');

respKeys     = [leftKey, rightKey];% left=1 right=2
expeKeys     = [escapeKey, continueKey, leftKey, rightKey, spaceKey];
allKeys      = [];
k = RestrictKeysForKbCheck(allKeys);

%% Init EYELINK: 

if video.eyeON
    % Check the connection is up and running:
    ip = '100.1.1.2';
    if ispc
        system(sprintf('ping %s -n 3', ip));
    elseif ismac
        system(sprintf('ping %s -c 5', ip));
    end
    
    % Initialize Eyelink computer:
    status = Eyelink('Initialize', 'PsychEyelinkDispatchCallback'); % ZERO IF OK
    if status ~= 0; error('could not initialize eye-tracker connection!'); end
    %status = Eyelink('InitializeDummy');
    [~,ev] = Eyelink('GetTrackerVersion');
    connected = Eyelink('IsConnected'); %1=ok, 2=broadcast, -1=dummy, 0=none
    fprintf('Connection to %s eye-tracker initialized = %d.\n',ev, connected);
    
    % INIT PTB EL
    el = EyelinkInitDefaults(video.window);
    %el.devicenumber = []; %see KbCheck for details of this value
    
    % UPDATE SET OPTIONS SCREEN:
        % calibration & validation options:
    Eyelink('Command', 'automatic_calibration_pacing = 1000');
    Eyelink('Command', 'randomize_calibration_order = YES');
    Eyelink('Command', 'enable_automatic_calibration = YES'); % force manual accept
    Eyelink('Command', 'select_eye_after_validation = NO');
    Eyelink('Command', 'cal_repeat_first_target = YES');
    Eyelink('Command', 'val_repeat_first_target = YES');
        
        % tracking options:
    %Eyelink('Command', 'search_limits = YES');
    %Eyelink('Command', 'move_limits = YES');
    %Eyelink('Command', 'mouse_simulation = YES');
    Eyelink('Command', 'pupil_size_diameter = NO');
        
        % events & data processing options:    - this is the cognitive configuration
    Eyelink('Command','select_parser_configuration = 0'); %this is the cognitive configuration
    Eyelink('Command','recording_parse_type = GAZE');
    Eyelink('Command','saccade_velocity_threshold = 30');    
    Eyelink('Command','saccade_acceleration_threshold = 8000');    
    Eyelink('Command','saccade_motion_threshold = 0.1'); 
    Eyelink('Command','saccade_pursuit_fixup = 60');    
    Eyelink('Command','fixation_update_interval = 0');  
    %Eyelink('Command', 'saccade_sensitivity = NORMAL');
    %Eyelink('Command', 'file_sample_filter = EXTRA');
    %Eyelink('Command','Analog_filter = STD');
    
        % analog output options:
    %Eyelink('Command', 'analog_output = 0');
    
        % File contents:   - defaults from manual adapted to my task:
    Eyelink('Command','file_event_data = GAZE, GAZERES, AREA, HREF, VELOCITY');
    Eyelink('Command','file_event_filter = LEFT, RIGHT, FIXATION, SACCADE, BLINK, MESSAGE');
    Eyelink('Command','file_sample_data = LEFT,RIGHT,GAZE,GAZERES,AREA,STATUS');       
    
        % UPDATE CAMERA SETUP SCREEN:
        % Tracking mode:
    Eyelink('Command', 'pupil_size_diameter = YES'); % = Valentin does diameters...
    Eyelink('Command', 'use_ellipse_fitter = NO'); % = centroid
    Eyelink('Command', 'sample_rate = 1000'); % Valentin does 500..? could be 1000???
    Eyelink('Command', 'elcl_tt_power = 2'); % 75%  by default
    Eyelink('Command', 'active_eye = LEFT'); % in front of camera (less noise)
    Eyelink('Command','binocular_enabled = NO'); % monoc
    Eyelink('Command','corneal_mode = YES'); % pupilCR
    Eyelink('Command','use_high_speed = YES'); % defines sampling rate... > redundant
    
            % SCREEN PARAMETERS (physical and pixels):
    Eyelink('Command', 'simulation_screen_distance = %d', video.eyelink.inputscreendistance_mm);
    video.eyelink.pixelcoords = {video.windowRect(1), video.windowRect(2), video.windowRect(3)-1, video.windowRect(4)-1}';
    Eyelink('Command', 'screen_pixel_coords = %d, %d, %d, %d', video.eyelink.pixelcoords{:});
    video.eyelink.physicalscreensize = {-round(video.eyelink.inputphys(1)/2) round(video.eyelink.inputphys(2)/2) round(video.eyelink.inputphys(1)/2) -round(video.eyelink.inputphys(2)/2)}';
    Eyelink('Command', 'screen_phys_coords = %d, %d, %d, %d', video.eyelink.physicalscreensize{:});
    % NB only update calibration parameters after display parameters
    % call calibration type afterwards
    
            % CALIBRATION PARAMETERS VALENTIN:
    el.backgroundcolour = video.lbgd;
    Eyelink('Command','calibration_type = HV5'); % HV5 or HV9
    Eyelink('Command','generate_default_targets = NO'); % YES:default or NO:custom
    cnt = [video.screenXpixels/2,video.screenYpixels/2];
    off = video.ppd*9; %off = video.ppd*8;
    pnt = zeros(5,2);
    pnt(1,:) = cnt;
    pnt(2,:) = cnt-[0,off];
    pnt(3,:) = cnt+[0,off];
    pnt(4,:) = cnt-[off,0];
    pnt(5,:) = cnt+[off,0];
    pnt = num2cell(reshape(pnt',[],1));
    %calib
    Eyelink('Command','calibration_samples = 6');
    Eyelink('Command','calibration_sequence = 0,1,2,3,4,5');
    Eyelink('Command','calibration_targets = %d,%d %d,%d %d,%d %d,%d %d,%d',pnt{:});
    %valid
    Eyelink('Command','validation_samples = 5');
    Eyelink('Command','validation_sequence = 0,1,2,3,4,5');
    Eyelink('Command','validation_targets = %d,%d %d,%d %d,%d %d,%d %d,%d',pnt{:});
    
    Eyelink('Command','calibration_type = HV5'); % HV5 or HV9
    
    el.displayCalResults = 1;
    EyelinkUpdateDefaults(el); % pass the values back to the Eyelink
    
end

%% TRY/CATCH/END - start display

try
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% test kb, hide cursor, stop spilling keys
    if ~video.testingON
        fprintf('Press [ESPACE] to check keyboard responsiveness...: ');
        if WaitKeyPress([],30) == 0
            fprintf('\n\n');
            error('No key press detected after 30 seconds.');
        else;fprintf('Good.\n\n');
        end
        HideCursor;
        FlushEvents;
        ListenChar(2);
    end %setting testing parameters
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% welcome screen:
    line1 = 'Debut de l''experience :\n\n\n';
    line2 = '    - Il y aura 4 blocs de 16 sequences de symboles,\n';
    line3 = '    - Des pauses sont prevues entre les blocs.\n    - D''abord : calibration de l''eyetracker';
    DrawFormattedText(video.window, [line1 line2 line3], video.screenYpixels * 0.25, video.screenYpixels * 0.25, black);
    spaceline = 'Appuyez sur [espace] pour commencer.';
    DrawFormattedText(video.window, spaceline, 'center', video.screenYpixels * 0.90, black);
    Screen('Flip', video.window);
    WaitKeyPress(spaceKey);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% loop over blocks
    fprintf('Starting looping over blocks... \n\n');
    blocks = unique([sqc.blocknb]);
    nblocks = numel(blocks);
    
    itotal = 0;
    for iblock = 1:nblocks
        
        thissqc = sqc([sqc.blocknb] == iblock);
        nseq = numel(thissqc);
        blockscores = NaN(1, nseq);
        
        % EYELINK: Calibrate & lock eye after calibration
        % EYELINK: Open EDF file, write header
        % EYELINK: Send SYNCHTIME trigger
        if video.eyeON
            
            k = RestrictKeysForKbCheck(allKeys);
            
            fprintf('Calibrate eyetracker for block %d\n', iblock);
            %%% calibrate:
            calibration = EyelinkDoTrackerSetup(el, el.ENTER_KEY); %el.ENTER_KEY shows the eye image on display screen
            [~, calmessage] = Eyelink('CalMessage');
            fprintf('\nCalibration resulted in %s \n', calmessage); %TODO what's this?
            
            %%% start file:
            edfname = sprintf('COX%03d_%d', subj, iblock);
            status = Eyelink('OpenFile', edfname);
            fprintf('\nOpening EDF for block %d returns: %d (zero is ok!) \n', iblock, status);
            Eyelink('Command', 'set_idle_mode');
            Eyelink('Message', 'DISPLAY_COORDS %ld %ld %ld %ld', video.eyelink.pixelcoords{:});
            Eyelink('Message', 'PHY_SCREEN_SIZE %ld %ld %ld %ld', video.eyelink.physicalscreensize{:});
            Eyelink('Message', 'PHY_SCREEN_DISTANCE = %d', video.eyelink.inputscreendistance_mm);
            Eyelink('Message',calmessage);
            Eyelink('command', 'draw_cross %d %d 1', video.xCenter, video.yCenter); % will stay as long as not changed...
            WaitSecs(0.1);
            Eyelink('StartRecording');
            WaitSecs(1); % record 1sec -> preproc Valentin removes 1'' before & after block    
        end
        
        k = RestrictKeysForKbCheck(expeKeys);
        
        % Begin block screen:
        fprintf('>>>Start block %.f/%.f----- %s<<<\n', iblock, nblocks, datestr(now, 'HH:MM:SS'));
        line = sprintf('DEBUT DU BLOC %.f/%.f !\n', iblock, nblocks);
        DrawFormattedText(video.window, line, 'center', video.screenYpixels * 0.30, black);
        spaceline = 'Appuyez sur [espace] pour commencer.';
        DrawFormattedText(video.window, spaceline, 'center', video.screenYpixels * 0.90, black);
        Screen('Flip', video.window);
        WaitKeyPress(spaceKey);
        
        % EYELINK: send synchtime message
        if video.eyeON
            Eyelink('Message', 'SYNCHTIME'); %TODO IS THAT NECESSARY???
        end
        
        if video.mouseON
            fprintf('mouseON break');
            closeEL(eyeON, subj, iblock);
            closePTB();
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% loop over sequences
        fprintf('Starting looping over sequences\n');
        for iseq = 1:nseq
            
            % Get current seq index in sqc
            isqc = find([sqc.blocknb] == iblock & [sqc.blkseqnb] == iseq);
            if mod(isqc,5) == 0; fprintf('isqc # %d/%d...\n', isqc, numel(sqc)); end
            
            % Draw sequence instructions
            targetsymbol = sqc(isqc).targetsymbol;
            othersymbol  = sqc(isqc).othersymbol;
            condition    = sqc(isqc).condition;
            targetcolor  = sqc(isqc).targetcolor;
            sidea  = randi(2,1); % get random sides of symbols for the instructions
            sideb  = 3 - sidea;
            if      targetcolor == 1; colorstr = 'bleu'; % get target color string
            elseif  targetcolor == 2; colorstr = 'orange';
            end
            if condition == 1 % get instructions string
                line = sprintf('Voici les symboles, essayez de tirer du %s !', colorstr);
            elseif condition == 2
                line = 'Voici les symboles, devinez leur couleur !';
            end
            DrawFormattedText(video.window, line, 'center', video.screenYpixels * 0.30, black);
            
            % Draw sequence symbols
            Screen('DrawDots', video.window, [video.xCenter video.yCenter], 15 , black, [0 0], 2);
            Screen('DrawTexture', video.window, shape_tex(1,targetsymbol), [], positions(sidea,:), 0, [], [], [0,0,0]);
            Screen('DrawTexture', video.window, shape_tex(2,targetsymbol), [], positions(sidea,:), 0, [], [], video.lbgd);
            Screen('DrawTexture', video.window, shape_tex(1,othersymbol), [], positions(sideb,:), 0, [], [], [0,0,0]);
            Screen('DrawTexture', video.window, shape_tex(2,othersymbol), [], positions(sideb,:), 0, [], [], video.lbgd);
            
            % Draw sequence progress at the bottom
            spaceline = 'Appuyez sur [espace] pour lancer la sequence.';
            progressline = sprintf('Sequence %.f/%.f dans ce bloc, continuez !\n', iseq, nseq);
            nextline = [progressline spaceline];
            DrawFormattedText(video.window, nextline, 'center', video.screenYpixels * 0.90, black);
            
            % Flip, wait for spacebar press
            Screen('DrawingFinished', video.window);
            Screen('Flip', video.window);
            WaitSecs(0.5);
            k = RestrictKeysForKbCheck(spaceKey);
            WaitKeyPress(spaceKey);
            
            % EYELINK: check recording
            % EYELINK: send begin sequence x message
            if video.eyeON
                error=Eyelink('checkrecording');
                if(error~=0)
                    error('checkrecording problem')
                end
                msg = sprintf('ISQC_instr_%03d', isqc);
                Eyelink('Message', msg);
            end
            
            % Flip fixation point for 1sec before starting
            Screen('DrawDots', video.window, [video.xCenter video.yCenter], 15 , black, [0 0], 2);
            Screen('Flip', video.window);
            WaitSecs(1);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% loop over samples
            nsamples = sqc(isqc).nsamples;
            seqpoints = NaN(1,nsamples);
            
            for isample = 1:nsamples
                
                k = RestrictKeysForKbCheck(expeKeys); %reset by clearall or RestrictKeysForKbCheck([])
                
                % Check if abort key is pressed
                if CheckKeyPress(escapeKey)
                    Priority(0);
                    FlushEvents;
                    ListenChar(0);
                    ShowCursor;
                    closeEL(eyeON, subj, iblock);
                    closePTB();
                    error('aborted by ESCAPE');
                end
                                
                % Draw sampling alternatives: draw "target" on "symbolside"
                targetside  = sqc(isqc).symbolsides(isample);
                otherside   = 3 - sqc(isqc).symbolsides(isample);
                Screen('DrawDots', video.window, [video.xCenter video.yCenter], 15 , black, [0 0], 2);
                Screen('DrawTexture', video.window, shape_tex(1,targetsymbol), [], positions(targetside,:), 0, [], [], [0,0,0]);
                Screen('DrawTexture', video.window, shape_tex(2,targetsymbol), [], positions(targetside,:), 0, [], [], video.lbgd);
                Screen('DrawTexture', video.window, shape_tex(1,othersymbol), [], positions(otherside,:) , 0, [], [], [0,0,0]);
                Screen('DrawTexture', video.window, shape_tex(2,othersymbol), [], positions(otherside,:) , 0, [], [], video.lbgd);
                Screen('DrawingFinished', video.window);
                vbldisplaytime = Screen('Flip', video.window);
                
                % EYELINK: flip sampling probe x
                if video.eyeON
                    msg = sprintf('ISMP_smpprobe_seq%03d_smp%02d', isqc, isample);
                    Eyelink('Message', msg);
                end
                
                % EYELINK: draw fixation cross on host PC??? TODO
                %Eyelink('Command','draw_cross', 'x', 'y', '15'); %'color (1 to 15));
                %Eyelink('Command','draw_cross', 'x', 'y', '15'); %'color (1 to 15));
                %PCEyelink('command', 'draw_cross %d %d 15', winWidth/2,winHeight/2);
                
                % Get sampling choice:
                while true
                    [keyIsDown,secs,keylist] = KbCheck();
                    if keyIsDown && sum(keylist) == 1 && (keylist(leftKey) == 1 || keylist(rightKey) == 1)
                        % EYELINK: resp sampling x
                        if video.eyeON
                            msg = sprintf('ISMP_resp_seq%03d_smp%02d', isqc, isample);
                            Eyelink('Message', msg);
                        end
                        rt = secs - vbldisplaytime;
                        response = find(respKeys == (find(keylist==1)));
                        break
                    end
                end%while awaiting for resp
                
                
                % Draw color feedback for 0.5 secs: sqc(isqc).rgb(1)= target // sqc(isqc).rgb(2)= other
                if response == targetside
                    samplecolor = sqc(isqc).rgb(isample,:,1); % 1 is target line
                    Screen('DrawDots', video.window, [video.xCenter video.yCenter], 15 , black, [0 0], 2);
                    Screen('DrawTexture', video.window, shape_tex(1,targetsymbol), [], positions(targetside,:), 0, [], [], [0,0,0]);
                    Screen('DrawTexture', video.window, shape_tex(2,targetsymbol), [], positions(targetside,:), 0, [], [], samplecolor);
                    seqpoints(isample) = sqc(isqc).colorslin(1,isample);
                elseif response == otherside
                    samplecolor = sqc(isqc).rgb(isample,:,2); % 2 is other color
                    Screen('DrawDots', video.window, [video.xCenter video.yCenter], 15 , black, [0 0], 2);
                    Screen('DrawTexture', video.window, shape_tex(1,othersymbol), [], positions(otherside,:) , 0, [], [], [0,0,0]);
                    Screen('DrawTexture', video.window, shape_tex(2,othersymbol), [], positions(otherside,:) , 0, [], [], samplecolor);
                    seqpoints(isample) = sqc(isqc).colorslin(2,isample);
                end
                Screen('Flip', video.window, secs+roundfp(0.1,0,video.ifi));
                
                % EYELINK: fb sampling x
                if video.eyeON
                    msg = sprintf('ISMP_colfb_seq%03d_smp%02d', isqc, isample);
                    Eyelink('Message', msg);
                end
                WaitSecs(0.5);
                
                % Draw fixation point for 0.5 secs
                Screen('DrawDots', video.window, [video.xCenter video.yCenter], 15 , black, [0 0], 2);
                Screen('Flip', video.window);
                WaitSecs(0.5);
                
                % Store results:
                itotal = itotal +1;
                sampling(itotal).subj = subj;
                sampling(itotal).condition = condition;
                sampling(itotal).mirroring = sqc(isqc).mirroring;
                sampling(itotal).blocknb = sqc(isqc).blocknb;
                sampling(itotal).blkseqnb = sqc(isqc).blkseqnb;
                sampling(itotal).sequencenb = iseq;
                sampling(itotal).samplenb = isample;
                sampling(itotal).nsamples = sqc(isqc).nsamples;
                sampling(itotal).trialnb = isqc;
                sampling(itotal).symbolpair = sqc(isqc).symbolpair;
                sampling(itotal).targetsym = targetsymbol;
                
                sampling(itotal).coloragency = sqc(isqc).coloragency;
                sampling(itotal).targetcolor = targetcolor;
                
                sampling(itotal).targetside = targetside;
                sampling(itotal).samplingchoice = response;
                sampling(itotal).samplingcolor = samplecolor;
                sampling(itotal).samplingrt = rt;
                
            end%FOR each sample
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            k = RestrictKeysForKbCheck(expeKeys);
            
            % End of each sequence: draw final probe while no answer yet:
            % Final probe color buttons
            goodanswer = sqc(isqc).probeside;% this is the side of the target color
            wronganswer = 3-goodanswer;
            if targetcolor == 1%blue
                bluebuttonside = goodanswer;
                redbuttonside = 3-bluebuttonside;
            elseif targetcolor == 2%red
                bluebuttonside = wronganswer;
                redbuttonside = 3-bluebuttonside;
            end
            
            % Final probe text; depending on condition,
            switch condition
                case 1
                    line = sprintf('Dans cette sequence,\navez-vous eu l''impression de tirer plutot du : ');
                    DrawFormattedText(video.window,  line, 'center', video.screenYpixels*0.20, black);
                case 2
                    line = sprintf('Dans cette sequence,\nde quelle couleur etait le symbole suivant : ');
                    DrawFormattedText(video.window,  line, 'center', video.screenYpixels*0.20, black);
                    Screen('DrawTexture', video.window, shape_tex(1,targetsymbol), [], centerpos, 0, [], [], [0,0,0]);
                    Screen('DrawTexture', video.window, shape_tex(2,targetsymbol), [], centerpos, 0, [], [], video.lbgd);
            end
            Screen('DrawTexture', video.window, rep_tex(1,redbuttonside), [], rep_pos(redbuttonside,:), 0, [], [], [0,0,0]);
            Screen('DrawTexture', video.window, rep_tex(2,redbuttonside), [], rep_pos(redbuttonside,:), 0, [], [], orange);
            Screen('DrawTexture', video.window, rep_tex(1, bluebuttonside), [], rep_pos(bluebuttonside,:), 0, [], [], [0,0,0]);
            Screen('DrawTexture', video.window, rep_tex(2, bluebuttonside), [], rep_pos(bluebuttonside,:), 0, [], [], blue);
            Screen('DrawingFinished', video.window);
            vbldisplaytime = Screen('Flip', video.window);
            
            % EYELINK: final probe seq x
            if video.eyeON
                msg = sprintf('ISQC_finalprobe_%03d', isqc);
                Eyelink('Message', msg);
            end
            
            % Wait for final response:
            while true
                [keyIsDown,secs,keylist] = KbCheck();
                if keyIsDown && (keylist(leftKey) == 1 ||  keylist(rightKey) == 1) && sum(keylist) == 1
                    % EYELINK: final resp seq x
                    if video.eyeON
                        msg = sprintf('ISQC_resp_%03d', isqc);
                        Eyelink('Message', msg);
                    end
                    respprobe = find(respKeys == (find(keylist==1)));
                    rtprobe  = secs - vbldisplaytime;
                    break
                end
            end
            
            % Redraw probe + show square around selected shape 0.8secs:
            switch condition
                case 2 ; line = sprintf('Dans cette sequence,\nde quelle couleur etait le symbole suivant : ');
                    Screen('DrawTexture', video.window, shape_tex(1,targetsymbol), [], centerpos, 0, [], [], [0,0,0]);
                    Screen('DrawTexture', video.window, shape_tex(2,targetsymbol), [], centerpos, 0, [], [], video.lbgd);
                case 1 ; line = sprintf('Dans cette sequence,\navez-vous eu l''impression de tirer plutot du : ');
            end
            DrawFormattedText(video.window,  line, 'center', video.screenYpixels*0.20, black);
            Screen('DrawTexture', video.window, rep_tex(1,redbuttonside), [], rep_pos(redbuttonside,:), 0, [], [], [0,0,0]);
            Screen('DrawTexture', video.window, rep_tex(2,redbuttonside), [], rep_pos(redbuttonside,:), 0, [], [], orange);
            Screen('DrawTexture', video.window, rep_tex(1,bluebuttonside), [], rep_pos(bluebuttonside,:), 0, [], [], [0,0,0]);
            Screen('DrawTexture', video.window, rep_tex(2,bluebuttonside), [], rep_pos(bluebuttonside,:), 0, [], [], blue);
            switch respprobe
                case goodanswer ; Screen('DrawTexture', video.window, fbsquare, [], fbpos(respprobe,:), [], [], [], black); %square
                case wronganswer ; Screen('DrawTexture', video.window, fbsquare, [], fbpos(respprobe,:), [], [], [], black); %square
            end
            Screen('DrawingFinished', video.window);
            Screen('Flip', video.window, secs+roundfp(0.1,0,video.ifi));
            WaitSecs(0.8);

            % EYELINK: fb seq x
            % EYELINK: end seq x
            if video.eyeON
                msg = sprintf('ISQC_fbsel_%03d', isqc);
                Eyelink('Message', msg); %WaitSecs(0.01);???
                msg = sprintf('ISQC_end_%03d', isqc);
                Eyelink('Message', msg);
            end
            
            Screen('FillRect', video.window, video.lbgd);
            Screen('Flip', video.window);
            WaitSecs(0.3);
            
            % Store final response: cuisine pour remettre info dans tous les samples: TODO deal()?
            for i = (itotal - (nsamples-1)) : itotal
                sampling(i).finaltargetside = goodanswer;
                sampling(i).finalresponse = respprobe;
                sampling(i).finalrt = rtprobe;
            end
            finalchoice(isqc).subj = subj;
            finalchoice(isqc).condition = condition;
            finalchoice(isqc).mirroring = sqc(isqc).mirroring;
            finalchoice(isqc).blocknb = sqc(isqc).blocknb;
            finalchoice(isqc).blkseqnb = sqc(isqc).blkseqnb;
            finalchoice(isqc).sequencenb = isqc;
            finalchoice(isqc).nsamples = sqc(isqc).nsamples;
            finalchoice(isqc).symbolpair = sqc(isqc).symbolpair;
            finalchoice(isqc).targetsymbol = targetsymbol;
            finalchoice(isqc).coloragency = sqc(isqc).coloragency;
            finalchoice(isqc).targetcolor = sqc(isqc).targetcolor;
            finalchoice(isqc).goodanswerside = sqc(isqc).probeside;
            finalchoice(isqc).finalresponse  = respprobe;
            finalchoice(isqc).finalrt        = rtprobe;
            
            % Store this sequence's points in blockscores structure
            switch condition
                case 1
                    switch targetcolor %colorslin = -1(blue) to +1(orange)
                        case 1; blockscores(iseq) = mean((seqpoints-1)/-2);%if targ blue rescale 1 to 0
                        case 2; blockscores(iseq) = mean((seqpoints+1)/2);  %if targ org rescale 0 to 1
                    end
                case 2
                    switch respprobe
                        case goodanswer; blockscores(iseq) = 1;
                        case wronganswer; blockscores(iseq) = 0;
                    end
            end
            
        end%FOR each sequence of samples in this block
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        

        
        % EYELINK: close EDF and transfer
        if video.eyeON
            %save
            WaitSecs(1);
            Eyelink('StopRecording');
            Eyelink('CloseFile');
            % flip to participant to wait:
            line = sprintf('... veuillez patienter pendant le transfert de vos donn�es pour ce bloc...');
            DrawFormattedText(video.window, line, 'center', video.screenYpixels * 0.5, black);
            Screen('Flip', video.window);
            %transfer
            edfpath = fullfile(pwd, '/', subjdir);
            gotit = Eyelink('ReceiveFile', edfname, edfpath, 1);
            fprintf('\nData transfered = %d Mo, file >>>%s.edf can be found in >>>%s. \n', gotit, edfname, edfpath);
            %convert to asc
            system(sprintf('edf2asc %s.edf', [edfpath '/' edfname]));
            fprintf('\nData converted from edf to asc file in local >>>%s. \n', edfpath);
        end
        
        % End of each block, insert pause & give score
        fprintf('Participant is on their 30'''' break!\nPress [000] and [y] to continue...\n')
        presSecs = 30:-1:0;
        line1 = sprintf('C''est la fin de ce bloc, bravo !\nVotre score pour ce bloc est de %0.f / 100!', mean(blockscores)*100);
        DrawFormattedText(video.window, line1, 'center', video.screenYpixels * 0.25, black);
        line2 = 'PRENEZ UNE PAUSE !\n\n';
        DrawFormattedText(video.window, line2, 'center', video.screenYpixels * 0.50, black);
        line3 = 'Quand vous etes pret-e a� continuer,\nallez chercher l''experimentateur.';
        DrawFormattedText(video.window, line3, 'center', video.screenYpixels * 0.75, black);
        Screen('DrawingFinished', video.window);
        Screen('Flip', video.window);
        k = RestrictKeysForKbCheck(spaceKey);
        WaitKeyPress(spaceKey);
        k = RestrictKeysForKbCheck(continueKey);
        WaitKeyPress(continueKey);
        
        
    end%FOR EACH BLOCK IN SQC STRUCTURE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % DRAW FINAL THANK YOU SCREEN
    line = 'Merci d''avoir participe, l''experience est finie !\n';
    endline = 'TERMINER -- [espace] -- TERMINER';
    DrawFormattedText(video.window, [line endline], 'center', video.yCenter, black);
    Screen('DrawingFinished', video.window);
    Screen('Flip', video.window);
    k = RestrictKeysForKbCheck(spaceKey);
    WaitKeyPress(spaceKey);
    sca;
    
    % write when ended and if crash:
    header.end = datestr(now,'yyyymmdd-HHMM');
    header.aborted = 0;
    save([subjdir '/' 'output' num2str(subj)], 'sampling', 'finalchoice', 'sqc', 'header', 'video')

    % close PTB & Eyelink
    closePTB();
    closeEL(eyeON, subj, iblock);    
    
    fprintf('\n\n\nWE''RE DONE\n\n\n');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% catching results so far if crashes
catch ME
           
    % clean save behaviour:
    header.end = datestr(now,'yyyymmdd-HHMM');
    header.aborted = 1;
    save([subjdir '/' 'outputsofar' num2str(subj)], 'sampling', 'finalchoice', 'sqc', 'header', 'video')
    
    closeEL(eyeON, subj, iblock);

    closePTB();
    rethrow(ME);
    
end% TRY/CATCH

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function closePTB() % clean close Psychtoolbox
Priority(0);
FlushEvents;
ListenChar(0);
ShowCursor;
Screen('CloseAll');
Screen('Preference', 'SkipSyncTests', 0);%restore synch testing
Screen('Preference','VisualDebuglevel',4);%restore warning verbosity
Screen('Preference','SuppressAllWarnings',0);%=1 == warning verbosity=0
k = RestrictKeysForKbCheck([]);
end

function closeEL(eyeON, subj, iblock) % save if recording and clean close Eyelink
if eyeON
    recording = Eyelink('Checkrecording');
    if recording == 0 % still recording rn
        %save
        WaitSecs(1);
        Eyelink('StopRecording');
        Eyelink('CloseFile');
        %transfer
        subjdir = sprintf('../Data/S%02d',subj);
        edfpath = fullfile(pwd, '/', subjdir);
        edfname = sprintf('COX%03d_%d', subj, iblock);
        gotit = Eyelink('ReceiveFile', edfname, edfpath, 1);
        fprintf('\nData transfered = %d Mo, file >>>%s.edf can be found in >>>%s. \n', gotit, edfname, edfpath);
        %convert to asc
        system(sprintf('edf2asc %s.edf', [edfpath '/' edfname]));
        fprintf('\nData converted from edf to asc file in local >>>%s. \n', edfpath);
    end
    Eyelink('Command','generate_default_targets = YES'); % YES:default or NO:custom
    Eyelink('Shutdown');
end
end

function [t] = roundfp(t,dt,videoifi)
% apply duration rounding policy for video flips
% where t  - desired (input)/rounded (output) duration (Cle:in secs??)
%       dt - desired uniform jitter on duration (default: none)
%       videoifi - Cle: interframeinterval
n = round(t/videoifi); % nb of frames in t duration
% apply uniform jitter
if nargin > 1 && dt > 0
    m = round(dt/videoifi);
    n = n+ceil((m*2+1)*rand)-(m+1);
end
% convert frames to duration
t = (n-0.5)*videoifi; % convert back to frames, but 0.5 frames before total duration
end


