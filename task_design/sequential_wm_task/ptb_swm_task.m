%% sequential working memory task

% stimulus presentation:
% in a circle with radius 1, and 0 degrees in location (1,0) [x=1, y=0],
% and runs counter clockwise, meaning 90 degrees in location (0,1)  [x=0, y=1]
% location 1 =   0 degrees
% location 2 =  45 degrees
% location 3 =  90 degrees
% location 4 = 135 degrees
% location 5 = 180 degrees
% location 6 = 225 degrees
% location 7 = 270 degrees
% location 8 = 315 degrees

% notes
% do not use PsychImaging(), because meg triggers do not work properly
% (only three channels out of 8)
% [window, windowRect] = PsychImaging('OpenWindow', screenNumber, black);
% things that should be changed when swithching from
% PsychImaging('OpenWindow',...) to Screen('Openwindow',...)
% 1) RGB values are 0-255 (instead of 0-1)

% response buttons Vpixx: only left hand and not lateral white button
% lateral (white) button in right hand is broken, sending signals as always on

% Clear the workspace and the screen
sca; close all; clearvars;

%% settings script

% settings experiment
debug = 1; % 
% 1 for coding in laptop (skip synchronization tests)
% 0 for experiment in MEG (keep synchronization tests)

little_window = 0; 
% display little window in the same screen: 1 true, 0 false

trigger_test = 0; 
% if 0, trigger is 1 pixel, 
% if 1 trigger is bigger (to be able to see it)

Vpixx_response = 2; 
% if 0, laptop keyboard responses (debugging) + trigger pixel mode off
% if 1, Vpixx response buttons (without dial) + trigger pixel mode on
    % response Vpixx if 1
    % left right: go left       
    % left green: go right
    % left white (lateral): select location
% if 2, Vpixx response buttons + dial         + trigger pixel mode on 

%% settings input

% task information
sequence_lengths = '34';
version_sequence_task = '1';

path_local = 'path_to_local';
path_meg = 'path_to_meg_computer';

if isfolder(path_local)
    path = path_local;
elseif isfolder(path_meg)
    path = path_meg;
end

% paths
path_input = [path '/scripts/sequential_wm_task/input'];
path_output = [path '/results/behavior'];

% subject information
subject = input(sprintf('\n\tSubject_id = '));
session = input(sprintf('\n\tsession_number = '));
number_trials = input(sprintf('\n\tnumber_of_trials = '));
block = input(sprintf('\n\tblock_number = '));

if session == 1
    randomization = 1;
elseif session == 2
    randomization = 2;
elseif session == 3
    randomization = 1;
elseif session == 4
    randomization = 2;
elseif session == 5
    randomization = 1;
elseif session == 6
    randomization = 2;
end

% create output folder
if subject < 10
    output_folder = [path_output '/sub_0' num2str(subject) '/sess_0' num2str(session) '/block_0' num2str(block)];
else
    output_folder = [path_output '/sub_' num2str(subject) '/sess_0' num2str(session) '/block_0' num2str(block)];
end

% if foldername exists, ask for resume session or start a new one
resume_session = 0;

if ~isfolder(output_folder)
    mkdir(output_folder);
else 
    disp(' ')
    disp('ERROR: existing folder with this subject ID and session');
    disp(' ')
    response_foldername = input(sprintf('Resume an existing session in the trial in which it was interrupted? yes (1), no (2): '));
    if response_foldername == 2
        disp(' ')
        disp('WARNING: Change subject ID, session and/or block to start a NEW session');
        disp(' ')
        return
    elseif response_foldername == 1
        disp(' ')
        disp('Resuming session at the trial in which it was interrupted')
        resume_session = 1;
        disp(' ')
    end
end

% new session
if resume_session == 0

    % set sequence trials
    load([path_input '/task_sequences_' num2str(number_trials) 'trials_'...
        sequence_lengths 'seq_v' version_sequence_task '_random' num2str(randomization) '_block' num2str(block) '.mat']) % sequences
    
    % save trial presentation after randomization
    save([output_folder '/task_trial_presentation.mat'], 'sequences');

    vector_trials = 1:size(sequences, 1);

    % empty matrix for responses
    responses = nan(size(sequences));

% resume session previously interrupted
elseif resume_session == 1

    % load sequences
    load([output_folder '/task_trial_presentation.mat']);

    % load responses
    load([output_folder '/responses_task_backup.mat']); % responses
    tmp = responses;
    tmp = tmp(:,1);
    tmp = tmp(~isnan(tmp));

    vector_trials = (length(tmp)+1):size(sequences, 1);

    % load dim fixation task
    load([output_folder '/task_trials_dimfixation_backup.mat']); % dim

end

%% settings dim flickering task (fixation point)

% settings 
dim_duration = 1; % second
luminance_change = [50 50 50]; % 10% of luminance change
frames_flickering = 2; % number of frames for flickering

if resume_session == 0
    
    % load trials with dim fixation point
    load([path_input '/dim_' num2str(number_trials) 'trials_'...
        sequence_lengths 'seq_v' version_sequence_task '_random' num2str(randomization) '_block' num2str(block) '.mat']);
    
    % save trials with dim fixation point
    save([output_folder '/task_trials_dimfixation.mat'], 'dim');

end

%% settings task

% % settings timing (seconds)
baseline_duration = 2.2; % always add 0.2 sec to baseline duration, to remove the first 200 msec, to unmix from response time
stim_duration = 0.3;
delay_duration = 4;
minimum_iti = 0.5;
isi_duration = 0.1;
feedback_duration = 1;

% settings stimulus feature
dotSizePix = 60; % size of stimulus
if little_window == 1
    excentricity = 200; % pixels
else 
    excentricity = 350; % pixels
end

% settings stimulus location
number_locations = 8; 
theta = linspace(0, 2*pi, number_locations+1);
theta = theta(1:(end-1));
% angle theta = 0 radians, in position (1,0), running counter clock wise
% [ pi / 2 radians in position (0,1) ]

% settings fixation point
size_fixation_point = 15; % size
color_fixation_point = [255 255 255];

% Pen width for the frames
penWidthPixels = 2;

% settings colours
grey  = [128 128 128];
black = [0 0 0];
white = [255 255 255];
red   = [255 0 0];
grey_dark = [80 80 80];

% Color is defined by red green and blue components (RGB). So we have three numbers which
% define our RGB values. The maximum number for each is 255 and the minimum
% 0. So, "full red" is [255 0 0]. "Full green" [0 255 0] and "full blue" [0 0
% 255]. 

% images feedback
tick_correct = [path_input '/tick_correct.png'];
tick_incorrect = [path_input '/tick_incorrect.png'];

% sequences training trials
sequences_train = [1 2 3 0 0;...
                   6 5 4 3 0;
                   7 6 5 6 7];

    
%% setting up Psychtoolbox

% skip synchronization test
if debug == 1
    Screen('Preference', 'SkipSyncTests', 1);
end

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer.
screens = Screen('Screens');

% To draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen.
screenNumber = max(screens);

% settings for small window (1) or full screen presentation
if little_window == 1

    % Start cordinate in pixels of our window. Note that setting both of these
    % to zero will make the window appear in the top right of the screen.
    startXpix = 1100;
    startYpix = 50;
    
    % Dimensions in pixels of our window in the X (left-right) and Y (up down) dimensions
    dimX = 800;
    dimY = 600;
    
    % Open (partial) window with black background  
    [window, windowRect] = Screen('OpenWindow', screenNumber, black,...
        [startXpix startYpix startXpix + dimX startYpix + dimY], [], [], [], [], [], kPsychGUIWindow);

else

    % Open (full) window with grey background  
    [window, windowRect] = Screen('OpenWindow', screenNumber, grey);

end

% Set up alpha-blending for smooth (anti-aliased) lines
% https://www.youtube.com/watch?v=7eFGY6JVTnc
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

% Set maximum priority level
topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);

% Get the size of the on screen window in pixels
[screenXpixels, screenYpixels] = Screen('WindowSize', window);

% Get the centre coordinate of the window in pixels
[xCenter, yCenter] = RectCenter(windowRect);

% The coordinates define the top left and bottom right coordinates of our rect
% [top-left-x top-left-y bottom-right-x bottom-right-y].
baseRect = [0 0 screenXpixels screenYpixels];

% Center the rectangle on the centre of the screen 
centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);

% Query the frame duration
ifi = Screen('GetFlipInterval', window);


%% settings of response buttons with keyboard

% unify key names across computers
KbName('UnifyKeyNames');

% The avaliable keys to press
escapeKey = KbName('ESCAPE');
upKey = KbName('UpArrow');
downKey = KbName('DownArrow');
leftKey = KbName('LeftArrow');
rightKey = KbName('RightArrow');
barKey = KbName('space');

%% settings dial

if Vpixx_response == 2

    LoadPsychHID 
    
    % matlab reads clockwise as 't' and conterclockwise as 'b'
    cwKey = KbName('t');
    ccwKey = KbName('b');
    
    % this value should be -1 for running on Windows in the Vpixx computer
    dialID = -1;
    
    % disable matlab listening keyboard input in command window 
    % (to avoid prompting a bunch of ts and bs from the dial)
    ListenChar(-1);
    
    % theta in degrees
    theta_deg = rad2deg(theta);
    
    % ts (or bs) per degrees of angle in the dial
    ts_or_bs_dial_angle = 5;

end

%% trigger and response setup (Vpixx)
% https://vpixx.com/vocal/pixelmode/

% triggers as color (RGB) of tope-left pixel in the screen
trigger_224 = [4  0 0]; % 224 meg channel
trigger_225 = [16 0 0]; % 225 meg channel
trigger_226 = [64 0 0]; % 226 meg channel
trigger_227 = [0  1 0]; % 227 meg channel
trigger_228 = [0  4 0]; % 228 meg channel
trigger_229 = [0 16 0]; % 229 meg channel
trigger_230 = [0 64 0]; % 230 meg channel
trigger_231 = [0 0  1]; % 231 meg channel

% trigger locations for stimulus presentation and response
trigger_location1 = trigger_224;
trigger_location2 = trigger_224 + trigger_225;
trigger_location3 = trigger_225;
trigger_location4 = trigger_226;
trigger_location5 = trigger_225 + trigger_226;
trigger_location6 = trigger_224 + trigger_225 + trigger_226;
trigger_location7 = trigger_227;
trigger_location8 = trigger_226 + trigger_227;

trigger_response_dim = trigger_228;
trigger_baseline_onset = trigger_229; % baseline onset
trigger_delay_onset = trigger_230; % delay onset
trigger_response_onset = trigger_231; % response onset, feedback onset

% size and location of trigger pixel
if trigger_test == 0
    baseRect_trigger = [0 0 1 1];
    centeredRect_trigger = CenterRectOnPointd(baseRect_trigger, 0.5, 0.5);
elseif trigger_test == 1
    baseRect_trigger = [0 0 50 50];
    centeredRect_trigger = CenterRectOnPointd(baseRect_trigger, 25, 25);
end

if Vpixx_response == 1 || Vpixx_response == 2

    Datapixx('Open');
    Datapixx('EnablePixelMode'); % trigger pixel mode
    Datapixx('RegWrRd'); % Synchronize DATAPixx registers to local register cache

end

%% prompt instructions

duration_text = 3;

if block == 1 && resume_session == 0

    % Draw grey background to the screen.
    Screen('FillRect', window, grey, centeredRect);
    
    % Trigger (black)
    Screen('FillRect', window, black, centeredRect_trigger);
    
    % Flip to the screen
    Screen('Flip', window);
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Now, the visual memory task', 'center', screenYpixels * 0.25, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
    
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Now, the visual memory task', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'You will see a red dot in different locations...', 'center', screenYpixels * 0.5, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
    
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Now, the visual memory task', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'You will see a red dot in different locations...', 'center', screenYpixels * 0.5, white);
    DrawFormattedText(window, '...and should reproduce the exact same sequence', 'center', screenYpixels * 0.75, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 45);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Always keep your gaze at the point in the center of the screen', 'center', screenYpixels * 0.33, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 45);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Always keep your gaze at the point in the center of the screen', 'center', screenYpixels * 0.33, white);
    DrawFormattedText(window, 'Press a button when the luminance of that point changes', 'center', screenYpixels * 0.66, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Please, do not move or blink', 'center', screenYpixels * 0.25, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Please, do not move or blink', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'Blink only while responding or between trials', 'center', screenYpixels * 0.5, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Please, do not move or blink', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'Blink only while responding or between trials', 'center', screenYpixels * 0.5, white);
    DrawFormattedText(window, 'We will start in a few seconds', 'center', screenYpixels * 0.75, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);

else

    % Draw grey background to the screen.
    Screen('FillRect', window, grey, centeredRect);
    
    % Trigger (black)
    Screen('FillRect', window, black, centeredRect_trigger);
    
    % Flip to the screen
    Screen('Flip', window);
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Now, we will resume the visual memory task', 'center', screenYpixels * 0.25, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
            
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Now, we will resume the visual memory task', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'Please, do not move or blink', 'center', screenYpixels * 0.5, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
            
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Now, we will resume the visual memory task', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'Please, do not move or blink', 'center', screenYpixels * 0.5, white);
    DrawFormattedText(window, 'Blink only while responding or between trials', 'center', screenYpixels * 0.75, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);

end

%% training trials

% this block is a copy of the 'loop over trials', but with 
% a) two sequences
% b) no triggers
% c) no response recording

% only at the beginning of the experiment
if block == 1 && resume_session == 0

    % set variable
    exit_task = 0;
    
    trials_train = size(sequences_train, 1);
    
    responses_train = nan(size(sequences_train));
    
    for train_i = 1:trials_train
       
        disp(' '); disp(' '); disp(' ');
        disp(['training trial ' num2str(train_i) ' out of ' num2str(trials_train)]);
        disp(' '); disp(' '); disp(' ');
    
        if exit_task == 1
    
            % enable matlab listen keyboard input in command window
            ListenChar
    
            return
    
        else
    
            % subset sequence (remove zeros)
            sequence = sequences_train(train_i,:);
            sequence = sequence(sequence ~= 0);
    
            %% ITI
            % only background, proceed with baseline after participant press a key
            % press button when ready
                
            % Draw grey background to the screen.
            Screen('FillRect', window, grey, centeredRect);
        
            % Trigger (black)
            Screen('FillRect', window, black, centeredRect_trigger);
        
            % Flip to the screen
            Screen('Flip', window);
        
            % minimum time for iti
            time = 0;
        
            while time < minimum_iti
        
                time = time + ifi;
        
            end
        
            % Draw text
            Screen('TextSize', window, 60);
            Screen('TextFont', window, 'Courier');
            DrawFormattedText(window, 'Blink now (if needed)', 'center', screenYpixels * 0.33, white);

            if Vpixx_response == 0
                DrawFormattedText(window, 'Press bar space when ready to start', 'center', screenYpixels * 0.66, white);
            elseif Vpixx_response == 1 || Vpixx_response == 2
                DrawFormattedText(window, 'Press any button when ready to start', 'center', screenYpixels * 0.66, white);
            end
            
            % Trigger (black)
            Screen('FillRect', window, black, centeredRect_trigger);
        
            % Flip to the screen
            Screen('Flip', window);
                            
            % set variable
            timeout = 0;
        
            if Vpixx_response == 1 || Vpixx_response == 2
        
                any_button = 0;
        
                % while timeout keeps being equal to zero
                while ~timeout
        
                    % Check the keyboard to see if a button has been pressed
                    % (only for escape key)
                    [keyIsDown,secs, keyCode] = KbCheck;
                
                    % Check the Vpixx response buttons to see if a button has been pressed
                    Datapixx('RegWrRd');
                    kbcheck = dec2bin(Datapixx('GetDinValues'));
                    
                    if kbcheck(end)   == '1' % left red
        
                        any_button = 1;
        
%                     elseif kbcheck(end-2) == '1' % left green
%         
%                         any_button = 1;
%         
%                     elseif kbcheck(end-3) == '1' % left blue
%         
%                         any_button = 1;
%         
%                     elseif kbcheck(end-1) == '1' % left yellow
%         
%                         any_button = 1;
                
                    end
        
                    % set variables if Vpixx response buttons and/or keyboard scape key are pressed
                    if any_button == 1 % Vpixx
                        timeout = 1;
                    elseif keyCode(escapeKey) % keyboard (only for espace key)
                        exit_task = 1;
                    end
                    
                end
        
                kbcheck = [];
        
            elseif Vpixx_response == 0
        
                % while timeout keeps being equal to zero
                while ~timeout
            
                    % Check the keyboard to see if a button has been pressed
                    [keyIsDown,secs, keyCode] = KbCheck;
            
                    % set variables if spacebar and/or scape keys are pressed
                    if keyCode(barKey)
                        timeout = 1;
                    elseif keyCode(escapeKey)
                        exit_task = 1;
                    end
                    
                end
        
            end
        
            %% baseline
                            
            % set variables
            timeout = 0; % logical variable to control time
            time = 0; % seconds
            trigger = 0;

            % Set color of background empty dots 
            ovalColorMatrix = repmat(white', 1, number_locations);
                            
            % Make a base Rect
            baseRect = [0 0 dotSizePix dotSizePix];
                
            % Center the rectangle on the centre of the screen to make the empty dots
            centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);
            
            % Vertex of empty dots (1-3 x axis; 2-4 y-axis) 
            ovalPositionMatrix = nan(4, number_locations);
            for pos_i = 1:number_locations
            
                ovalPositionMatrix(1, pos_i) = centeredRect(1) + excentricity * cos(theta(pos_i));
                ovalPositionMatrix(3, pos_i) = centeredRect(3) + excentricity * cos(theta(pos_i));
            
                ovalPositionMatrix(2, pos_i) = centeredRect(2) - excentricity * sin(theta(pos_i));
                ovalPositionMatrix(4, pos_i) = centeredRect(4) - excentricity * sin(theta(pos_i));
            
            end
            
            % while timeout keeps being equal to zero
            while ~timeout
        
                % exit task if scape bar is pressed
                if exit_task == 0
                    [keyIsDown,secs, keyCode] = KbCheck;
                    if keyCode(escapeKey)
                        exit_task = 1;
                    end
                end
            
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);

                % Draw empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
        
                % draw pixel trigger
                if ~trigger
        
                    % Trigger (black) - no triggers during training trials
                    Screen('FillRect', window, black, centeredRect_trigger);
    
                    trigger = 1;
        
                else
        
                    % Trigger (black)
                    Screen('FillRect', window, black, centeredRect_trigger);
        
                end
            
                % Flip to the screen
                Screen('Flip', window);
                            
                % Increment the time
                time = time + ifi;
            
                % update timeout (when baseline duration is completed)
                if time >= baseline_duration
                    timeout = 1;
                end
            
            end
            
            %% stimulus presentation
                    
            % exit task if scape bar is pressed
            if exit_task == 0
                [keyIsDown,secs, keyCode] = KbCheck;
                if keyCode(escapeKey)
                    exit_task = 1;
                end
            end
        
            % Set color of background empty dots 
            ovalColorMatrix = repmat(white', 1, number_locations);
        
            % Set dot size in pixels
            dotSizePixMatrix = repmat(dotSizePix, 1, number_locations);
                    
            % Make a base Rect
            baseRect = [0 0 dotSizePix dotSizePix];
                
            % Center the rectangle on the centre of the screen to make the empty dots
            centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);
            
            % Vertex of empty dots (1-3 x axis; 2-4 y-axis) 
            ovalPositionMatrix = nan(4, number_locations);
            for pos_i = 1:number_locations
            
                ovalPositionMatrix(1, pos_i) = centeredRect(1) + excentricity * cos(theta(pos_i));
                ovalPositionMatrix(3, pos_i) = centeredRect(3) + excentricity * cos(theta(pos_i));
            
                ovalPositionMatrix(2, pos_i) = centeredRect(2) - excentricity * sin(theta(pos_i));
                ovalPositionMatrix(4, pos_i) = centeredRect(4) - excentricity * sin(theta(pos_i));
            
            end
        
            % coordinates of coloured dots
            dotPositionMatrix = nan(2, number_locations); % rows (x-position, y-position), columns (locations)
            for pos_i = 1:number_locations
                dotPositionMatrix(1, pos_i) = xCenter + excentricity * cos(theta(pos_i));
                dotPositionMatrix(2, pos_i) = yCenter - excentricity * sin(theta(pos_i));
            end
    
            % set variables
            trigger = 0;
            
            % display sequence of dots
            for stim_i = 1:length(sequence)
        
                % set variables
                timeout_stimulus = 0;
                time = 0; % seconds
        
                % while timeout_stimulus keeps being zero
                while ~timeout_stimulus
            
                    % set color, size and position of stimuli
                    dotSizePixMatrix_stim = dotSizePixMatrix(sequence(stim_i));
                    dotPositionMatrix_stim = dotPositionMatrix(:, sequence(stim_i));
                
                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);
            
                    % Draw empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
        
                    % Draw stimulus dot to the screen. 
                    Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, red, [], 2);
        
                    % draw pixel trigger
                    if ~trigger
            
                        % Trigger (colour) - no trigger during training trials
                        if sequence(stim_i) == 1
                            trigger_stim = trigger_location1;
                        elseif sequence(stim_i) == 2
                            trigger_stim = trigger_location2;
                        elseif sequence(stim_i) == 3
                            trigger_stim = trigger_location3;
                        elseif sequence(stim_i) == 4
                            trigger_stim = trigger_location4;
                        elseif sequence(stim_i) == 5
                            trigger_stim = trigger_location5;
                        elseif sequence(stim_i) == 6
                            trigger_stim = trigger_location6;
                        elseif sequence(stim_i) == 7
                            trigger_stim = trigger_location7;
                        elseif sequence(stim_i) == 8
                            trigger_stim = trigger_location8;
                        end
                            
                        Screen('FillRect', window, black, centeredRect_trigger);
    
                        trigger = 1;
        
                    else
        
                        % Trigger (black)
                        Screen('FillRect', window, black, centeredRect_trigger);
        
                    end
        
                    % Flip to the screen. 
                    Screen('Flip', window);
                
                    % Increment the time
                    time = time + ifi;
        
                    % update timeout_stimulus variable
                    if time > stim_duration
                        timeout_stimulus = 1;
                    end
        
                end
        
                timeout_isi = 0;
        
                while ~timeout_isi
        
                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);
                
                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
            
                    % trigger
                    Screen('FillRect', window, black, centeredRect_trigger);
                    
                    % Flip to the screen. 
                    Screen('Flip', window);
            
                    % Increment the time
                    time = time + ifi;
            
                    % update timeout_isi variable
                    if time > isi_duration + stim_duration
            
                        timeout_isi = 1;
                        
                    end
        
        
                end
        
            end
            
            %% delay
            
            % set variables
            timeout = 0; % logical variable to control time
            time = 0; % seconds
    
            trigger = 0;
    
            dim_button = 0;
    
            flickering_count = 0;
            
            % while timeout keeps being zero
            while ~timeout
        
                % exit task if scape bar is pressed
                if exit_task == 0
                    [keyIsDown,secs, keyCode] = KbCheck;
                    if keyCode(escapeKey)
                        exit_task = 1;
                    end
                end
    
        
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);

                % Draw empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
        
                % draw pixel trigger 
                if ~trigger
        
                    % Trigger (colour) - no trigger during training trials
                    Screen('FillRect', window, black, centeredRect_trigger);
    
                    trigger = 1;
        
                else
        
                    % Trigger (black)
                    Screen('FillRect', window, black, centeredRect_trigger);
        
                end
            
                % Flip to the screen
                Screen('Flip', window);
            
                % Increment the time
                time = time + ifi;
            
                % update timeout
                if time >= delay_duration
                    timeout = 1;
                end
            
            end
            
            %% response
        
            if Vpixx_response == 2
    
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);
                
                % empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                        
                % Trigger (black)
                Screen('FillRect', window, black, centeredRect_trigger);
                
                % Flip to the screen
                Screen('Flip', window);
    
                % set variables
                empty_kbcheck   =  '111111111111110000000000';
                last_kbcheck    =  '000111111111110000000000';
                left_red_button =  '111111111111110000000001';
    
%                 ready2start = 0;
                
%                 % start response as soon as stop pressing any button
%                 while ~ready2start
%     
%                     % check the Vpixx response buttons
%                     Datapixx('RegWrRd');
%                     kbcheck = dec2bin(Datapixx('GetDinValues'));
%     
%                     if strcmp(kbcheck, empty_kbcheck) == 1
%                         ready2start = 1;
%                     end
%     
%                 end
                                                    
                % Set the color of filled dot of selected spatial location (grey)
                dotColorMatrix = repmat(grey', 1, number_locations);
                        
                % start at random location (not first location in sequence)
                random_position = 0;
                while ~random_position
                    location_response = randi([1 length(sequence)], 1);
                    if location_response ~= sequence(1)
                        random_position = 1;
                    end
                end
    
                % location to degrees
                if location_response == 1
                    dial_degrees = 0;
                elseif location_response == 2
                    dial_degrees = 45;
                elseif location_response == 3
                    dial_degrees = 90;
                elseif location_response == 4
                    dial_degrees = 135;
                elseif location_response == 5
                    dial_degrees = 180;
                elseif location_response == 6
                    dial_degrees = 225;
                elseif location_response == 7
                    dial_degrees = 270;
                elseif location_response == 8
                    dial_degrees = 315;
                end
    
                % set first location
                dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                                    
                % set variables
                number_responses = 0;
    
                % while # responses is below length of sequence
                while (number_responses < length(sequence))
                
                    % set variables
                    select_location = 0;
                    update_location = 0;
        
                    % color response position
                    dotColorMatrix_stim = grey_dark;
                    
                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);
                    
                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                    
                    % Draw the dot to the screen.
                    Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
                
                    % Trigger (black)
                    Screen('FillRect', window, black, centeredRect_trigger);
                    
                    % Flip to the screen
                    Screen('Flip', window);
                
                    % Check the dial
                    [keyIsDown, secs, keyCode] = KbCheck(dialID);
    
                    if keyCode(cwKey)
    
                        % update dial location (degrees)
                        dial_degrees = dial_degrees - ts_or_bs_dial_angle;
    
                        update_location = 1;
    
                    elseif keyCode(ccwKey)
    
                        % update dial location (degrees)
                        dial_degrees = dial_degrees + ts_or_bs_dial_angle;
    
                        update_location = 1;
    
                    end
    
                    % keep in range [0 360)
                    if ~(dial_degrees < 360)
                        dial_degrees = dial_degrees - 360;
                    elseif dial_degrees < 0
                        dial_degrees = dial_degrees + 360;
                    end
    
                    % check the Vpixx response buttons
                    Datapixx('RegWrRd');
                    kbcheck = dec2bin(Datapixx('GetDinValues'));
    
                    % reset all responses as left red
                    if kbcheck(end)   == '1'
        
                        kbcheck = left_red_button;
        
%                     elseif kbcheck(end-2) == '1' % left green
%         
%                         kbcheck = left_red_button;
%         
%                     elseif kbcheck(end-3) == '1' % left blue
%         
%                         kbcheck = left_red_button;
%         
%                     elseif kbcheck(end-1) == '1' % left yellow
%         
%                         kbcheck = left_red_button;
                
                    end
    
                    % if respone is different than previous response 
                    if strcmp(kbcheck, last_kbcheck) ~= 1
    
                        last_kbcheck = kbcheck;
            
                        if kbcheck(end)   == '1' % left red - % SELECT LOCATION
            
                            select_location = 1;
            
%                         elseif kbcheck(end-2) == '1' % left green
%             
%                             select_location = 1;
%             
%                         elseif kbcheck(end-3) == '1' % left blue
%             
%                             select_location = 1;
%             
%                         elseif kbcheck(end-1) == '1' % left yellow
%             
%                             select_location = 1;
                    
                        end
    
                    end
    
                        
                    % update location
                    if update_location == 1
    
                        % update location
                        if dial_degrees < (0 + 22.5) || dial_degrees > (360 - 22.5)
                            location_response = 1;
    
                        elseif dial_degrees < (45 + 22.5) && dial_degrees > (45 - 22.5)
                            location_response = 2;
    
                        elseif dial_degrees < (90 + 22.5) && dial_degrees > (90 - 22.5)
                            location_response = 3;
    
                        elseif dial_degrees < (135 + 22.5) && dial_degrees > (135 - 22.5)
                            location_response = 4;
    
                        elseif dial_degrees < (180 + 22.5) && dial_degrees > (180 - 22.5)
                            location_response = 5;
    
                        elseif dial_degrees < (225 + 22.5) && dial_degrees > (225 - 22.5)
                            location_response = 6;
    
                        elseif dial_degrees < (270 + 22.5) && dial_degrees > (270 - 22.5)
                            location_response = 7;
    
                        elseif dial_degrees < (315 + 22.5) && dial_degrees > (315 - 22.5)
                            location_response = 8;
    
                        end
    
                        % update location of dark_grey dot
                        dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                        dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
    
                        % Draw fixation point
                        Screen('DrawDots', window, [xCenter yCenter],...
                            size_fixation_point, color_fixation_point, [], 2);
                    
                        % empty backgrond dots
                        Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                            
                        % Draw the dot to the screen.
                        Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
        
                        % Trigger (black)
                        Screen('FillRect', window, black, centeredRect_trigger);
    
                        % Flip to the screen
                        Screen('Flip', window);
    
                    % select location
                    elseif select_location == 1
        
                        % change color of selected location to red
                        dotColorMatrix_stim = red;
            
                        % save response
                        number_responses = number_responses + 1;
                        responses_train(train_i, number_responses) = location_response;
    
                        % Draw fixation point
                        Screen('DrawDots', window, [xCenter yCenter],...
                            size_fixation_point, color_fixation_point, [], 2);
                    
                        % empty backgrond dots
                        Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                            
                        % Draw the dot to the screen.
                        Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
        
                        % Trigger (color) - no triggers during training trials
                        if location_response == 1
                            trigger_resp = trigger_224;
                        elseif location_response == 2
                            trigger_resp = trigger_225;
                        elseif location_response == 3
                            trigger_resp = trigger_226;
                        elseif location_response == 4
                            trigger_resp = trigger_227;
                        elseif location_response == 5
                            trigger_resp = trigger_228;
                        elseif location_response == 6
                            trigger_resp = trigger_229;
                        elseif location_response == 7
                            trigger_resp = trigger_230;
                        elseif location_response == 8
                            trigger_resp = trigger_231;
                        end                        
    
                        Screen('FillRect', window, black, centeredRect_trigger);
    
                        % Flip to the screen
                        Screen('Flip', window);
    
                        % Draw fixation point
                        Screen('DrawDots', window, [xCenter yCenter],...
                            size_fixation_point, color_fixation_point, [], 2);
                    
                        % empty backgrond dots
                        Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                            
                        % Draw the dot to the screen.
                        Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
        
                        % Trigger (black)
                        Screen('FillRect', window, black, centeredRect_trigger);
    
                        % Flip to the screen
                        Screen('Flip', window);
    
                        WaitSecs(0.2);
    
                    end
    
                end
                                
                % change color of selected location to red
                dotColorMatrix_stim = grey;
    
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);
            
                % empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                    
                % Trigger
                Screen('FillRect', window, black, centeredRect_trigger);
    
                % Flip to the screen
                Screen('Flip', window);
            
                % pause after all responses were provided
                WaitSecs(0.5);
    
            elseif Vpixx_response == 1
                                                    
                % Set the color of filled dot of selected spatial location (grey)
                dotColorMatrix = repmat(grey', 1, number_locations);
                        
                % start at random location (not first location in sequence)
                random_position = 0;
                while ~random_position
                    location_response = randi([1 length(sequence)], 1);
                    if location_response ~= sequence(1)
                        random_position = 1;
                    end
                end
                
                % set first location
                dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                                    
                % set variables
                number_responses = 0;
    
                % empty response
                empty_kbcheck = [];
                last_kbcheck = [];
                 
                % while # responses is below length of sequence
                while (number_responses < length(sequence))
                
                    % set variables
                    updated_response = 0;
                    wait_location = 0;
        
                    % color response position
                    dotColorMatrix_stim = grey_dark;
                    
                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);
                    
                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                    
                    % Draw the dot to the screen.
                    Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
                
                    % Trigger (black)
                    Screen('FillRect', window, black, centeredRect_trigger);
                    
                    % Flip to the screen
                    Screen('Flip', window);
                
                    % Check the keyboard to see if a button has been pressed
                    [keyIsDown, secs, keyCode] = KbCheck;
        
                    % check the buttonbox response
                    Datapixx('RegWrRd');
                    kbcheck = dec2bin(Datapixx('GetDinValues'));
    
                    % empty/last response
                    if isempty(empty_kbcheck)
                        empty_kbcheck = zeros(1, length(kbcheck));
                    end
    
                    if isempty(last_kbcheck)
                        last_kbcheck = ones(1, length(kbcheck));
                    end
                    
                    % exit task if scape bar has been pressed
                    if keyCode(escapeKey)
                        
                        if exit_task == 0
                            exit_task = 1;
                        end
        
                    end
    
                    % if kbcheck is different than last_kbcheck OR kbcheck is empty
                    if sum(kbcheck ~= last_kbcheck) > 0 || sum(kbcheck ~= empty_kbcheck) == 0
    
                        % last button pressed
                        last_kbcheck = kbcheck;
    
                        % response button Vpixx
                        if kbcheck(end)   == '1' % left red - GO LEFT
        
                            if updated_response == 0
            
                                location_response = location_response + 1;
                                if location_response == 9
                                    location_response = 1;
                                end
                                dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                                dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                    
                                updated_response = 1;
        
                            end
            
                        elseif kbcheck(end-2) == '1' % left green - GO RIGHT
            
                            if updated_response == 0
        
                                location_response = location_response - 1;
                                if location_response == 0
                                    location_response = 8;
                                end
                                dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                                dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                    
                                updated_response = 1;
        
                            end
            
                        elseif kbcheck(end-4) == '1' % left white - SELECT LOCATION
        
                            if wait_location == 0
            
                                dotColorMatrix_stim = red;
                    
                                % save response
                                number_responses = number_responses + 1;
                                responses_train(train_i, number_responses) = location_response;
        
                                wait_location = 1;
        
                            end
            
                        end
                    
                        % Draw fixation point
                        Screen('DrawDots', window, [xCenter yCenter],...
                            size_fixation_point, color_fixation_point, [], 2);
                    
                        % empty backgrond dots
                        Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                            
                        % Draw the dot to the screen.
                        Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
                            
                        if wait_location == 1 
            
                            % Trigger (color) - no triggers during training trials
    
                            if location_response == 1
                                trigger_resp = trigger_224;
                            elseif location_response == 2
                                trigger_resp = trigger_225;
                            elseif location_response == 3
                                trigger_resp = trigger_226;
                            elseif location_response == 4
                                trigger_resp = trigger_227;
                            elseif location_response == 5
                                trigger_resp = trigger_228;
                            elseif location_response == 6
                                trigger_resp = trigger_229;
                            elseif location_response == 7
                                trigger_resp = trigger_230;
                            elseif location_response == 8
                                trigger_resp = trigger_231;
                            end                        
    
                            Screen('FillRect', window, black, centeredRect_trigger);
            
                        else
            
                            % Trigger (black)
                            Screen('FillRect', window, black, centeredRect_trigger);
            
                        end
                        
                        % Flip to the screen
                        Screen('Flip', window);
                    
                        % update variables                
                        if wait_location == 1                
                            WaitSecs(0.2);
                        end
                                    
                    end
    
                end
            
                % pause after all responses were provided
                WaitSecs(0.5);
        
            elseif Vpixx_response == 0
        
                % Set the color of filled dot of selected spatial location (grey)
                dotColorMatrix = repmat(grey', 1, number_locations);
                        
                % start at random location (not first location in sequence)
                random_position = 0;
                while ~random_position
                    location_response = randi([1 length(sequence)], 1);
                    if location_response ~= sequence(1)
                        random_position = 1;
                    end
                end
                
                % set first location
                dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                                    
                % set variables
                number_responses = 0;
    
                % empty response
                empty_keyCode = [];
                last_keyCode = [];
                 
                % while # responses is below length of sequence
                while (number_responses < length(sequence))
                
                    % set variables
                    updated_response = 0;
                    wait_location = 0;
        
                    % color response position
                    dotColorMatrix_stim = grey_dark;
                    
                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);
                    
                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                    
                    % Draw the dot to the screen.
                    Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
                
                    % Trigger (black)
                    Screen('FillRect', window, black, centeredRect_trigger);
                    
                    % Flip to the screen
                    Screen('Flip', window);
                
                    % Check the keyboard to see if a button has been pressed
                    [keyIsDown,secs, keyCode] = KbCheck;
    
                    % empty/last response
                    if isempty(empty_keyCode)
                        empty_keyCode = zeros(1, length(keyCode));
                    end
    
                    if isempty(last_keyCode)
                        last_keyCode = ones(1, length(keyCode));
                    end
    
                    % if keyCode is different than last_keyCode OR keyCode is empty
                    if sum(keyCode ~= last_keyCode) > 0 || sum(keyCode ~= empty_keyCode) == 0
    
                        % last button pressed
                        last_keyCode = keyCode;
                        
                        % exit task if scape bar has been pressed
                        if keyCode(escapeKey)
                            
                            exit_task = 1;
                            
                        % update location based on response 
                        elseif keyCode(upKey)
                            
                            if updated_response == 0
            
                                location_response = location_response + 1;
                                if location_response == 9
                                    location_response = 1;
                                end
                                dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                                dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                    
                                updated_response = 1;
        
                            end
                            
                        % update location based on response
                        elseif keyCode(downKey)
                            
                            if updated_response == 0
        
                                location_response = location_response - 1;
                                if location_response == 0
                                    location_response = 8;
                                end
                                dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                                dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                    
                                updated_response = 1;
        
                            end
                    
                        % highlight selected location in red
                        elseif keyCode(barKey)
                    
                            if wait_location == 0
            
                                dotColorMatrix_stim = red;
                    
                                % save response
                                number_responses = number_responses + 1;
                                responses_train(train_i, number_responses) = location_response;
        
                                wait_location = 1;
        
                            end
                            
                        end
                    
                        % Draw fixation point
                        Screen('DrawDots', window, [xCenter yCenter],...
                            size_fixation_point, white, [], 2);
                    
                        % empty backgrond dots
                        Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                            
                        % Draw the dot to the screen.
                        Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
                            
                        if wait_location == 1 
            
                            % Trigger (color) - no triggers during training trials
    
                            if location_response == 1
                                trigger_resp = trigger_224;
                            elseif location_response == 2
                                trigger_resp = trigger_225;
                            elseif location_response == 3
                                trigger_resp = trigger_226;
                            elseif location_response == 4
                                trigger_resp = trigger_227;
                            elseif location_response == 5
                                trigger_resp = trigger_228;
                            elseif location_response == 6
                                trigger_resp = trigger_229;
                            elseif location_response == 7
                                trigger_resp = trigger_230;
                            elseif location_response == 8
                                trigger_resp = trigger_231;
                            end                        
    
                            Screen('FillRect', window, black, centeredRect_trigger);
            
                        else
            
                            % Trigger (black)
                            Screen('FillRect', window, black, centeredRect_trigger);
            
                        end
                        
                        % Flip to the screen
                        Screen('Flip', window);
                    
                        % update variables                
                        if wait_location == 1                
                            WaitSecs(0.2);
                        end
    
                    end
                                
                end
                                    
                % pause after all responses were provided
                pause(0.4);
                            
            end
        
            %% feedback
        
            % if correct/incorrect responses, load image of correct/incorrect tick
            responses_trial = responses_train(train_i, 1:length(sequence));
        
            if sum(responses_trial == sequence) == length(sequence)
                image = imread(tick_correct);
            else
                image = imread(tick_incorrect);
            end
            
            % Make the image into a texture
            imageTexture = Screen('MakeTexture', window, image);
        
            % Trigger (colour) - no triggers during training trials
            Screen('FillRect', window, black, centeredRect_trigger);
        
            % Flip screen
            Screen('Flip', window);
        
            time = 0;
        
            while time < feedback_duration
        
                time = time + ifi;
                
                % Draw grey background to the screen.
                Screen('FillRect', window, grey, centeredRect);
        
                % Draw the image to the screen
                Screen('DrawTexture', window, imageTexture, [], [], 0);
                    
                % Flip screen
                Screen('Flip', window);
        
            end
    
        end
            
    end

end

%% loop over trials

% set variable
exit_task = 0;

correct_responses = 0;

for seq_i = vector_trials
   
    disp(' '); disp(' '); disp(' ');
    disp(['trial ' num2str(seq_i) ' out of ' num2str(size(sequences, 1))]);
    if seq_i > 1
        disp(['Correct responses: ' num2str(correct_responses)]);
    end
    disp(' '); disp(' '); disp(' ');

    if exit_task == 1

        % enable matlab listen keyboard input in command window
        ListenChar

        return

    else

        % subset sequence (remove zeros)
        sequence = sequences(seq_i,:);
        sequence = sequence(sequence ~= 0);

        %% ITI
        % only background, proceed with baseline after participant press a key
        % press button when ready
            
        % Draw grey background to the screen.
        Screen('FillRect', window, grey, centeredRect);
    
        % Trigger (black)
        Screen('FillRect', window, black, centeredRect_trigger);
    
        % Flip to the screen
        Screen('Flip', window);
    
        % minimum time for iti
        time = 0;
    
        while time < minimum_iti
    
            time = time + ifi;
    
        end
    
        % Draw text
        Screen('TextSize', window, 60);
        Screen('TextFont', window, 'Courier');
        DrawFormattedText(window, 'Blink now (if needed)', 'center', screenYpixels * 0.33, white);

        if Vpixx_response == 0
            DrawFormattedText(window, 'Press bar space when ready to start', 'center', screenYpixels * 0.66, white);
        elseif Vpixx_response == 1 || Vpixx_response == 2
            DrawFormattedText(window, 'Press any button when ready to start', 'center', screenYpixels * 0.66, white);
        end
        
        % Trigger (black)
        Screen('FillRect', window, black, centeredRect_trigger);
    
        % Flip to the screen
        Screen('Flip', window);
    
        % set variable
        timeout = 0;
    
        if Vpixx_response == 1 || Vpixx_response == 2
    
            any_button = 0;
    
            % while timeout keeps being equal to zero
            while ~timeout
    
                % Check the keyboard to see if a button has been pressed
                % (only for escape key)
                [keyIsDown,secs, keyCode] = KbCheck;
            
                % Check the Vpixx response buttons to see if a button has been pressed
                Datapixx('RegWrRd');
                kbcheck = dec2bin(Datapixx('GetDinValues'));
                
                if kbcheck(end)   == '1' % left red
    
                    any_button = 1;
    
%                 elseif kbcheck(end-2) == '1' % left green
%     
%                     any_button = 1;
%     
%                 elseif kbcheck(end-3) == '1' % left blue
%     
%                     any_button = 1;
%     
%                 elseif kbcheck(end-1) == '1' % left yellow
%     
%                     any_button = 1;
        
                end
    
                % set variables if Vpixx response buttons and/or keyboard scape key are pressed
                if any_button == 1 % Vpixx
                    timeout = 1;
                elseif keyCode(escapeKey) % keyboard (only for espace key)
                    exit_task = 1;
                end
                
            end
    
            kbcheck = [];
    
        elseif Vpixx_response == 0
    
            % while timeout keeps being equal to zero
            while ~timeout
        
                % Check the keyboard to see if a button has been pressed
                [keyIsDown,secs, keyCode] = KbCheck;
        
                % set variables if spacebar and/or scape keys are pressed
                if keyCode(barKey)
                    timeout = 1;
                elseif keyCode(escapeKey)
                    exit_task = 1;
                end
                
            end
    
        end
    
        %% baseline
                        
        % set variables
        timeout = 0; % logical variable to control time
        time = 0; % seconds
        trigger = 0;

        % Set color of background empty dots 
        ovalColorMatrix = repmat(white', 1, number_locations);
                        
        % Make a base Rect
        baseRect = [0 0 dotSizePix dotSizePix];
            
        % Center the rectangle on the centre of the screen to make the empty dots
        centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);
        
        % Vertex of empty dots (1-3 x axis; 2-4 y-axis) 
        ovalPositionMatrix = nan(4, number_locations);
        for pos_i = 1:number_locations
        
            ovalPositionMatrix(1, pos_i) = centeredRect(1) + excentricity * cos(theta(pos_i));
            ovalPositionMatrix(3, pos_i) = centeredRect(3) + excentricity * cos(theta(pos_i));
        
            ovalPositionMatrix(2, pos_i) = centeredRect(2) - excentricity * sin(theta(pos_i));
            ovalPositionMatrix(4, pos_i) = centeredRect(4) - excentricity * sin(theta(pos_i));
        
        end
        
        % while timeout keeps being equal to zero
        while ~timeout
    
            % exit task if scape bar is pressed
            if exit_task == 0
                [keyIsDown,secs, keyCode] = KbCheck;
                if keyCode(escapeKey)
                    exit_task = 1;
                end
            end
        
            % Draw fixation point
            Screen('DrawDots', window, [xCenter yCenter],...
                size_fixation_point, color_fixation_point, [], 2);

            % Draw empty backgrond dots
            Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
    
            % draw pixel trigger
            if ~trigger
    
                % Trigger (colour)
                Screen('FillRect', window, trigger_baseline_onset, centeredRect_trigger);

                trigger = 1;
    
            else
    
                % Trigger (black)
                Screen('FillRect', window, black, centeredRect_trigger);
    
            end
        
            % Flip to the screen
            Screen('Flip', window);
                        
            % Increment the time
            time = time + ifi;
        
            % update timeout (when baseline duration is completed)
            if time >= baseline_duration
                timeout = 1;
            end
        
        end
        
        %% stimulus presentation
                
        % exit task if scape bar is pressed
        if exit_task == 0
            [keyIsDown,secs, keyCode] = KbCheck;
            if keyCode(escapeKey)
                exit_task = 1;
            end
        end
    
        % Set color of background empty dots 
        ovalColorMatrix = repmat(white', 1, number_locations);
    
        % Set dot size in pixels
        dotSizePixMatrix = repmat(dotSizePix, 1, number_locations);
                
        % Make a base Rect
        baseRect = [0 0 dotSizePix dotSizePix];
            
        % Center the rectangle on the centre of the screen to make the empty dots
        centeredRect = CenterRectOnPointd(baseRect, xCenter, yCenter);
        
        % Vertex of empty dots (1-3 x axis; 2-4 y-axis) 
        ovalPositionMatrix = nan(4, number_locations);
        for pos_i = 1:number_locations
        
            ovalPositionMatrix(1, pos_i) = centeredRect(1) + excentricity * cos(theta(pos_i));
            ovalPositionMatrix(3, pos_i) = centeredRect(3) + excentricity * cos(theta(pos_i));
        
            ovalPositionMatrix(2, pos_i) = centeredRect(2) - excentricity * sin(theta(pos_i));
            ovalPositionMatrix(4, pos_i) = centeredRect(4) - excentricity * sin(theta(pos_i));
        
        end
    
        % coordinates of coloured dots
        dotPositionMatrix = nan(2, number_locations); % rows (x-position, y-position), columns (locations)
        for pos_i = 1:number_locations
            dotPositionMatrix(1, pos_i) = xCenter + excentricity * cos(theta(pos_i));
            dotPositionMatrix(2, pos_i) = yCenter - excentricity * sin(theta(pos_i));
        end
        
        % display sequence of dots
        for stim_i = 1:length(sequence)

            % set variables
            trigger = 0;
    
            % set variables
            timeout_stimulus = 0;
            time = 0; % seconds
    
            % while timeout_stimulus keeps being zero
            while ~timeout_stimulus
        
                % set color, size and position of stimuli
                dotSizePixMatrix_stim = dotSizePixMatrix(sequence(stim_i));
                dotPositionMatrix_stim = dotPositionMatrix(:, sequence(stim_i));
            
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);
        
                % Draw empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
    
                % Draw stimulus dot to the screen. 
                Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, red, [], 2);
    
                % draw pixel trigger
                if ~trigger
        
                    % Trigger (colour)
                    if sequence(stim_i) == 1
                        trigger_stim = trigger_location1;
                    elseif sequence(stim_i) == 2
                        trigger_stim = trigger_location2;
                    elseif sequence(stim_i) == 3
                        trigger_stim = trigger_location3;
                    elseif sequence(stim_i) == 4
                        trigger_stim = trigger_location4;
                    elseif sequence(stim_i) == 5
                        trigger_stim = trigger_location5;
                    elseif sequence(stim_i) == 6
                        trigger_stim = trigger_location6;
                    elseif sequence(stim_i) == 7
                        trigger_stim = trigger_location7;
                    elseif sequence(stim_i) == 8
                        trigger_stim = trigger_location8;
                    end
                        
                    Screen('FillRect', window, trigger_stim, centeredRect_trigger);

                    trigger = 1;
    
                else
    
                    % Trigger (black)
                    Screen('FillRect', window, black, centeredRect_trigger);
    
                end
    
                % Flip to the screen. 
                Screen('Flip', window);
            
                % Increment the time
                time = time + ifi;
    
                % update timeout_stimulus variable
                if time > stim_duration
                    timeout_stimulus = 1;
                end
    
            end
    
            timeout_isi = 0;
    
            while ~timeout_isi
    
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);
            
                % empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
        
                % trigger
                Screen('FillRect', window, black, centeredRect_trigger);
                
                % Flip to the screen. 
                Screen('Flip', window);
        
                % Increment the time
                time = time + ifi;
        
                % update timeout_isi variable
                if time > isi_duration + stim_duration
        
                    timeout_isi = 1;
                    
                end
    
    
            end
    
        end
        
        %% delay
        
        % set variables
        timeout = 0; % logical variable to control time
        time = 0; % seconds

        trigger = 0;

        dim_button = 0;

        dim_button_done = 0;

        flickering_count = 0;
        
        % while timeout keeps being zero
        while ~timeout
    
            % exit task if scape bar is pressed
            if exit_task == 0
                [keyIsDown,secs, keyCode] = KbCheck;
                if keyCode(escapeKey)
                    exit_task = 1;
                end
            end

            if Vpixx_response == 1 || Vpixx_response == 2

                % Check the Vpixx response buttons to see if a button has been pressed
                Datapixx('RegWrRd');
                kbcheck = dec2bin(Datapixx('GetDinValues'));
    
                if kbcheck(end)   == '1' % left red
    
                    dim_button = 1;
    
%                 elseif kbcheck(end-2) == '1' % left green
%     
%                     dim_button = 1;
%     
%                 elseif kbcheck(end-3) == '1' % left blue
%     
%                     dim_button = 1;
%     
%                 elseif kbcheck(end-1) == '1' % left yellow
%     
%                     dim_button = 1;
        
                end

                % trigger when response to dim flickering (only first frame)
    
                if dim_button == 1 && dim_button_done == 0
    
                    % record response
                    dim.response(seq_i) = 1;

                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);

                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);

                    % Trigger (colour) - response dim
                    Screen('FillRect', window, trigger_response_dim, centeredRect_trigger);

                    % Flip to the screen
                    Screen('Flip', window);

                    % update time
                    time = time + ifi;

                    dim_button_done = 1;
    
                else

                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);

                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);

                    % Trigger (black)
                    Screen('FillRect', window, black, centeredRect_trigger);

                    % Flip to the screen
                    Screen('Flip', window);

                    % update time
                    time = time        + ifi;

                end

            end
    
            % flicker fixation point when required

            if dim.dim_task_trial(seq_i) == 0
            
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);

                % Draw empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);

            elseif dim.dim_task_trial(seq_i) == 1

                if time > dim.dim_onset(seq_i) && time < dim.dim_onset(seq_i) + dim_duration

                    flickering_count = flickering_count + 1;

                    % flickering
                    if flickering_count > 0

                        % Draw fixation point with luminance change
                        Screen('DrawDots', window, [xCenter yCenter],...
                            size_fixation_point, color_fixation_point - luminance_change, [], 2);

                        % Draw empty backgrond dots
                        Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);

                        if flickering_count == frames_flickering
                            flickering_count = -10;
                        end

                    elseif flickering_count <= 0

                        % Draw fixation point without luminance change
                        Screen('DrawDots', window, [xCenter yCenter],...
                            size_fixation_point, color_fixation_point, [], 2);

                        % Draw empty backgrond dots
                        Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);

                    end

                else

                    % Draw fixation point with luminance change
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);

                    % Draw empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);

                end

            end
    
            % draw pixel trigger 
            if ~trigger
    
                % Trigger (colour)
                Screen('FillRect', window, trigger_delay_onset, centeredRect_trigger);

                trigger = 1;
    
            else
    
                % Trigger (black)
                Screen('FillRect', window, black, centeredRect_trigger);
    
            end
        
            % Flip to the screen
            Screen('Flip', window);
        
            % Increment the time
            time = time + ifi;
        
            % update timeout
            if time >= delay_duration
                timeout = 1;
            end
        
        end
        
        %% response

        % Draw fixation point
        Screen('DrawDots', window, [xCenter yCenter],...
            size_fixation_point, color_fixation_point, [], 2);

        % Trigger (colour) - event onset
        Screen('FillRect', window, trigger_response_onset, centeredRect_trigger);

        % Flip to the screen
        Screen('Flip', window);
    
        if Vpixx_response == 2

            % Draw fixation point
            Screen('DrawDots', window, [xCenter yCenter],...
                size_fixation_point, color_fixation_point, [], 2);
            
            % empty backgrond dots
            Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                    
            % Trigger (black)
            Screen('FillRect', window, black, centeredRect_trigger);
            
            % Flip to the screen
            Screen('Flip', window);

            % set variables
            empty_kbcheck   =  '111111111111110000000000';
            last_kbcheck    =  '000111111111110000000000';
            left_red_button =  '111111111111110000000001';

%             ready2start = 0;
%             
%             % start response as soon as stop pressing any button
%             while ~ready2start
% 
%                 % check the Vpixx response buttons
%                 Datapixx('RegWrRd');
%                 kbcheck = dec2bin(Datapixx('GetDinValues'));
% 
%                 if strcmp(kbcheck, empty_kbcheck) == 1
%                     ready2start = 1;
%                 end
% 
%             end
                                                
            % Set the color of filled dot of selected spatial location (grey)
            dotColorMatrix = repmat(grey', 1, number_locations);
                    
            % start at random location (not first location in sequence)
            random_position = 0;
            while ~random_position
                location_response = randi([1 length(sequence)], 1);
                if location_response ~= sequence(1)
                    random_position = 1;
                end
            end

            % location to degrees
            if location_response == 1
                dial_degrees = 0;
            elseif location_response == 2
                dial_degrees = 45;
            elseif location_response == 3
                dial_degrees = 90;
            elseif location_response == 4
                dial_degrees = 135;
            elseif location_response == 5
                dial_degrees = 180;
            elseif location_response == 6
                dial_degrees = 225;
            elseif location_response == 7
                dial_degrees = 270;
            elseif location_response == 8
                dial_degrees = 315;
            end

            % set first location
            dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
            dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                                
            % set variables
            number_responses = 0;

            % while # responses is below length of sequence
            while (number_responses < length(sequence))
            
                % set variables
                select_location = 0;
                update_location = 0;
    
                % color response position
                dotColorMatrix_stim = grey_dark;
                
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);
                
                % empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                
                % Draw the dot to the screen.
                Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
            
                % Trigger (black)
                Screen('FillRect', window, black, centeredRect_trigger);
                
                % Flip to the screen
                Screen('Flip', window);
            
                % Check the dial
                [keyIsDown, secs, keyCode] = KbCheck(dialID);

                if keyCode(cwKey)

                    % update dial location (degrees)
                    dial_degrees = dial_degrees - ts_or_bs_dial_angle;

                    update_location = 1;

                elseif keyCode(ccwKey)

                    % update dial location (degrees)
                    dial_degrees = dial_degrees + ts_or_bs_dial_angle;

                    update_location = 1;

                end

                % keep in range [0 360)
                if ~(dial_degrees < 360)
                    dial_degrees = dial_degrees - 360;
                elseif dial_degrees < 0
                    dial_degrees = dial_degrees + 360;
                end

                % check the Vpixx response buttons
                Datapixx('RegWrRd');
                kbcheck = dec2bin(Datapixx('GetDinValues'));

                % reset all responses as left red
                if kbcheck(end)   == '1'
    
                    kbcheck = left_red_button;
    
%                 elseif kbcheck(end-2) == '1' % left green
%     
%                     kbcheck = left_red_button;
%     
%                 elseif kbcheck(end-3) == '1' % left blue
%     
%                     kbcheck = left_red_button;
%     
%                 elseif kbcheck(end-1) == '1' % left yellow
%     
%                     kbcheck = left_red_button;
            
                end

                % if respone is different than previous response 
                if strcmp(kbcheck, last_kbcheck) ~= 1

                    last_kbcheck = kbcheck;
        
                    if kbcheck(end)   == '1' % left red - % SELECT LOCATION
        
                        select_location = 1;
        
%                     elseif kbcheck(end-2) == '1' % left green
%         
%                         select_location = 1;
%         
%                     elseif kbcheck(end-3) == '1' % left blue
%         
%                         select_location = 1;
%         
%                     elseif kbcheck(end-1) == '1' % left yellow
%         
%                         select_location = 1;
            
                    end

                end

                    
                % update location
                if update_location == 1

                    % update location
                    if dial_degrees < (0 + 22.5) || dial_degrees > (360 - 22.5)
                        location_response = 1;

                    elseif dial_degrees < (45 + 22.5) && dial_degrees > (45 - 22.5)
                        location_response = 2;

                    elseif dial_degrees < (90 + 22.5) && dial_degrees > (90 - 22.5)
                        location_response = 3;

                    elseif dial_degrees < (135 + 22.5) && dial_degrees > (135 - 22.5)
                        location_response = 4;

                    elseif dial_degrees < (180 + 22.5) && dial_degrees > (180 - 22.5)
                        location_response = 5;

                    elseif dial_degrees < (225 + 22.5) && dial_degrees > (225 - 22.5)
                        location_response = 6;

                    elseif dial_degrees < (270 + 22.5) && dial_degrees > (270 - 22.5)
                        location_response = 7;

                    elseif dial_degrees < (315 + 22.5) && dial_degrees > (315 - 22.5)
                        location_response = 8;

                    end

                    % update location of dark_grey dot
                    dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                    dotPositionMatrix_stim = dotPositionMatrix(:, location_response);

                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);
                
                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                        
                    % Draw the dot to the screen.
                    Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
    
                    % Trigger (black)
                    Screen('FillRect', window, black, centeredRect_trigger);

                    % Flip to the screen
                    Screen('Flip', window);

                % select location
                elseif select_location == 1
    
                    % change color of selected location to red
                    dotColorMatrix_stim = red;
        
                    % save response
                    number_responses = number_responses + 1;
                    responses(seq_i, number_responses) = location_response;

                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);
                
                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                        
                    % Draw the dot to the screen.
                    Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
    
                    % Trigger (color) - depending on selected location
                    if location_response == 1
                        trigger_resp = trigger_location1;
                    elseif location_response == 2
                        trigger_resp = trigger_location2;
                    elseif location_response == 3
                        trigger_resp = trigger_location3;
                    elseif location_response == 4
                        trigger_resp = trigger_location4;
                    elseif location_response == 5
                        trigger_resp = trigger_location5;
                    elseif location_response == 6
                        trigger_resp = trigger_location6;
                    elseif location_response == 7
                        trigger_resp = trigger_location7;
                    elseif location_response == 8
                        trigger_resp = trigger_location8;
                    end                        

                    Screen('FillRect', window, trigger_resp, centeredRect_trigger);

                    % Flip to the screen
                    Screen('Flip', window);

                    WaitSecs(0.2);

                end

            end
                            
            % Draw fixation point
            Screen('DrawDots', window, [xCenter yCenter],...
                size_fixation_point, color_fixation_point, [], 2);
        
            % empty backgrond dots
            Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                
            % Trigger
            Screen('FillRect', window, black, centeredRect_trigger);

            % Flip to the screen
            Screen('Flip', window);
        
            % pause after all responses were provided
            WaitSecs(0.5);

            % save temporary responses at the end of trial (backup)
            save([output_folder '/responses_task_backup.mat'], 'responses');

            % save trials with dim fixation point
            save([output_folder '/task_trials_dimfixation_backup.mat'], 'dim');

        elseif Vpixx_response == 1
                                                
            % Set the color of filled dot of selected spatial location (grey)
            dotColorMatrix = repmat(grey', 1, number_locations);
                    
            % start at random location (not first location in sequence)
            random_position = 0;
            while ~random_position
                location_response = randi([1 length(sequence)], 1);
                if location_response ~= sequence(1)
                    random_position = 1;
                end
            end
            
            % set first location
            dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
            dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                                
            % set variables
            number_responses = 0;

            % empty response
            empty_kbcheck = [];
            last_kbcheck = [];
             
            % while # responses is below length of sequence
            while (number_responses < length(sequence))
            
                % set variables
                updated_response = 0;
                wait_location = 0;
    
                % color response position
                dotColorMatrix_stim = grey_dark;
                
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);
                
                % empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                
                % Draw the dot to the screen.
                Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
            
                % Trigger (black)
                Screen('FillRect', window, black, centeredRect_trigger);
                
                % Flip to the screen
                Screen('Flip', window);
            
                % Check the keyboard to see if a button has been pressed
                [keyIsDown, secs, keyCode] = KbCheck;
    
                % check the buttonbox response
                Datapixx('RegWrRd');
                kbcheck = dec2bin(Datapixx('GetDinValues'));

                % empty/last response
                if isempty(empty_kbcheck)
                    empty_kbcheck = zeros(1, length(kbcheck));
                end

                if isempty(last_kbcheck)
                    last_kbcheck = ones(1, length(kbcheck));
                end
                
                % exit task if scape bar has been pressed
                if keyCode(escapeKey)
                    
                    if exit_task == 0
                        exit_task = 1;
                    end
    
                end

                % if kbcheck is different than last_kbcheck OR kbcheck is empty
                if sum(kbcheck ~= last_kbcheck) > 0 || sum(kbcheck ~= empty_kbcheck) == 0

                    % last button pressed
                    last_kbcheck = kbcheck;

                    % response button Vpixx
                    if kbcheck(end)   == '1' % left red - GO LEFT
    
                        if updated_response == 0
        
                            location_response = location_response + 1;
                            if location_response == 9
                                location_response = 1;
                            end
                            dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                            dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                
                            updated_response = 1;
    
                        end
        
                    elseif kbcheck(end-2) == '1' % left green - GO RIGHT
        
                        if updated_response == 0
    
                            location_response = location_response - 1;
                            if location_response == 0
                                location_response = 8;
                            end
                            dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                            dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                
                            updated_response = 1;
    
                        end
        
                    elseif kbcheck(end-4) == '1' % left white - SELECT LOCATION
    
                        if wait_location == 0
        
                            dotColorMatrix_stim = red;
                
                            % save response
                            number_responses = number_responses + 1;
                            responses(seq_i, number_responses) = location_response;
    
                            wait_location = 1;
    
                        end
        
                    end
                
                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, color_fixation_point, [], 2);
                
                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                        
                    % Draw the dot to the screen.
                    Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
                        
                    if wait_location == 1 
        
                        % Trigger (color)

                        if location_response == 1
                            trigger_resp = trigger_location1;
                        elseif location_response == 2
                            trigger_resp = trigger_location2;
                        elseif location_response == 3
                            trigger_resp = trigger_location3;
                        elseif location_response == 4
                            trigger_resp = trigger_location4;
                        elseif location_response == 5
                            trigger_resp = trigger_location5;
                        elseif location_response == 6
                            trigger_resp = trigger_location6;
                        elseif location_response == 7
                            trigger_resp = trigger_location7;
                        elseif location_response == 8
                            trigger_resp = trigger_location8;
                        end                        

                        Screen('FillRect', window, trigger_resp, centeredRect_trigger);
        
                    else
        
                        % Trigger (black)
                        Screen('FillRect', window, black, centeredRect_trigger);
        
                    end
                    
                    % Flip to the screen
                    Screen('Flip', window);
                
                    % update variables                
                    if wait_location == 1                
                        WaitSecs(0.2);
                    end
                                
                end

            end
        
            % pause after all responses were provided
            WaitSecs(0.5);

            % save temporary responses at the end of trial (backup)
            save([output_folder '/responses_task_backup.mat'], 'responses');

            % save trials with dim fixation point
            save([output_folder '/task_trials_dimfixation_backup.mat'], 'dim');
    
        elseif Vpixx_response == 0
    
            % Set the color of filled dot of selected spatial location (grey)
            dotColorMatrix = repmat(grey', 1, number_locations);
                    
            % start at random location (not first location in sequence)
            random_position = 0;
            while ~random_position
                location_response = randi([1 length(sequence)], 1);
                if location_response ~= sequence(1)
                    random_position = 1;
                end
            end
            
            % set first location
            dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
            dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                                
            % set variables
            number_responses = 0;

            % empty response
            empty_keyCode = [];
            last_keyCode = [];
             
            % while # responses is below length of sequence
            while (number_responses < length(sequence))
            
                % set variables
                updated_response = 0;
                wait_location = 0;
    
                % color response position
                dotColorMatrix_stim = grey_dark;
                
                % Draw fixation point
                Screen('DrawDots', window, [xCenter yCenter],...
                    size_fixation_point, color_fixation_point, [], 2);
                
                % empty backgrond dots
                Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                
                % Draw the dot to the screen.
                Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
            
                % Trigger (black)
                Screen('FillRect', window, black, centeredRect_trigger);
                
                % Flip to the screen
                Screen('Flip', window);
            
                % Check the keyboard to see if a button has been pressed
                [keyIsDown,secs, keyCode] = KbCheck;

                % empty/last response
                if isempty(empty_keyCode)
                    empty_keyCode = zeros(1, length(keyCode));
                end

                if isempty(last_keyCode)
                    last_keyCode = ones(1, length(keyCode));
                end

                % if keyCode is different than last_keyCode OR keyCode is empty
                if sum(keyCode ~= last_keyCode) > 0 || sum(keyCode ~= empty_keyCode) == 0

                    % last button pressed
                    last_keyCode = keyCode;
                    
                    % exit task if scape bar has been pressed
                    if keyCode(escapeKey)
                        
                        exit_task = 1;
                        
                    % update location based on response 
                    elseif keyCode(upKey)
                        
                        if updated_response == 0
        
                            location_response = location_response + 1;
                            if location_response == 9
                                location_response = 1;
                            end
                            dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                            dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                
                            updated_response = 1;
    
                        end
                        
                    % update location based on response
                    elseif keyCode(downKey)
                        
                        if updated_response == 0
    
                            location_response = location_response - 1;
                            if location_response == 0
                                location_response = 8;
                            end
                            dotSizePixMatrix_stim = dotSizePixMatrix(location_response);
                            dotPositionMatrix_stim = dotPositionMatrix(:, location_response);
                
                            updated_response = 1;
    
                        end
                
                    % highlight selected location in red
                    elseif keyCode(barKey)
                
                        if wait_location == 0
        
                            dotColorMatrix_stim = red;
                
                            % save response
                            number_responses = number_responses + 1;
                            responses(seq_i, number_responses) = location_response;
    
                            wait_location = 1;
    
                        end
                        
                    end
                
                    % Draw fixation point
                    Screen('DrawDots', window, [xCenter yCenter],...
                        size_fixation_point, white, [], 2);
                
                    % empty backgrond dots
                    Screen('FrameOval', window, ovalColorMatrix, ovalPositionMatrix, penWidthPixels);
                        
                    % Draw the dot to the screen.
                    Screen('DrawDots', window, dotPositionMatrix_stim, dotSizePixMatrix_stim, dotColorMatrix_stim, [], 2);
                        
                    if wait_location == 1 
        
                        % Trigger (color)

                        if location_response == 1
                            trigger_resp = trigger_location1;
                        elseif location_response == 2
                            trigger_resp = trigger_location2;
                        elseif location_response == 3
                            trigger_resp = trigger_location3;
                        elseif location_response == 4
                            trigger_resp = trigger_location4;
                        elseif location_response == 5
                            trigger_resp = trigger_location5;
                        elseif location_response == 6
                            trigger_resp = trigger_location6;
                        elseif location_response == 7
                            trigger_resp = trigger_location7;
                        elseif location_response == 8
                            trigger_resp = trigger_location8;
                        end                        

                        Screen('FillRect', window, trigger_resp, centeredRect_trigger);
        
                    else
        
                        % Trigger (black)
                        Screen('FillRect', window, black, centeredRect_trigger);
        
                    end
                    
                    % Flip to the screen
                    Screen('Flip', window);
                
                    % update variables                
                    if wait_location == 1                
                        WaitSecs(0.2);
                    end

                end
                            
            end
                                
            % pause after all responses were provided
            pause(0.4);
                    
            % save temporary responses at the end of trial (backup)
            save([output_folder '/responses_task_backup.mat'], 'responses');

            % save trials with dim fixation point
            save([output_folder '/task_trials_dimfixation_backup.mat'], 'dim');
    
        end
    
        %% feedback
    
        % if correct/incorrect responses, load image of correct/incorrect tick
        responses_trial = responses(seq_i, 1:length(sequence));
    
        if sum(responses_trial == sequence) == length(sequence)
            image = imread(tick_correct);
            correct_responses = correct_responses + 1;
        else
            image = imread(tick_incorrect);
        end
        
        % Make the image into a texture
        imageTexture = Screen('MakeTexture', window, image);
    
        % Trigger (colour)
        Screen('FillRect', window, trigger_response_onset, centeredRect_trigger);
    
        % Flip screen
        Screen('Flip', window);
    
        time = 0;
    
        while time < feedback_duration
    
            time = time + ifi;
            
            % Draw grey background to the screen.
            Screen('FillRect', window, grey, centeredRect);
    
            % Draw the image to the screen
            Screen('DrawTexture', window, imageTexture, [], [], 0);
                
            % Flip screen
            Screen('Flip', window);
    
        end

    end
        
end

if Vpixx_response == 1 || Vpixx_response == 2

    % close pixel mode
    Datapixx('DisablePixelMode');
    Datapixx('RegWrRd');
    Datapixx('Close');

end

% save responses at the end of the task
save([output_folder '/responses_task.mat'], 'responses');

% save trials with dim fixation point
save([output_folder '/task_trials_dimfixation.mat'], 'dim');

%% prompt text end of the block

if block == 1 || block == 2

    % Draw grey background to the screen.
    Screen('FillRect', window, grey, centeredRect);
    
    % Trigger (black)
    Screen('FillRect', window, black, centeredRect_trigger);
    
    % Flip to the screen
    Screen('Flip', window);
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Now, time for resting a bit', 'center', screenYpixels * 0.25, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
            
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Now, time for resting a bit', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'Please, do not move', 'center', screenYpixels * 0.5, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
            
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'Now, time for resting a bit', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'Please, do not move', 'center', screenYpixels * 0.5, white);
    DrawFormattedText(window, 'We will resume the task in 30 seconds', 'center', screenYpixels * 0.75, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);

else % last block

    % Draw grey background to the screen.
    Screen('FillRect', window, grey, centeredRect);
    
    % Trigger (black)
    Screen('FillRect', window, black, centeredRect_trigger);
    
    % Flip to the screen
    Screen('Flip', window);
    
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'The task is completed', 'center', screenYpixels * 0.25, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
            
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'The task is completed', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'Please, do not move or blink', 'center', screenYpixels * 0.5, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);
            
    % Draw text in the middle of the screen in Courier in white
    Screen('TextSize', window, 60);
    Screen('TextFont', window, 'Courier');
    DrawFormattedText(window, 'The task is completed', 'center', screenYpixels * 0.25, white);
    DrawFormattedText(window, 'Please, do not move', 'center', screenYpixels * 0.5, white);
    DrawFormattedText(window, 'We will get you out of the MEG in 30 seconds', 'center', screenYpixels * 0.75, white);
    Screen('FillRect', window, black, centeredRect_trigger);
    Screen('Flip', window);
    WaitSecs(duration_text);

end


%% Close all

% enable matlab listen keyboard input in command window
ListenChar
