%% Trajectory Shaping Genetic Algorithm 
function check = GA_ATO_EXEC_LaunchFlyoutwPtchDwn(select_no, distro_option_output)

fclose all;
% clc;
clearvars -except select_no distro_option_output     %clear all;

%% UI Inputs:
str1    = {'Select preferred GA Selection option:'};%; 'single selection'};
S1      = {'Roulette-Wheel';'Tournament-Selection'; 'Rank-Selection'};
%select_output  = listdlg('PromptString', str1, 'ListSize', [200 200], 'ListString', S1, 'SelectionMode','single');

select = S1{select_no};
%select = S1{3};

str2    = {'Select preferred GA Allele Crossover option:'};%; 'single selection'};
S2      = {'Uniform'; 'Single-Segment'; 'Double-Segment'};
% distro_option_output = listdlg('PromptString', str2, 'ListSize', [200 200], 'ListString', S2, 'SelectionMode','single');

distro_option = S2{distro_option_output};
%distro_option = S2{1};

%% Global List for Computations:
global input_list; input_list = {};
global dice_roll_int
global max_ballast max_coast_time elite_operator 
global maxQe_add Qe_StartBase Az_StartBase maxAz_sub
global payload %orbitMarker TperiodSecs 
global univ_marker univ_bit
global target_LMO_alt;
global target_LMO_incl;
global max_pitchdwn max_pitchdwn_no pitch_marker

%% Mars Target:

%target_LMO_alt  = 343.0; % km 
%target_LMO_incl = 25.0;  % deg

target_LMO_alt  = 380.0; % km 
target_LMO_incl = 27.0;  % deg

%% Universal markers for Qe, Az decimal quarters:   
univ_marker = [0 0; 1 0.25;  2 0.50; 3 0.75];
univ_bit    = {'00';'01';'10';'11'};

%% Payload (in kg):
 payload = 16;

%% Orbit Determine Triggers for Main Body:
%orbitMarker = 0;    TperiodSecs = 0.0;
%highFit = 0.0;
%mutateHelp = 20.00;

elite_operator = 'yes'; 
%Dec2Bin: 7 ... 3 digits, 15 ... 4 digits, 31 ... 5 digits, 63 ... 6
%digits, 127 ... 7 digits, 255 ... 8 digits, 511 ... 9 digits

%% Min Qe 
 Qe_StartBase   = 58.0;   
 %max_pitchdwn = 0.20;
 max_pitchdwn_no = 63;

 maxQe_add      = 31.0;     % 5 bit
 Az_StartBase   = 90.0; 
 maxAz_sub      = 63.0;     % 6 bit
 max_ballast    = 0.0;   % 7 bit
 max_coast_time = 511;      % 9 bit
                 % Qe dec   % 2 bit
                 % Az dec   % 2 bit
                   % ===== 31 bit total =====

pitch_marker = rescale(1:1:max_pitchdwn_no, 0, 0.2);                   

 %% Launch Vehicle Mass Metrics
dice_roll_int    = 2; 
PopulationNumber = 100;
J_COST_MAX       = 0.00;

%% Main Body Below to Run Multiple GA 300 generation runs

    %Clock-DTG Retrieve
    DTG = clock;
    DTG_string = [num2str(DTG(2)),'_', num2str(DTG(3)),'_', num2str(DTG(1)),'_', ...
    num2str(DTG(4)),'_', num2str(DTG(5)),'_', num2str(round(DTG(6)))];
    disp(DTG_string);
    
    %% Mat String for saving results:
    mat_string = ['GA_SpecsMAT','_',num2str(PopulationNumber),'Pop','_','Qe_start',num2str(Qe_StartBase), '_', num2str(payload),'kg', ...
        '_',select,'_', distro_option, '_Dice_1outof',num2str(dice_roll_int) ,'_','DTG_', DTG_string];

    %% Inhabit Initial Population:

    J_COST_MAX = 100; 
    score_record = 20;
    while max(score_record) <= 250

        for ii = 1:1:PopulationNumber

            az_sub_guess = round(90 - asind(cosd(target_LMO_incl)/cosd(18.38)));

            Qe_add     = (randi([15, maxQe_add],   1));        % This is for the Qe-add to be added to the Qe-base .. 

            pitchdwn   = pitch_marker(randi([1, max_pitchdwn_no], 1));

            Az_sub     = az_sub_guess;
            coast2_add = (randi([400, max_coast_time], 1));
            ballast    = (randi([0, 0],    1));       
            
            Qe_dec_int = (randi([0,3],1));
            Az_dec_int = (randi([0,3],1));
            Coast2_dec_int = (randi([0,3],1));

            if Qe_dec_int == 0;     Qe_dec = 0.00;            elseif Qe_dec_int == 1; Qe_dec = 0.25;
            elseif Qe_dec_int == 2; Qe_dec = 0.50;            else;                   Qe_dec = 0.75;
            end

            if Az_dec_int == 0;     Az_dec = 0.00;            elseif Az_dec_int == 1; Az_dec = 0.25;
            elseif Az_dec_int == 2; Az_dec = 0.50;            else;                   Az_dec = 0.75;
            end

            if Coast2_dec_int == 0;     Coast2_dec = 0.00;    elseif Coast2_dec_int == 1; Coast2_dec = 0.25;
            elseif Coast2_dec_int == 2; Coast2_dec = 0.50;    else;                       Coast2_dec = 0.75;
            end

            % add the Qe addition plus Qe decrement (0.25 or 0.50 or 0.75):
            Qe_integer = Qe_StartBase + Qe_add;      
            Qe = Qe_integer + Qe_dec;

            coast2 = coast2_add + Coast2_dec; 

            % subtract the Az subtraction value: 
            Az = Az_StartBase - (Az_sub + Az_dec);      

            x_inputs(1) = Qe;                       x_inputs(2) = Az;           
            x_inputs(3) = coast2;                   x_inputs(4) = pitchdwn;
            
            disp(x_inputs);
            %x_inputs = [82.00         70.00        489.00          0.18];
            output = GA_OrbitInsert_Call3DOF_QeAzCstTimePitch_LaunchFlyout(x_inputs);

            score       = output(1);  % this is a quantitative representation of the traj 'miss', bigger J, worse orbit                 
            e           = output(2);           
            perigee_alt = output(3);                   
            incl        = output(4);    
            apogee_alt  = output(6);
            sma         = output(7);

            if output(5) > 360; output(5) = output(5)-360; end
            raan        = output(5);

            score_record(ii) = perigee_alt;

            TRAJPOP_INDV(ii).Qe = Qe;
            TRAJPOP_INDV(ii).pitchdown = pitchdwn;
            TRAJPOP_INDV(ii).Az = Az;
            TRAJPOP_INDV(ii).coast2 = coast2;
            TRAJPOP_INDV(ii).Ballast = ballast;
            TRAJPOP_INDV(ii).J = 1/score; 
            TRAJPOP_INDV(ii).e = e;
            TRAJPOP_INDV(ii).perigee_alt = perigee_alt; 
            TRAJPOP_INDV(ii).apogee_alt  = apogee_alt;
            TRAJPOP_INDV(ii).incl = incl;
            TRAJPOP_INDV(ii).raan = raan;
            TRAJPOP_INDV(ii).sma  = sma;

            J_COST(ii) = 1/score;

            J_COST_AVG = mean(J_COST);
            J_COST_MAX = max(J_COST);

            disp('========================================================================');
            fprintf('Initial Population Mean Fitness Value: %f\n ',J_COST_AVG); 
            disp('========================================================================');
            fprintf('Initial Population Max Fitness Value: %f\n  ',J_COST_MAX); 
            disp('========================================================================');
            disp(TRAJPOP_INDV(ii));

        end
    end

    ga_run_count = 1;
    
    while ga_run_count <= 200
        fprintf('GA Run Count: %i\n', ga_run_count);
        disp(ga_run_count);

        %[J_cost_avg, J_cost_max(ga_run_count), TrajPopIndivNew ] = ...
        %    GA_TrajShapeMainBody_QeAzWtCst(TrajPopIndiv, select, distro_option, ga_run_count);

        [J_COST_AVG, J_COST_MAX(ga_run_count), TRAJPOP_INDV_NEW ] = ...
            GA_TRAJ_MAIN_BODY_QeAzWtCstPitchDwn(TRAJPOP_INDV, select, distro_option, ga_run_count);

        %score_array(ga_run_count) = J_cost_max;
        %Jcost_super = rescale(score_array);

        GA_Specs{ga_run_count} = {TRAJPOP_INDV_NEW  J_COST_AVG, J_COST_MAX};
        eval(['save',' ',mat_string,' ','GA_Specs']);
      
        TRAJPOP_INDV = TRAJPOP_INDV_NEW;
        disp('================================================');        disp('================================================');
        disp('GA-run complete');
        fprintf('GA-run cost-average: %0.4f\n ',J_COST_AVG); 
        disp('================================================');
        fprintf('GA-run max: %0.4f\n ',J_COST_MAX); 
        disp('================================================');        disp('================================================');

        if numel(J_COST_MAX) > 50
            if J_COST_MAX(ga_run_count) == J_COST_MAX(ga_run_count-30)
                break;
            end
        end

        ans1 = GA_Specs{ga_run_count}{1};
        disp(ans1(1));
        ga_run_count = ga_run_count + 1;

    end

clock_display = clock(); 
disp(clock_display); 



