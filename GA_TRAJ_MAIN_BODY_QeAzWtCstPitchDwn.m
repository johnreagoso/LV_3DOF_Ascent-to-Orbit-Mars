function [J_COST_AVG, J_COST_MAX, TRAJ_POPINDV_FINAL] = GA_TRAJ_MAIN_BODY_QeAzWtCstPitchDwn(TRAJPOP_INDV, select_option, distro_option, ga_run_count)
global Qe_StartBase maxQe_add
global Az_StartBase maxAz_sub
global elite_operator
global max_Az Az_base
global max_coast_time
global univ_marker univ_bit
global max_pitchdwn_no pitch_marker


%% Chromosome Length Determinator:
%Dec2Bin: 7 ... 3 digits, 15 ... 4 digits, 31 ... 5 digits, 63 ... 6
%digits, 127 ... 7 digits, 255 ... 8 digits, ... 511 9 digits

Qe_add_integer_no = numel(dec2bin(maxQe_add));
Qe_dec_no         = numel(dec2bin(3)); 

Az_sub_integer_no = numel(dec2bin(maxAz_sub));
Az_dec_no         = numel(dec2bin(3));

coast2_no       = numel(dec2bin(max_coast_time));
coast_dec_no        = numel(dec2bin(3));

pitchdwn_no     = numel(dec2bin(max_pitchdwn_no));

chrom_length = Qe_add_integer_no + Qe_dec_no + Az_sub_integer_no + Az_dec_no + coast2_no + coast_dec_no + pitchdwn_no ; 

%% Select Optimal Parent Population (Rank/Elite Status): Bigger than better.. 
if strcmp(select_option,'Roulette-Wheel') % Roulette Wheel Parent Determination Logic:
 
    fitness_look        = [TRAJPOP_INDV(:).J]; % build an array of the fitness of each adult in population
    [fitness, index]    = sort(fitness_look, 'descend');

    %% Output:
    % disp('highest fitness for this generation is:');
    % disp(fitness(1));
    % disp(TrajPopIndiv(index(1)));    
    
    new = fitness./sum(fitness);
    segment = round(new*100000);
    
    %% Record highest fitness for each generation:
    
     highFit(ga_run_count) = fitness(1);
    
    %% Highest Fitness Plotting
   % hold on;
   % plot(int32(ga_run_count), fitness(1),'o','MarkerSize', 4, 'MarkerEdgeColor','b','MarkerFaceColor','b');

    %% Roulett Wheel Spinning
    try
    [ROULETTE_INDX_GRP_A, ROULETTE_INDX_GRP_B] = RouletteWheel(segment, index, TRAJPOP_INDV);
    catch
        disp('here');
    end

elseif strcmp(select_option,'Rank-Selection') % Rank Selection Determination Logic:
   
    fitness_look     = [TRAJPOP_INDV(:).J]; % build an array of the fitness of each adult in population
    [fitness, index] = sort(fitness_look, 'descend'); 

    for count = 1:1:numel(fitness)    
        F_rank(count) =  min(fitness) + (numel(fitness)-count)*(max(fitness)-min(fitness))/(numel(fitness));      
        % F(count) = max(fitness) - (max(fitness)- min(fitness))*(count-1)/(numel(fitness)); 
    end
        Pi      = F_rank./sum(F_rank);
        segment = round(Pi*1000);

    [ROULETTE_INDX_GRP_A, ROULETTE_INDX_GRP_B] = RouletteWheel(segment, index, TRAJPOP_INDV);

elseif strcmp(select_option,'Tournament-Selection')
   % disp('tournament');
    fitness_look     = [TRAJPOP_INDV(:).J]; % build an array of the fitness of each adult in population
    [~, index] = sort(fitness_look, 'descend'); 
    
    for mm = 1:1:numel(TRAJPOP_INDV)/2-1 
        ParentCandidateA_index = randi(numel(TRAJPOP_INDV));
        ParentCandidateB_index = randi(numel(TRAJPOP_INDV));      
        
        ParentAscore =  TRAJPOP_INDV(ParentCandidateA_index).J;
        ParentBscore =  TRAJPOP_INDV(ParentCandidateB_index).J;
        
        if ParentAscore > ParentBscore
          ROULETTE_INDX_GRP_A(mm) = ParentCandidateA_index;   
        elseif ParentBscore > ParentAscore
          ROULETTE_INDX_GRP_A(mm) = ParentCandidateB_index;
        else
            roll = randi(2);
            if roll ==1
                ROULETTE_INDX_GRP_A(mm) = ParentCandidateA_index;   
            elseif roll ==2
                ROULETTE_INDX_GRP_A(mm) = ParentCandidateB_index;
            else 
                ROULETTE_INDX_GRP_A(mm) = ParentCandidateA_index;    
            end
        end
    end  

    for mm = 1:1:numel(TRAJPOP_INDV)/2-1 
        ParentCandidateA_index = randi(numel(TRAJPOP_INDV));
        ParentCandidateB_index = randi(numel(TRAJPOP_INDV));      

        ParentAscore =  TRAJPOP_INDV(ParentCandidateA_index).J;
        ParentBscore =  TRAJPOP_INDV(ParentCandidateB_index).J;

    if ParentAscore > ParentBscore
      ROULETTE_INDX_GRP_B(mm) = ParentCandidateA_index;   
    elseif ParentBscore > ParentAscore
      ROULETTE_INDX_GRP_B(mm) = ParentCandidateB_index;
    else
        roll = randi(2);
        if roll ==1
            ROULETTE_INDX_GRP_B(mm) = ParentCandidateA_index;   
        elseif roll ==2
            ROULETTE_INDX_GRP_B(mm) = ParentCandidateB_index;
        else 
            ROULETTE_INDX_GRP_B(mm) = ParentCandidateA_index;  
        end
    end
    end  
end

%% Elite Fitness Guarantee Logic- Tom Brady/Gisele Code Logic
%  Ensures that the top-fit parent reproduces itself for each iteration. 

if strcmp(select_option,'Roulette-Wheel')           % Roulette Wheel Parent Determination Logic:
    ROULETTE_INDX_GRP_A(1) = index(1);
    ROULETTE_INDX_GRP_B(1) = index(2);

elseif strcmp(select_option,'Rank-Selection')       % Rank Selection Determination Logic:
    ROULETTE_INDX_GRP_A(1) = index(1);
    ROULETTE_INDX_GRP_B(1) = index(2);

elseif strcmp(select_option,'Tournament-Selection') % Tournament Determination Logic:
    ROULETTE_INDX_GRP_A(end+1) = index(1);
    ROULETTE_INDX_GRP_B(end+1) = index(2);
end

%% New Child Mating (x2):
count = 1;

for kk = 1:1:numel(ROULETTE_INDX_GRP_A)

        try
            mate1_index         = ROULETTE_INDX_GRP_A(kk);
            
            Qe_integer_chromo1  = dec2bin(floor(TRAJPOP_INDV(mate1_index).Qe - Qe_StartBase),  Qe_add_integer_no);  %4 -> 15, 5 -> 31
            qe_dec_input        = find(univ_marker(:,2) == round(TRAJPOP_INDV(mate1_index).Qe - floor(TRAJPOP_INDV(mate1_index).Qe),2), 1, 'first'); 
            Qe_dec_chromo1      = univ_bit{qe_dec_input};
            
            pitch_idx          = find(pitch_marker(:) == TRAJPOP_INDV(mate1_index).pitchdown);
            pitchdwn_chromo1    = dec2bin(pitch_idx, pitchdwn_no);
            
            Az_integer_chromo1  = dec2bin(floor(Az_StartBase - TRAJPOP_INDV(mate1_index).Az),  Az_sub_integer_no);  %4 -> 15, 5 -> 31
            az_dec_input        = find(univ_marker(:,2) == round(TRAJPOP_INDV(mate1_index).Az - floor(TRAJPOP_INDV(mate1_index).Az),2), 1, 'first'); 
            Az_dec_chromo1      = univ_bit{az_dec_input};
         
            %coast2_chromo1      = dec2bin(TrajPopIndiv(mate1_index).coast2, coast2_no);
            coast2_chromo1      = dec2bin(floor(TRAJPOP_INDV(mate1_index).coast2), coast2_no);
            
            coast_dec_input     = find(univ_marker(:,2) == round(TRAJPOP_INDV(mate1_index).coast2  - floor(TRAJPOP_INDV(mate1_index).coast2),2), 1, 'first'); 
            coast_dec_chromo1   = univ_bit{coast_dec_input};
                        
        catch
            disp('prob here line 136:');
        end
        
            mate1_chromo        = [Qe_integer_chromo1, Qe_dec_chromo1, Az_integer_chromo1, Az_dec_chromo1, coast2_chromo1, coast_dec_chromo1, pitchdwn_chromo1];

        try
            mate2_index         = ROULETTE_INDX_GRP_B(kk);
            
            Qe_integer_chromo2  = dec2bin(floor(TRAJPOP_INDV(mate2_index).Qe - Qe_StartBase),  Qe_add_integer_no);  %4 -> 15, 5 -> 31
            qe_dec_input           = find(univ_marker(:,2) == round(TRAJPOP_INDV(mate2_index).Qe - floor(TRAJPOP_INDV(mate2_index).Qe),2), 1, 'first'); 
            Qe_dec_chromo2      = univ_bit{qe_dec_input};

            pitch_idx          = find(pitch_marker(:) == TRAJPOP_INDV(mate2_index).pitchdown);
            pitchdwn_chromo2    = dec2bin(pitch_idx, pitchdwn_no);
            
            Az_integer_chromo2  = dec2bin(floor(Az_StartBase - TRAJPOP_INDV(mate2_index).Az),  Az_sub_integer_no);  %4 -> 15, 5 -> 31
            az_dec_input           = find(univ_marker(:,2) == round(TRAJPOP_INDV(mate2_index).Az - floor(TRAJPOP_INDV(mate2_index).Az),2), 1, 'first'); 
            Az_dec_chromo2      = univ_bit{az_dec_input};
         
            % coast2_chromo2    = dec2bin(TrajPopIndiv(mate2_index).coast2, coast2_no);
            coast2_chromo2      = dec2bin(floor(TRAJPOP_INDV(mate2_index).coast2), coast2_no);

            coast_dec_input     = find(univ_marker(:,2) == round(TRAJPOP_INDV(mate2_index).coast2  - floor(TRAJPOP_INDV(mate2_index).coast2),2), 1, 'first'); 
            coast_dec_chromo2   = univ_bit{coast_dec_input};
           % ballast_chromo2    = dec2bin(TrajPopIndiv(mate2_index).Ballast, ballast_no);
        
        catch
            disp('prob here');
        end

            mate2_chromo        = [Qe_integer_chromo2, Qe_dec_chromo2, Az_integer_chromo2, Az_dec_chromo2, coast2_chromo2, coast_dec_chromo2, pitchdwn_chromo2];

    if numel(mate2_chromo) ~= 25 || numel(mate1_chromo) ~=25
  %      disp('here');
    end

    %% Allele Distribution
    if strcmp(distro_option, 'Uniform')
        try  % Uniform 
        for jj = 1:1:numel(mate1_chromo)-1
            child1(1,jj)   = mate1_chromo(jj);
            child1(1,jj+1) = mate2_chromo(jj+1);
        end
        catch
            disp('problem here- line201');
            pause;
        end

    elseif strcmp(distro_option, 'Single-Segment')
        
        chromo_split = randi([2, numel(mate1_chromo)-1]);
        mate1_chromo_split = mate1_chromo(1:chromo_split);
        mate2_chromo_split = mate2_chromo(1+chromo_split:end);
        child1 = [mate1_chromo_split, mate2_chromo_split];

    elseif strcmp(distro_option, 'Double-Segment')

        chromo_split = randi([2, numel(mate1_chromo)-1]);
        mate1_chromo_split = mate1_chromo(1:chromo_split);
        mate2_chromo_split = mate2_chromo(1+chromo_split:end);
        child1 = [mate1_chromo_split, mate2_chromo_split];

    end

    %% Mutation Operations
    % Dice Roll: 
    if kk ~= 1
                 
%         if ga_run_count > mutateHelp && mean(highFit(ga_run_count - mutateHelp:ga_run_count)) == highFit(ga_run_count)
%             disp('Mutation increased !!');
%             
%             for iii = 1:1:numel(child1)
%                 dice = randi(numel(child1));
%                 
%                 if dice == 1 || dice == 2 || dice == 3 % if dice-roll equals '1' from previous line, proceed with mutating this specific bit % if dice-roll equals '1' from previous line, proceed with mutating this specific bit
%             %   disp('mutation');
%                     if strcmp(child1(iii),'1') == 1
%                         child1(iii) = '0';
%                     else
%                         child1(iii) = '1';
%                     end
%                 end
%             end
%                         
%         else
            
            for iii = 1:1:numel(child1)
                dice = randi(numel(child1));

                if dice == 1 || dice == 2 || dice == 3  % if dice-roll equals '1' from previous line, proceed with mutating this specific bit
            %   disp('mutation');
                    if strcmp(child1(iii),'1') == 1
                        child1(iii) = '0';
                    else
                        child1(iii) = '1';
                    end
                end

            end
    %    end
    %     mutate_dice_roll = randi(dice_roll_int);
    %     mutate_yes = 1;
    % 
    %     if mutate_yes == mutate_dice_roll  %yes
    %       %  disp('mutation-1 occuring');
    %         mutate_location  = randi(numel(child1));
    % 
    %         if strcmp(child1(mutate_location),'1') == 1
    %            child1(mutate_location) = '0';
    %         else
    %            child1(mutate_location) = '1';
    %         end
    %     
    %         mutate_location  = randi(numel(child1));
    % 
    %         if strcmp(child1(mutate_location),'1') == 1
    %            child1(mutate_location) = '0';
    %         else
    %            child1(mutate_location) = '1';
    %        
    %         end
    %     end
    end

    %% Allele Distribution
    if strcmp(distro_option,'Uniform')

        try  % Uniform 
        for jj = 1:1:numel(mate1_chromo)-1
            child2(1,jj)   = mate2_chromo(jj);
            child2(1,jj+1) = mate1_chromo(jj+1);
        end
        catch
            disp('problem here');
            pause;
        end
    elseif strcmp(distro_option, 'Single-Segment')

        chromo_split = randi([2, numel(mate1_chromo)-1]);
        mate1_chromo_split = mate1_chromo(1:chromo_split);
        mate2_chromo_split = mate2_chromo(1+chromo_split:end);
        child2 = [mate1_chromo_split, mate2_chromo_split];

    elseif strcmp(distro_option, 'Double-Segment')

        chromo_split = randi([2, numel(mate1_chromo)-1]);
        mate1_chromo_split = mate1_chromo(1:chromo_split);
        mate2_chromo_split = mate2_chromo(1+chromo_split:end);
        child2 = [mate1_chromo_split, mate2_chromo_split];

    end

    %% Mutation Operations
    if kk ~= 1
%%
%         if ga_run_count > mutateHelp && double(mean(highFit(ga_run_count - mutateHelp:ga_run_count))) == double(highFit(ga_run_count))
%                    
%             disp('Mutation increased !!');
%                    
%             for iii = 1:1:numel(child2)
%                 dice = randi(numel(child2));
% 
%                 if dice == 1 || dice == 2 || dice == 3 % if dice-roll equals '1' from previous line, proceed with mutating this specific bit
%             %   disp('mutation');
%                     if strcmp(child2(iii),'1') == 1
%                         child2(iii) = '0';
%                     else
%                         child2(iii) = '1';
%                     end
%                 end
%             end
%                         
%         else
            
            for iii = 1:1:numel(child2)
                dice = randi(numel(child2));

                if dice == 1  % if dice-roll equals '1' from previous line, proceed with mutating this specific bit
            %   disp('mutation');
                    if strcmp(child2(iii),'1') == 1
                        child2(iii) = '0';
                    else
                        child2(iii) = '1';
                    end
                end
            end
    %    end    

%%
    %     mutate_dice_roll = randi(dice_roll_int);
    %     mutate_yes = 1;
    % 
    %     if mutate_yes == mutate_dice_roll  %yes
    %       %  disp('mutation-2 occuring');
    %         mutate_location  = randi(numel(child2));
    % 
    %         if strcmp(child2(mutate_location),'1') == 1
    %            child2(mutate_location) = '0';
    %         else
    %            child2(mutate_location) = '1';
    %         end
    %         
    %          mutate_location  = randi(numel(child2));
    % 
    %         if strcmp(child2(mutate_location),'1') == 1
    %            child2(mutate_location) = '0';
    %         else
    %            child2(mutate_location) = '1';
    %         end  
    %     end

    end

% Generate new children:
%% Child-1 production:
    Qe_add = bin2dec(child1(1:5));     
    Qe_dec = univ_marker(bin2dec(child1(6:7))+1,2);    
    Az_sub = bin2dec(child1(8:13));
    Az_dec = univ_marker(bin2dec(child1(14:15))+1,2);
    
    coast2     = bin2dec(child1(16:24));
    coast2_dec = univ_marker(bin2dec(child1(25:26))+1,2);
    
    pitch_idx = bin2dec(child2(27:32));
    if pitch_idx == 0; pitch_idx = 1; end

    pitchdown = pitch_marker(pitch_idx);

    ChildPop(count).Qe = Qe_StartBase + Qe_add + Qe_dec;
    ChildPop(count).Az = Az_StartBase - abs(Az_sub + Az_dec);
    ChildPop(count).coast2 = coast2 + coast2_dec;
    ChildPop(count).pitchdown = pitchdown;

    
%% Child-2 production:
    count = count + 1;

    Qe_add = bin2dec(child2(1:5));     
    Qe_dec = univ_marker(bin2dec(child2(6:7))+1,2);    
    Az_sub = bin2dec(child2(8:13));
    Az_dec = univ_marker(bin2dec(child2(14:15))+1,2);
    
    coast2     = bin2dec(child2(16:24));
    coast2_dec = univ_marker(bin2dec(child2(25:26))+1,2);
    
    pitch_idx = bin2dec(child2(27:32));
    if pitch_idx == 0; pitch_idx = 1; end
    
    pitchdown = pitch_marker(pitch_idx);
       
    ChildPop(count).Qe        = Qe_StartBase + Qe_add + Qe_dec;
    ChildPop(count).Az        = Az_StartBase - abs(Az_sub + Az_dec);
    ChildPop(count).coast2    = coast2 + coast2_dec;
    ChildPop(count).pitchdown = pitchdown;

    count = count + 1;
end

disp('new population complete');
for kk = 1:1:numel(ChildPop)

    TRAJ_POPINDV_FINAL(kk).Qe      = ChildPop(kk).Qe;
    TRAJ_POPINDV_FINAL(kk).Az      = ChildPop(kk).Az;
    TRAJ_POPINDV_FINAL(kk).coast2  = ChildPop(kk).coast2;
    TRAJ_POPINDV_FINAL(kk).pitchdown = ChildPop(kk).pitchdown;
    
    x_inputs(1) = TRAJ_POPINDV_FINAL(kk).Qe;
    x_inputs(2) = TRAJ_POPINDV_FINAL(kk).Az;
    x_inputs(3) = TRAJ_POPINDV_FINAL(kk).coast2;  
    x_inputs(4) = TRAJ_POPINDV_FINAL(kk).pitchdown;

    % disp(x_inputs);
    output = GA_OrbitInsert_Call3DOF_QeAzCstTimePitch_LaunchFlyout(x_inputs);
    
    J           = output(1);              
    e           = output(2);
    perigee_alt = output(3);        
    incl        = output(4); 

    if output(5) > 360;        output(5) = output(5)-360;    end     
    raan        = output(5);
    
    apogee_alt  = output(6);        
    sma         = output(7);

    TRAJ_POPINDV_FINAL(kk).J           = 1/J;
    TRAJ_POPINDV_FINAL(kk).e           = e;
    TRAJ_POPINDV_FINAL(kk).perigee_alt = perigee_alt;
    TRAJ_POPINDV_FINAL(kk).apogee_alt  =  apogee_alt;    
    TRAJ_POPINDV_FINAL(kk).incl        = incl;    
    TRAJ_POPINDV_FINAL(kk).raan        = raan;  
    TRAJ_POPINDV_FINAL(kk).sma         = sma;  

    J_cost(kk) = 1/J;
    
end

%% Elite Operator: 
% Ensures that the fittest population member of the preceding generation is carried over with 
% the subsequent generation population.

    if strcmp(elite_operator,'yes') == 1 && TRAJPOP_INDV(index(1)).J > TRAJ_POPINDV_FINAL(1).J
        try
            TRAJ_POPINDV_FINAL(1) = TRAJPOP_INDV(index(1));      
        catch
            disp('here');
        end
        
        J_cost(1) = TRAJPOP_INDV(index(1)).J;
    else
        roll                     = randi([2,numel(TRAJ_POPINDV_FINAL)],1);
        TRAJ_POPINDV_FINAL(roll) = TRAJPOP_INDV(index(1));
        J_cost(roll)             = TRAJ_POPINDV_FINAL(roll).J;
    end

J_COST_AVG = mean(J_cost);
J_COST_MAX = max(J_cost);

end
 
function [dummy1, dummy2] = RouletteWheel(segment, index, TrajPopIndiv)

%   Build Roulette wheel tabs:    
    array(1,1) = 0.0;
    array(1,2) = segment(1);

    for ii = 2:1:numel(segment)
        array(ii,1) = array(ii-1,2) + 1;
        array(ii,2) = array(ii,1) + segment(ii);   
    end

for kk = 1:1:numel(TrajPopIndiv)/2
    rnd_integer = randi(sum(segment)); % dice roll

    for jj = 1:1:numel(array)
        if  rnd_integer >= array(jj,1) && rnd_integer <= array(jj,2)
            break;
        end
    end
    roulette_index_GrpA(kk) = index(jj);
end
    dummy1 =  roulette_index_GrpA;

for kk = 1:1:numel(TrajPopIndiv)/2
    rnd_integer = randi(sum(segment)); % dice roll
    
    for jj = 1:1:numel(array)
        if  rnd_integer >= array(jj,1) && rnd_integer <= array(jj,2)
            break;
        end
    end
    roulette_index_GrpB(kk) = index(jj);
end
    dummy2 =  roulette_index_GrpB;
end























