%% Sample response processing
% developed by Hamid Karimi-Rouzbahani on 15/June/2022

clc
clear all;

subjects=[10] ; % subjects you want the include in analysis
percentage_target_cond=[0.06 0.12 0.24 0.5]; % Frequency of targets across conditions
Testing_blocks=[1:10]; % blocks you want to include in analysis
Blocks_per_condition = 10;

dirs=dir();
% or Determine where the data is stored on PC
% dirs=dir('C:\');


%% Data preparation
for Subj=subjects
    for blk=Testing_blocks
        for cndss=percentage_target_cond
            correct_reaction_times_att=0;
            
            for i=3:size(dirs,1)
                if strcmp(dirs(i).name(end-3:end),'.mat')
                    if strcmp(dirs(i).name(1:19+length(num2str(Subj))+length(num2str(blk))+1),['Subj_',num2str(Subj),'_Blk_',num2str(blk),'_Freq_',sprintf('%.2f',cndss)])
                        load(dirs(i).name);
                        Targ_Freq_Condition_blk=str2double(dirs(i).name(end-19:end-16));
                    end
                end
            end
            
            mean_sampling_time=1./60;
            for dot_num=1:Num_moving_dots*Trials_per_block
                tr=ceil(dot_num./Num_moving_dots);
                dot_in_trial=dot_num-(tr-1).*Num_moving_dots;
                
                if ~isempty(find(key_pressed1(dot_in_trial,:,tr),1))
                    
                    key_press_sample=find(key_pressed1(dot_in_trial,:,tr), 1, 'first');
                    if isnan(distance_traj1(dot_num,key_press_sample))
                        distance_traj1(dot_num,key_press_sample)=3000;
                    end
                    dist_relative_to_boundary(dot_in_trial,tr)=distance_traj1(dot_num,key_press_sample)-hitting_border_distance;
                else
                    dist_relative_to_boundary(dot_in_trial,tr)=nan;
                end
                distance_change_per_sample(dot_in_trial,tr)=(distance_traj1(dot_num,appearance_time(dot_in_trial,tr)+10)-distance_traj1(dot_num,appearance_time(dot_in_trial,tr)+20))./(11);
                
                if ~isempty(find(key_pressed2(dot_in_trial,:,tr),1))
                    
                    key_press_sample2=find(key_pressed2(dot_in_trial,:,tr), 1, 'first' );
                    if isnan(distance_traj2(dot_num,key_press_sample2))
                        distance_traj2(dot_num,key_press_sample)=3000;
                    end
                    dist_relative_to_boundary2(dot_in_trial,tr)=distance_traj2(dot_num,key_press_sample2)-hitting_border_distance;
                else
                    dist_relative_to_boundary2(dot_in_trial,tr)=nan;
                end
                distance_change_per_sample2(dot_in_trial,tr)=(distance_traj2(dot_num,appearance_time2(dot_in_trial,tr)+10)-distance_traj2(dot_num,appearance_time2(dot_in_trial,tr)+20))./(11);
            end
            
            
            distance_change_per_sample(distance_change_per_sample<0)=mean(distance_change_per_sample(distance_change_per_sample>0));
            distance_change_per_sample2(distance_change_per_sample2<0)=mean(distance_change_per_sample2(distance_change_per_sample2>0));
            
            reaction_times=((-dist_relative_to_boundary)./distance_change_per_sample).*mean_sampling_time;
            reaction_times2=((-dist_relative_to_boundary2)./distance_change_per_sample2).*mean_sampling_time;
            %% Behavioural Performance
            
            % attended
            tp_att=0;
            tn_att=0;
            fp_F_att=0;
            fp_S_att=0;
            fp_T_att=0;
            fn_att=0;
            
            g=0;
            for dot_num=1:Num_moving_dots*Trials_per_block
                tr=ceil(dot_num./Num_moving_dots);
                dot_in_trial=dot_num-(tr-1).*Num_moving_dots;
                
                
                if sum(dot_in_trial==top_events(:,tr))==1 && dot_color(dot_in_trial,tr)==Cued_color_in_block(Subj,blk)
                    g=g+1;
                    if isnan(reaction_times(dot_in_trial,tr)) && (top_events(tr)~=top_targets(tr))
                        tn_att=tn_att+1;    % number of non-target events with no resp;
                    elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                        fp_F_att=fp_F_att+1;    % number of non-target events with fast resp;
                    elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)~=top_targets(tr) && (reaction_times(dot_in_trial,tr)>=0)
                        fp_S_att=fp_S_att+1;    % number of non-target events with Slow resp;
                    elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && (reaction_times(dot_in_trial,tr)<0)
                        fp_T_att=fp_T_att+1;    % number of target events with Too early resp;
                    elseif isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr)
                        fn_att=fn_att+1;    % number of target events with no resp;
                    elseif ~isnan(reaction_times(dot_in_trial,tr)) && top_events(tr)==top_targets(tr) && reaction_times(dot_in_trial,tr)>0
                        tp_att=tp_att+1;    % number of target events with resp;
                        correct_reaction_times_att=correct_reaction_times_att+reaction_times(dot_in_trial,tr);
                    end
                end
                
                if sum(dot_in_trial==top_events2(:,tr))==1 && dot_color2(dot_in_trial,tr)==Cued_color_in_block(Subj,blk)
                    g=g+1;
                    if isnan(reaction_times2(dot_in_trial,tr)) && (top_events2(tr)~=top_targets2(tr))
                        tn_att=tn_att+1;
                    elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)
                        fp_F_att=fp_F_att+1;
                    elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)~=top_targets2(tr) && (reaction_times2(dot_in_trial,tr)>=0)
                        fp_S_att=fp_S_att+1;
                    elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && (reaction_times2(dot_in_trial,tr)<0)% || reaction_times2(dot_in_trial,tr)>time_to_touch_the_obstacle2(dot_in_trial,tr))
                        fp_T_att=fp_T_att+1;
                    elseif isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr)
                        fn_att=fn_att+1;
                    elseif ~isnan(reaction_times2(dot_in_trial,tr)) && top_events2(tr)==top_targets2(tr) && reaction_times2(dot_in_trial,tr)>0 %&& reaction_times2(dot_in_trial,tr)<time_to_touch_the_obstacle2(dot_in_trial,tr)
                        tp_att=tp_att+1;
                        correct_reaction_times_att=correct_reaction_times_att+reaction_times2(dot_in_trial,tr);
                    end
                end
            end
            
            fp_att=fp_F_att+fp_S_att+fp_T_att;
            correct_reaction_times_att=correct_reaction_times_att./tp_att;
            % Removed the unattended dots for simplicity of the data
            
            [~,cond]=ismember(cndss,percentage_target_cond);
            if Targ_Freq_Condition_blk==cndss
                % Accuracy
                Data{cond,1}(blk,Subj)=(tp_att+tn_att)./(sum(top_events>0)+sum(top_events2>0));
                
                % Hit rate
                Data{cond,2}(blk,Subj)=(tp_att)./(tp_att+fn_att);
                
                % True negative rate
                Data{cond,3}(blk,Subj)=(tn_att)./(tn_att+fp_att);
                
                % False alarm
                Data{cond,4}(blk,Subj)=(fp_att)./(fp_att+tn_att);
                
                % Miss
                Data{cond,5}(blk,Subj)=(fn_att)./(tp_att+fn_att);
                
                % Dprime
                Data{cond,6}(blk,Subj)=Data{cond,2}(blk,Subj)-Data{cond,4}(blk,Subj);
                
                % Reaction time
                Data{cond,7}(blk,Subj)=correct_reaction_times_att;
            else
                Data{cond,1}(blk,Subj)=nan;
                Data{cond,2}(blk,Subj)=nan;
                Data{cond,3}(blk,Subj)=nan;
                Data{cond,4}(blk,Subj)=nan;
                Data{cond,5}(blk,Subj)=nan;
                Data{cond,6}(blk,Subj)=nan;
                Data{cond,7}(blk,Subj)=nan;
            end
        end
    end
end
%% Saving data as Excel file for analysis
for Subj=subjects
    for cond=1:size(Data,1)
        if ~ismember(Subj,subjects)
            Data{cond,1}(:,Subj)=nan;
            Data{cond,2}(:,Subj)=nan;
            Data{cond,3}(:,Subj)=nan;
            Data{cond,4}(:,Subj)=nan;
            Data{cond,5}(:,Subj)=nan;
            Data{cond,6}(:,Subj)=nan;
            Data{cond,7}(:,Subj)=nan;
        end
    end
    cond=1;
    False_alarm_condition=Data{cond,4}(:,Subj); % False_alarm in condition
    False_alarm_rate=nanmean(False_alarm_condition);
    Reaction_time_condition=Data{cond,7}(:,Subj); % reaction time in condition
    Mean_Reaction_time=nanmean(Reaction_time_condition);
    
    HR=[False_alarm_condition;nan(5,1);False_alarm_rate];
    RT=[Reaction_time_condition;nan(5,1);Mean_Reaction_time];
    T = table(HR,RT);
    T.Properties.VariableNames = {['False_alarm_target_freq_',num2str(percentage_target_cond(cond)*100)] ['RT_target_freq_',num2str(percentage_target_cond(cond)*100)]};
    Ttotal=T;
    Data_csv_total=[HR RT];
    
    for cond=2:size(Data,1)
        
        False_alarm_condition=Data{cond,4}(:,Subj); % False_alarm in condition
        False_alarm_rate=nanmean(False_alarm_condition);
        Reaction_time_condition=Data{cond,7}(:,Subj); % reaction time in condition
        Mean_Reaction_time=nanmean(Reaction_time_condition);
        HR=[False_alarm_condition;nan(5,1);False_alarm_rate];
        RT=[Reaction_time_condition;nan(5,1);Mean_Reaction_time];
        T = table(HR,RT);
        T.Properties.VariableNames = {['False_alarm_target_freq_',num2str(percentage_target_cond(cond)*100)] ['RT_target_freq_',num2str(percentage_target_cond(cond)*100)]};
        Ttotal=[Ttotal T];
        Data_csv_total=horzcat(Data_csv_total,[HR RT]);
    end
    
    filename = ['MoM_data_Levels_manual_counter_balanced_FA.xlsx']; % Change the name to anything you prefer
    writetable(Ttotal,filename,'Sheet',['Subj_' num2str(Subj)])
    
    %     filename = ['MoM_data_Levels.txt']; % Change the name to anything you prefer
    %     row=1;
    %     col=1;
    %     csvwrite(filename,Data_csv_total,row,col)
    
    [Subj]
end

%% Plotting some results
RT=0; % 1 for reaction time and 0 for hit rate
if RT==0
    dataA1=(1-Data{1,4}(:,subjects))*100;
    dataB1=(1-Data{2,4}(:,subjects))*100;
    dataC1=(1-Data{3,4}(:,subjects))*100;
    dataD1=(1-Data{4,4}(:,subjects))*100;
else
    dataA1=Data{1,7}(:,subjects)*1000;
    dataB1=Data{2,7}(:,subjects)*1000;
    dataC1=Data{3,7}(:,subjects)*1000;
    dataD1=Data{4,7}(:,subjects)*1000;
end

Mean1=nanmean(dataA1,2);
Mean1=Mean1(~isnan(Mean1));

Mean2=nanmean(dataB1,2);
Mean2=Mean2(~isnan(Mean2));

Mean3=nanmean(dataC1,2);
Mean3=Mean3(~isnan(Mean3));

Mean4=nanmean(dataD1,2);
Mean4=Mean4(~isnan(Mean4));

figure;
Shad1=plot([1:length(Mean1)],Mean1,'linewidth',3);
hold on;
Shad2=plot([1:length(Mean2)],Mean2,'linewidth',3);
Shad3=plot([1:length(Mean3)],Mean3,'linewidth',3);
Shad4=plot([1:length(Mean4)],Mean4,'linewidth',3);
xlabel('Block #')
if RT==0
    ylabel({'Percentage of false alarms (%)'})
else
    ylabel('Reaction time (ms)')
end
legend([Shad1,Shad2,Shad3,Shad4],{['Target Freq. = ',num2str(percentage_target_cond(1))],['Target Freq  = ',num2str(percentage_target_cond(2))],['Target Freq  = ',num2str(percentage_target_cond(3))],['Target Freq  = ',num2str(percentage_target_cond(4))]},'location','northwest','edgecolor','none')
