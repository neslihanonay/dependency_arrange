%%This code is modified version used in Bisby and colleagues study (2018).
Original source is here; https://figshare.com/articles/software/Model_code/5240752?file=8956867


T_D = readtable('subj109_disney.xlsx', 'Range', 'B1:G41') % range covers all data(including both emotional and neutral)
T_mat_D=table2array(T_D);
T_N_double = str2double (T_mat_D);
res_D = T_N_double


pairs = [[1 2] ;[3 4] ;[5 6]; [1 6] ;[2 4]; [3 5]]

dependency_all_pairs_D = zeros(6,3)
for ii = 1:size(pairs,1)
    pair = pairs(ii,:)
    
%function [dep] = dependency(res,pair,guess,c)

%%function to calculate behavioural dependency measure
% Aidan J Horner 04/2015

% input:
% res       =   MxN matrix of M 'events' and N 'retrieval trials'
% pair      =   pair of retrieval trials to calculate dependency 
%               (e.g., cue location retrieve object and cue location retrieve person)
%               example = [1 2] - uses first two columns to build contigency table for analysis
% optional inputs:
% guess     =   include guessing in Dependent Model (1 or 0)
%               default = 1
% c         =   number of choices (i.e., c-alternative forced choice) - for
%               estimating level of guessing
%               default = 6

% output:
% dep       =   dependency measure for [data independent_model dependent_model]

%% housekeeping

%if nargin   < 3
    guess   = 1;                                            % set 'guess' to 1 if not defined by user
    c       = 4;                                            % set 'c' to 6 if not defined by user
%end

%% calculate dependency for data
 
res2        = res_D(:,pair);                                  % create column Mx2 matrix for retrieval trials defined by 'pair'
dep(1)      = sum(sum(res2,2)~=1)/size(res2,1)       % calculate dependency for data

%% calculate dependency for independent model

acc         = mean(res2,1);                                 % calculate accuracy for each retrieval type
dep(2)      = ((acc(1)*acc(2))+((1-acc(1))*(1-acc(2))));    % calculate dependenct for independent model

%% calculate dependency for dependent model

cont        = nan(size(res2,1),2,2);                        % create matrix for dependent model probabilities
g           = (1-mean(res_D(:)))*(c/(c-1));                   % calculate level of guessing
b           = mean(res_D); %b(:,pair) = nan;                   % calculate average performance
for i       = 1:size(res2,1)                                % loop through all event   
    a       = res_D(i,:); %a(:,pair) = nan;                    % calculate event specific performance
    E       = mean(a)/mean(b);                        % calculate ratio of event / average performance (episodic factor)
    for p = 1:2;
        if E*acc(p)>1
            P(p) = 1;
        else
            if guess == 1
                P(p) = (E*(acc(p)-(g/c)))+(g/c);
            elseif guess == 0
                P(p) = E*acc(p);
            end
        end
    end
    cont(i,1,1) = P(1)*P(2);
    cont(i,1,2) = (1-P(1))*P(2);
    cont(i,2,1) = P(1)*(1-P(2));
    cont(i,2,2) = (1-P(1))*(1-P(2));
end
cont2       = squeeze(sum(cont));                           % create contingency table
if size(cont2) == [2,1]
    dep(3)      = cont2(1,1)/sum(cont2(:));
else
    dep(3)      = (cont2(1,1)+cont2(2,2))/sum(cont2(:))        % calculate dependency for dependent model
end

dependency_all_pairs_D(ii,:) = dep(1,:)

end 

T_N = readtable('subj109_disney.xlsx', 'Range', 'B1:G21') % range covers neutral data
T_mat_N = table2array(T_N);
T_N_double = str2double (T_mat_N);
res_N = T_N_double

pairs = [[1 2] ;[3 4] ;[5 6]; [1 6] ;[2 4]; [3 5]]

dependency_all_pairs_N = zeros(6,3)
for ii = 1:size(pairs,1)
    pair = pairs(ii,:)
    
%function [dep] = dependency(res,pair,guess,c)

%%function to calculate behavioural dependency measure
% Aidan J Horner 04/2015

% input:
% res       =   MxN matrix of M 'events' and N 'retrieval trials'
% pair      =   pair of retrieval trials to calculate dependency 
%               (e.g., cue location retrieve object and cue location retrieve person)
%               example = [1 2] - uses first two columns to build contigency table for analysis
% optional inputs:
% guess     =   include guessing in Dependent Model (1 or 0)
%               default = 1
% c         =   number of choices (i.e., c-alternative forced choice) - for
%               estimating level of guessing
%               default = 6

% output:
% dep       =   dependency measure for [data independent_model dependent_model]

%% housekeeping

%if nargin   < 3
    guess   = 1;                                            % set 'guess' to 1 if not defined by user
    c       = 4;                                            % set 'c' to 6 if not defined by user
%end

%% calculate dependency for data
 
res2        = res_N(:,pair);                                  % create column Mx2 matrix for retrieval trials defined by 'pair'
dep(1)      = sum(sum(res2,2)~=1)/size(res2,1)       % calculate dependency for data

%% calculate dependency for independent model

acc         = mean(res2,1);                                 % calculate accuracy for each retrieval type
dep(2)      = ((acc(1)*acc(2))+((1-acc(1))*(1-acc(2))));    % calculate dependenct for independent model

%% calculate dependency for dependent model

cont        = nan(size(res2,1),2,2);                        % create matrix for dependent model probabilities
g           = (1-mean(res_N(:)))*(c/(c-1));                   % calculate level of guessing
b           = mean(res_N); %b(:,pair) = nan;                   % calculate average performance
for i       = 1:size(res2,1)                                % loop through all event   
    a       = res_N(i,:); %a(:,pair) = nan;                    % calculate event specific performance
    E       = mean(a)/mean(b);                        % calculate ratio of event / average performance (episodic factor)
    for p = 1:2;
        if E*acc(p)>1
            P(p) = 1;
        else
            if guess == 1
                P(p) = (E*(acc(p)-(g/c)))+(g/c);
            elseif guess == 0
                P(p) = E*acc(p);
            end
        end
    end
    cont(i,1,1) = P(1)*P(2);
    cont(i,1,2) = (1-P(1))*P(2);
    cont(i,2,1) = P(1)*(1-P(2));
    cont(i,2,2) = (1-P(1))*(1-P(2));
end
cont2       = squeeze(sum(cont));                           % create contingency table
if size(cont2) == [2,1]
    dep(3)      = cont2(1,1)/sum(cont2(:));
else
    dep(3)      = (cont2(1,1)+cont2(2,2))/sum(cont2(:))        % calculate dependency for dependent model
end

dependency_all_pairs_N(ii,:) = dep(1,:)

end 

T_Em = readtable('subj109_disney.xlsx', 'Range', 'B21:G41') % range covers emotional data
T_mat_Em = table2array(T_Em); 
T_Em_double = str2double (T_mat_Em);
res_Em = T_Em_double

pairs = [[1 2] ;[3 4] ;[5 6]; [1 6] ;[2 4]; [3 5]]

dependency_all_pairs_Em = zeros(6,3)
for ii = 1:size(pairs,1)
    pair = pairs(ii,:)
    
%function [dep] = dependency(res,pair,guess,c)

%%function to calculate behavioural dependency measure
% Aidan J Horner 04/2015

% input:
% res       =   MxN matrix of M 'events' and N 'retrieval trials'
% pair      =   pair of retrieval trials to calculate dependency 
%               (e.g., cue location retrieve object and cue location retrieve person)
%               example = [1 2] - uses first two columns to build contigency table for analysis
% optional inputs:
% guess     =   include guessing in Dependent Model (1 or 0)
%               default = 1
% c         =   number of choices (i.e., c-alternative forced choice) - for
%               estimating level of guessing
%               default = 6

% output:
% dep       =   dependency measure for [data independent_model dependent_model]

%% housekeeping

%if nargin   < 3
    guess   = 1;                                            % set 'guess' to 1 if not defined by user
    c       = 4;                                            % set 'c' to 6 if not defined by user
%end

%% calculate dependency for data
 
res2        = res_Em(:,pair);                                  % create column Mx2 matrix for retrieval trials defined by 'pair'
dep(1)      = sum(sum(res2,2)~=1)/size(res2,1)       % calculate dependency for data

%% calculate dependency for independent model

acc         = mean(res2,1);                                 % calculate accuracy for each retrieval type
dep(2)      = ((acc(1)*acc(2))+((1-acc(1))*(1-acc(2))));    % calculate dependenct for independent model

%% calculate dependency for dependent model

cont        = nan(size(res2,1),2,2);                        % create matrix for dependent model probabilities
g           = (1-mean(res_Em(:)))*(c/(c-1));                   % calculate level of guessing
b           = mean(res_Em); %b(:,pair) = nan;                   % calculate average performance
for i       = 1:size(res2,1)                                % loop through all event   
    a       = res_Em(i,:); %a(:,pair) = nan;                    % calculate event specific performance
    E       = mean(a)/mean(b);                        % calculate ratio of event / average performance (episodic factor)
    for p = 1:2;
        if E*acc(p)>1
            P(p) = 1;
        else
            if guess == 1
                P(p) = (E*(acc(p)-(g/c)))+(g/c);
            elseif guess == 0
                P(p) = E*acc(p);
            end
        end
    end
    cont(i,1,1) = P(1)*P(2);
    cont(i,1,2) = (1-P(1))*P(2);
    cont(i,2,1) = P(1)*(1-P(2));
    cont(i,2,2) = (1-P(1))*(1-P(2));
end
cont2       = squeeze(sum(cont));                           % create contingency table
if size(cont2) == [2,1]
    dep(3)      = cont2(1,1)/sum(cont2(:));
else
    dep(3)      = (cont2(1,1)+cont2(2,2))/sum(cont2(:))        % calculate dependency for dependent model
end

dependency_all_pairs_Em(ii,:) = dep(1,:)

M_D = mean(dependency_all_pairs_D);
M_Em = mean(dependency_all_pairs_Em);
M_N = mean(dependency_all_pairs_N);
double1 = [M_D;M_Em;M_N]

rowNames = {'all_data','emotional','neutral'};
colNames = {'data_dependency','independent_dependency','dependent_dependency'};
subj109_model = array2table(double1,'RowNames',rowNames,'VariableNames',colNames)
writetable(subj109_model,'subj109_model_disney.xlsx')

end
 
