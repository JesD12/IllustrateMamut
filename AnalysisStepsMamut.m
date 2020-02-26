function [t_test_scores]=AnalysisStepsMamut(Data)
%% extract the data and calculate intervals
allsteps = Data.jmp_Length_corr;
x = Data.jmp_corr(:,1);
x_mean = mean(x);
x_sd = std(x);
y = Data.jmp_corr(:,2);
y_mean = mean(y);
y_sd = std(y);
z = Data.jmp_corr(:,3);
z_mean = mean(z);
z_sd = std(z);
tracks = unique(Data.TrackID);

% calcualted in the order: 5%, 1%,0.1%,0.01%,0.0001%
boundaries = [0.025, 0.975; 0.005, 0.995; 0.0005, 0.9995;0.00005, 0.99995; 0.0000005, 0.9999995];
x_boundaries = icdf('Normal',boundaries,x_mean,x_sd);
y_boundaries = icdf('Normal',boundaries,y_mean,y_sd);
z_boundaries = icdf('Normal',boundaries,z_mean,z_sd);

%% plot the histogram of the steps
% figure(1);
% clf
% subplot(2,3,1:3);
% histogram(allsteps);
% subplot(2,3,4);
% histogram(x)
% subplot(2,3,5);
% histogram(y)
% subplot(2,3,6);
% histogram(z)

%% calculate t-test-scores
t_test_scores = table('Size',[length(tracks) 5],'VariableTypes',{'string','double','double','double','double'},...
                      'VariableNames',{'track','t_test_x','t_test_y','t_test_z','negative_axes'});
t_test_scores.track=tracks;
xnumbers = zeros(length(tracks),5);
ynumbers = zeros(length(tracks),5);
znumbers = zeros(length(tracks),5);

for k =1:length(tracks)
    
    [t1,t_test_scores.t_test_x(k)] = ttest2(x,Data.jmp_corr(Data.TrackID == tracks(k),1));
    [t2,t_test_scores.t_test_y(k)] = ttest2(y,Data.jmp_corr(Data.TrackID == tracks(k),2));
    [t3,t_test_scores.t_test_z(k)] = ttest2(z,Data.jmp_corr(Data.TrackID == tracks(k),3));
    t_test_scores.anormal_axes(k) = max([t1, t2, t3]);
    t_test_scores.Steps_in_tracks(k) = sum(Data.TrackID == tracks(k));

    for j = 1:5
       xnumbers(k,j)=sum(Data.jmp_corr(Data.TrackID == tracks(k),1) < x_boundaries(j,1) | Data.jmp_corr(Data.TrackID == tracks(k),1) > x_boundaries(j,2));%...
                       %/length(Data.jmp_corr(Data.TrackID == tracks(k),1));
       ynumbers(k,j)=sum(Data.jmp_corr(Data.TrackID == tracks(k),2) < y_boundaries(j,1) | Data.jmp_corr(Data.TrackID == tracks(k),2) > y_boundaries(j,2));%...
                       %/length(Data.jmp_corr(Data.TrackID == tracks(k),2));
       znumbers(k,j)=sum(Data.jmp_corr(Data.TrackID == tracks(k),3) < z_boundaries(j,1) | Data.jmp_corr(Data.TrackID == tracks(k),3) > z_boundaries(j,2));%...
                       %/length(Data.jmp_corr(Data.TrackID == tracks(k),3));
    end


end

t_test_scores.xnumbers = xnumbers;
t_test_scores.ynumbers = ynumbers;
t_test_scores.znumbers = znumbers;

%% filter out any track with negative track ID (as it is a mistake)
%t_test_scores = t_test_scores(t_test_scores.track > 0,:);
