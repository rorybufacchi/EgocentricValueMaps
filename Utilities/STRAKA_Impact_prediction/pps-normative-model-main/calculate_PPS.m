function [ output_args ] = calculate_PPS( input_args )
%calculate_PPS Summary of this function goes here
%   Calculate and plot predictions for distances xTs. 
% The figure contains means and 25th/75th percentils.  

%maximal distance for which a prediction is calculated
max_dist=100;
body_pos=0;

%experiment parameters (baseline parameters)
vT=-25;
sigma_v=20;
sigma_x=2.5;

FN=5;

%same in all experiments
FP=1;
r=2;
deltaT=0.5;

%number of repetitions of an experiment
number_samples=1000;

%for better performance predicted values are limited to {0, 0.05, 0.1, .., 1} 
impact_decisions=[0:0.05:1];

sigma_displacement=deltaT*sigma_v;


figure
%positions for which predictions are calculated
xTs=0:5:max_dist;


%init
pred_impacts_means=-1*ones(numel(xTs),1);
position_GT_counter=1;


for xT=xTs
%point estimates
%prevent estimation of the object inside the body
xT_hat=max(0.1,xT+sigma_x.*randn(number_samples,1));
vT_hat=vT+sigma_v.*randn(number_samples,1);

%displacement means
mean_displ=deltaT.*vT_hat;

%init
pred_impacts_arr=-1*ones(number_samples,1);


for i=1:number_samples
%future position estimation and Bayesian decision
pred_impact=PPS_prediction_bayes(xT_hat(i), sigma_x, mean_displ(i), sigma_displacement, impact_decisions, FP, FN, r, body_pos);
pred_impacts_arr(i)=pred_impact;
end

pred_impacts_sorted=sort(pred_impacts_arr);
%plot 25th and 75th percentils 
plot([xT,xT],[pred_impacts_sorted(750), pred_impacts_sorted(250)],'-+b', 'MarkerSize',10)
hold on 
pred_impacts_means(position_GT_counter)=mean(pred_impacts_arr);


position_GT_counter=position_GT_counter+1;
end
%plot means
plot((xTs),pred_impacts_means,'-ob','LineWidth',2)

xlabel('Distance (in cm)','FontWeight','bold','FontSize',24)
ylabel('y^*_{pred}','FontWeight','bold','FontSize',24)
xlim([-1,max_dist])
ylim([0,1])
ax = gca;
ax.FontSize = 15; 
end