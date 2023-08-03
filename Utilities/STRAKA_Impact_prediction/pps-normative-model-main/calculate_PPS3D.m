function [ output_args ] = calculate_PPS3D( input_args )
% Calculate and plot predictions of 3D model for distances xTs. 
% The figure contains means and 25th/75th percentils.  

%maximal distance for which a prediction is calculated
max_dist=100;
min_dist=0.1;

%experiment parameters (baseline parameters)
vT=[-25, 0, 0];
sigma_v=[20, 5, 5];
sigma_x=[2.5, 5, 5];

body_size=2*[25,25]; %[ysize,zsize]
deltaT=0.5;

FN=5;

%same in all experiments
FP=1;


%number of repetitions of an experiment
number_samples=1000;

%for better performance predicted values are limited to {0, 0.05, 0.1, .., 1} 
impact_decisions=[0:0.05:1];

sigma_displacement=deltaT*sigma_v;


figure
%positions for which predictions are calculated
xTs=[min_dist,5:5:max_dist;0*(0:5:max_dist);0*(0:5:max_dist)]';


%init
pred_impacts_means=-1*ones(size(xTs,1),1);
position_GT_counter=1;


for xT=xTs'
%point estimates
xT_hat=zeros(number_samples,3);
vT_hat=zeros(number_samples,3);
%prevent estimation of the object inside the body
xT_hat(:,1)=max(min_dist, xT(1)+sigma_x(1).*randn(number_samples,1));
% xT_hat(:,1)=xT(1)+sigma_x(1).*randn(number_samples,1);
xT_hat(:,2)=xT(2)+sigma_x(2).*randn(number_samples,1);
xT_hat(:,3)=xT(3)+sigma_x(3).*randn(number_samples,1);

vT_hat(:,1)=vT(1)+sigma_v(1).*randn(number_samples,1);
vT_hat(:,2)=vT(2)+sigma_v(2).*randn(number_samples,1);
vT_hat(:,3)=vT(3)+sigma_v(3).*randn(number_samples,1);

%displacement means
mean_displ=deltaT*vT_hat;

%init
pred_impacts_arr=-1*ones(number_samples,1);


for i=1:number_samples
%future position estimation and Bayesian decision
pred_impact=PPS_prediction_bayes3D(xT_hat(i,:), sigma_x, mean_displ(i,:), sigma_displacement, impact_decisions, FP, FN, body_size);
pred_impacts_arr(i)=pred_impact;
end

pred_impacts_sorted=sort(pred_impacts_arr);
%plot 25th and 75th percentils 
plot([xT(1),xT(1)],[pred_impacts_sorted(round(0.75*number_samples)), pred_impacts_sorted(round(0.25*number_samples))],'-+b', 'MarkerSize',10)
hold on 
pred_impacts_means(position_GT_counter)=mean(pred_impacts_arr);


position_GT_counter=position_GT_counter+1;
end
%plot means
plot((xTs(:,1)),pred_impacts_means,'-ob','LineWidth',2)

xlabel('Distance (in cm)','FontWeight','bold','FontSize',24)
ylabel('y^*_{pred}','FontWeight','bold','FontSize',24)
xlim([-1,max_dist])
ylim([0,1])
ax = gca;
ax.FontSize = 15; 
end