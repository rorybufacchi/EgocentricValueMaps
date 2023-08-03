function [ output_args ] = FN_std_PPS_properties( input_args )
%FN_STD_PPS_PROPERTIES Summary of this function goes here
%   Plot S1 Fig 
max_dist=100;
deltaT=0.5;

vT=-25;
sigma_v_arr=[2.5:2.5:30];


FP=1;
FNs=[1,2,5,10,17.5,25,50,75,100];


r=2;

sigma_x=0.0;

xTs=0:1:max_dist;

beginning_threshold=0.01;

number_samples=1000;
beginnings=-1*ones(numel(FNs),numel(sigma_v_arr));
slopes=-1*ones(numel(FNs),numel(sigma_v_arr));

for i=1:numel(FNs)
    FN=FNs(i)
    for j=1:numel(sigma_v_arr)
       sigma_v = sigma_v_arr(j);
       pred_tact_acts_means=get_PPS_impact_means_distance(vT, sigma_v, FP, FN, r, number_samples, deltaT, sigma_x, xTs);
       beginning_idx = find_PPS_beginning(pred_tact_acts_means,beginning_threshold);
       beginnings(i,j)=xTs(beginning_idx);
       slopes(i,j)=calculate_slope(xTs,pred_tact_acts_means);
    end  
    
end
figure
surf(sigma_v_arr,FNs,beginnings)
colorbar
xlabel('\sigma_v (in cm/s)')
ylabel('FN')
zlabel('beginning (in cm)')
figure
surf(sigma_v_arr,FNs,-slopes)
colorbar
xlabel('\sigma_v (in cm/s)')
ylabel('FN')
zlabel('slope (in cm^{-1})')
end

function beginning_idx = find_PPS_beginning(pred_tact_acts_means,beginning_threshold)
for i=numel(pred_tact_acts_means):-1:1
    if pred_tact_acts_means(i)>beginning_threshold
        beginning_idx=i;
        return
    end
end
end

function slope=calculate_slope(xTs,pred_impact_means)
    ymax=max(pred_impact_means);
    ymin=min(pred_impact_means);
    center=(ymax-ymin)/2;
    interv=0.15;
    y1=center+interv;
    y2=center-interv;

    for i=1:numel(pred_impact_means)-1
        if pred_impact_means(i)>y1 & pred_impact_means(i+1)<=y1
            %line between the two points
        [P1,S] = polyfit([pred_impact_means(i),pred_impact_means(i+1)],[xTs(i),xTs(i+1)],1);
        x1=P1(1)*y1+P1(2);
        end
        if pred_impact_means(i)>y2 & pred_impact_means(i+1)<=y2
        [P2,S] = polyfit([pred_impact_means(i),pred_impact_means(i+1)],[xTs(i),xTs(i+1)],1);
        x2=P2(1)*y2+P2(2);
        end
    end
        slope=(y2-y1)/(x2-x1);
        return 
end

function pred_impact_means=get_PPS_impact_means_distance(vT, sigma_v, FP, FN, r, number_samples, deltaT, sigma_x, xTs)
%get_PPS_bayes_means_distance Summary of this function goes here
%  for each distance it returns mean predicted impact
body_pos=0;

impact_decisions=[0:0.05:1];

sigma_displacement=deltaT*sigma_v;

%inits
pred_impact_means=-1*ones(numel(xTs),1);
position_GT_counter=1;


for xT=xTs
%xT_hat=xT+sigma_x.*randn(number_samples,1);
xT_hat=max(0.1,xT+sigma_x.*randn(number_samples,1));
vT_hat=vT+sigma_v.*randn(number_samples,1);

mean_displ=deltaT.*vT_hat;


pred_impacts_arr=-1*ones(number_samples,1);

for i=1:number_samples
pred_impact=PPS_prediction_bayes(xT_hat(i), sigma_x, mean_displ(i), sigma_displacement, impact_decisions, FP, FN, r, body_pos);
pred_impacts_arr(i)=pred_impact;
end

pred_impact_means(position_GT_counter)=mean(pred_impacts_arr);

position_GT_counter=position_GT_counter+1;
end

end


