function pred_impacts=PPS_prediction_bayes3D(xT_hat, sigma_x, mean_displacement, sigma_displacement, impact_decisions, FP, FN, body_size)
% The function implements 3D future position estimation, hit probability calculation 
% which is followed by calculation of the Bayesian decision/prediction.
% The prediction is returned.

pred_impacts=[];


pred_impact=predict_impact(xT_hat, sigma_x, mean_displacement, sigma_displacement, impact_decisions, FP, FN, body_size);
pred_impacts=[pred_impacts,pred_impact];



end

function prob=get_impact_prob(xT_hat, sigma_x, mean_displacement, sigma_displacement, body_size)
%3D future position estimation
new_position_mean=xT_hat+mean_displacement;
%sigma for each position coordinate
new_position_sigma=sqrt(sigma_x.^2+sigma_displacement.^2);

%numerical solution - sampling
%generate samples for the numerical estimation of the probability
num_samples=10000;
new_position_gen_samples=zeros(num_samples,3);
new_position_gen_samples(:,1) = new_position_mean(1)+new_position_sigma(1).*randn(num_samples,1);
new_position_gen_samples(:,2) = new_position_mean(2)+new_position_sigma(2).*randn(num_samples,1);
new_position_gen_samples(:,3) = new_position_mean(3)+new_position_sigma(3).*randn(num_samples,1);
positive=0;
negative=0;

for i=1:num_samples
    new_pos_sample_x=new_position_gen_samples(i,1);
    new_pos_sample_y=new_position_gen_samples(i,2);
    new_pos_sample_z=new_position_gen_samples(i,3);
    %the samples is in front of the body
    if new_pos_sample_x>0
    negative=negative+1;
    continue
    end
    
    if ((xT_hat(2) + 0.5*body_size(1))/xT_hat(1))*new_pos_sample_x-0.5*body_size(1) <= new_pos_sample_y ...
            && new_pos_sample_y <= ((xT_hat(2) - 0.5*body_size(1))/xT_hat(1))*new_pos_sample_x+0.5*body_size(1) ...
       && ((xT_hat(3) + 0.5*body_size(2))/xT_hat(1))*new_pos_sample_x-0.5*body_size(2) <= new_pos_sample_z ...
            && new_pos_sample_z <= ((xT_hat(3) - 0.5*body_size(2))/xT_hat(1))*new_pos_sample_x+0.5*body_size(2) 
        positive=positive+1;
    else
        negative=negative+1;
    end
end

prob=positive/(positive+negative);



end


function loss=calc_loss(impact_state, impact_decision, FP, FN, r)
    loss=FP*(max(0,impact_decision-impact_state).^r)+...
    FN*(max(0,impact_state-impact_decision).^r);
end


function impact=predict_impact(xT_hat, sigma_x, mean_displacement, sigma_displacement, impact_decisions, FP, FN, body_size)
%decision with minimal expected loss is selected
r=2;
min_val=1000000;
min_dec=-1;
    impact_prob=get_impact_prob(xT_hat, sigma_x, mean_displacement, sigma_displacement, body_size);
    hit_state=1;
    for impact_decision=impact_decisions    
    val=impact_prob*calc_loss(hit_state, impact_decision, FP, FN, r)...
        +(1-impact_prob)*calc_loss(hit_state-1, impact_decision, FP, FN, r);
    if val<min_val
        min_val=val;
        min_dec=impact_decision;
    end
    end
impact=min_dec;
end
