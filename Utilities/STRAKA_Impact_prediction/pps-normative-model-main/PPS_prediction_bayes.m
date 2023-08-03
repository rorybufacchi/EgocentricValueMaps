function pred_impacts=PPS_prediction_bayes(xT_hat, sigma_x, mean_displacement, sigma_displacement, impact_decisions, FP, FN, r, body_pos)
%PPS_prediction_bayes Summary of this function goes here
%   The function implements future position estimation 
% which is followed by calculation of the Bayesian decision/prediction.
% The prediction is returned.

pred_impacts=[];


pred_impact=predict_impact(xT_hat, sigma_x, mean_displacement, sigma_displacement, impact_decisions, FP, FN, r, body_pos);
pred_impacts=[pred_impacts,pred_impact];



end

function prob=get_impact_prob(xT_hat, sigma_x, mean_displacement, sigma_displacement, body_pos)
prob=cdf('norm',body_pos,xT_hat+mean_displacement,sqrt(sigma_x.^2+sigma_displacement.^2));
end


function loss=calc_loss(impact_state, impact_decision, FP, FN, r)
    loss=FP*(max(0,impact_decision-impact_state).^r)+...
    FN*(max(0,impact_state-impact_decision).^r);
end


function impact=predict_impact(xT_hat, sigma_x, mean_displacement, sigma_displacement, impact_decisions, FP, FN, r, body_pos)
%decision with minimal expected loss is selected
min_val=1000000;
min_dec=-1;
    impact_prob=get_impact_prob(xT_hat, sigma_x, mean_displacement, sigma_displacement, body_pos);
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
