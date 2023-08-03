# pps-normative-model
A normative model of peripersonal space as impact prediction. The model is implemented in Matlab.

## calculate_PPS.m
Run this code to get a figure with dependency between the mean of 1000 predicted
tactile activations (for each distance) and distance xT (in centimeters) of the
stimuli from the body. Moreover, 25th and 75th percentils are plotted. Feel free to change the parameters of the experiment and compare the PPS profile with empirical data.

## calculate_PPS3D.m
A 3D version of **calculate_PPS.m**.

## PPS_prediction_bayes.m
This function implements future position estimation which is followed by calculation of the Bayesian decision/prediction.

## PPS_prediction_bayes3D.m
A 3D version of **PPS_prediction_bayes.m**.

## FN_std_PPS_properties.m
Plot S8 Fig (for one prediction time step size).

