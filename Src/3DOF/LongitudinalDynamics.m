function error = LongitudinalDynamics(predicted_AMatrix_result, simulation_output_states,simulation_error_states, simulation_accutators_result, predicted_BMatrix_result)
% Computes weighted RMSE between linear state-space derivative predictions 
% (x_dot = A*x + B*u) and measured state derivatives.

    x_dot = predicted_AMatrix_result*simulation_output_states + predicted_BMatrix_result*simulation_accutators_result;
    Weights = [0.5 0.2 0.1 0.1 0.1]';
    error = rmse(simulation_error_states,x_dot, Weight=Weights);
end