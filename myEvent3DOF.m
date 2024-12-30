
function [value, isterminal, direction] = myEvent3DOF(~, state_vector, ~)

        Rnorm = norm([state_vector(1); state_vector(2); state_vector(3)]);

        % Event we are looking for %Once you hit zero and then go negative, terminate the propagation...
            value       = Rnorm - Mars_GenPhysCons.RE_EQ; 
            
            isterminal = 1;  %if we want the integration to stop when this event it's 1.. . 
            direction  =-1;   

end
