
function [value, isterminal, direction] = myEvent3DOF_LMO(~, state_vector, ~)

%         vehicleObjTemp.E  = state_vector(1);  
%         vehicleObjTemp.F  = state_vector(2);  
%         vehicleObjTemp.G  = state_vector(3);
%         vehicleObjTemp.dE = state_vector(4); 
%         vehicleObjTemp.dF = state_vector(5);  
%         vehicleObjTemp.dG = state_vector(6);
        
        %[vehicleObjTemp] = ConvertVehicleStructureTrial_Mars(vehicleObjTemp, 'edm', 'km');    

        Rnorm = norm([state_vector(1); state_vector(2); state_vector(3)]);

        % Event we are looking for %Once you hit zero and then go negative, terminate the propagation...
        value       = Rnorm - Mars_GenPhysCons.RE_EQ; 
        
        isterminal = 1;  %if we want the integration to stop when this event it's 1.. . 
        direction  = -1;   

end
