% Get trellis 
function [Prev_State,Prev_Ip,Outputs_prev]= Get_Trellis()
Prev_State = [1,2,3,4; 1,2,3,4; 1,2,3,4 ; 1,2,3,4]; % previous state
Prev_Ip = [1,1,1,1;2,2,2,2;3,3,3,3;4,4,4,4]; % previous ip
Outputs_prev = [1,9,5,13;3,11,7,15;2,10,6,14;4,12,8,16]; % branch indices
end