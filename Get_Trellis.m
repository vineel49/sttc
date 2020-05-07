% Get trellis 
function [P_State,P_Ip,Ga_Inx]= Get_Trellis()
P_State = [1,2,3,4; 1,2,3,4; 1,2,3,4 ; 1,2,3,4]; % previous state
P_Ip = [1,1,1,1;2,2,2,2;3,3,3,3;4,4,4,4]; % previous ip
Ga_Inx = [1,9,5,13;3,11,7,15;2,10,6,14;4,12,8,16]; % branch indices
end