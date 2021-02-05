% Viterbi decoding - soft-input, hard-output
function [dec_ip]=Viterbi_alg(branch_metric,frame_size,decoding_delay)
% Get trellis
[Prev_State,Prev_Ip,Outputs_prev]= Get_Trellis();  
num_states = 4;  % for the considered encoder
% General Initialization
ip=0; % input initialization
dec_ip = zeros(1,frame_size-decoding_delay);
survivor_node = zeros(num_states,frame_size);
survivor_ip = zeros(num_states,frame_size);
path_metric = zeros(num_states,frame_size+1);
P_State_trans = Prev_State'; %transpose of P_State matrix - see trellis.
P_Ip_trans = Prev_Ip'; %transpose of P_Ip matrix - see trellis.
index_temp = [0;1*4;2*4;3*4]; %for ip=4; for linear indexing.
% Algorithm starts
for sym_cnt=  1:frame_size  
[path_metric(:,sym_cnt+1),index] = min([path_metric(Prev_State(:,1),sym_cnt)+ branch_metric(Outputs_prev(:,1),sym_cnt) ...
    path_metric(Prev_State(:,2),sym_cnt)+ branch_metric(Outputs_prev(:,2),sym_cnt) ...
    path_metric(Prev_State(:,3),sym_cnt)+ branch_metric(Outputs_prev(:,3),sym_cnt) ...
    path_metric(Prev_State(:,4),sym_cnt)+ branch_metric(Outputs_prev(:,4),sym_cnt)],[],2);
   
survivor_node(:,sym_cnt) = P_State_trans(index+index_temp);
survivor_ip(:,sym_cnt) = P_Ip_trans(index+index_temp); 
% Back tracing -- 
if (sym_cnt>decoding_delay) 
[~,trace_bf] = min(path_metric(:,sym_cnt+1));
for bk_cnt= 1 : decoding_delay+1
ip = survivor_ip(trace_bf,sym_cnt+1-bk_cnt);
trace_bf = survivor_node(trace_bf,sym_cnt+1-bk_cnt);    
end
dec_ip(sym_cnt-decoding_delay)=ip;
end % for if statement
end % for forward recursiion
dec_ip = dec_ip-1;
end % for function