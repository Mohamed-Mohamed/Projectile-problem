function [ e, Max_e, std_e, mean_e, RMS_e, e_abs, Max_abs_e, std_abs_e, mean_abs_e, RMS_abs_e ] = ERROR ( A1,A2 )
% This function is used to get error of A1 relative to A2
%% Coded by
% Mohamed Mohamed El-Sayed Atyya
% mohamed.atyya94@eng-st.cu.edu.eg
%% inputs:
% A1        : first solution
% A2        : second solution
% NOTE : A1 AND A2 MUST BE IN THE SAME SIZE 
%% outputs:
% e                         : Error
% Max_e                : Max. Error
% std_e                  : standard diviation of error
% mean_e              : mean of error
% RMS_e                : root mean squre of error 
% e_abs                 : abs. Error
% Max_abs_e         : Max. abs. Error
% std_abs_e          : standard diviation of abs. error
% mean_abs_e      : mean of abs. error
% RMS_abs_e        : root mean squre of abs. error 
%% Function body
e=A1-A2;
Max_e=max(e);
std_e=std(e);
mean_e=mean(e);
RMS_e=rms(e);
% -----------------------------------------------------------------------------------------------------------------------------------------------------------
e_abs=abs(A1-A2);
Max_abs_e=max(abs(e_abs));
std_abs_e=std(e_abs);
mean_abs_e=mean(e_abs);
RMS_abs_e=rms(e_abs);
end