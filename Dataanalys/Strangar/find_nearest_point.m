function [ Q1, I ] = find_nearest_point( Q0, Q)

[m, I]=min( nansum( ( repmat(Q0, 1,size(Q,2),1) -Q).^2, 3) );

Q1=Q(1, I, :);

end
% 
% function [D_sq] = dist(Q0, Q)
%     tmp=
%     D_sq= nansum( (repmat(Q0, 1, size(Q, 2),1);-Q).^2, 3);
%     %clearvars tmp
% end