% function [value, isterminal, direction] = MyEvent(t,y)
%     value = [y(1) - 0.001;
%               y(5)-0.3;
%               y(10)-0.55];% scalar event function
%     isterminal = [0;0;0] ;          % 0 = do NOT stop, 1 = stop
%     direction = [+1;+1;+1];           % +1 upward crossing, -1 downward, 0 both
% end
%%
function [value, isterminal, direction] = event_cross_poincare(y,indcomp,val,dir)
    value = y(indcomp) - val;% scalar event function
    isterminal = 0;          % 0 = do NOT stop, 1 = stop
    direction = dir    ;       % +1 upward crossing, -1 downward, 0 both
end