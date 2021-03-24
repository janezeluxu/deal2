function [variable,Residue] = Newton_solve(meshData,IEN,vertexData,IBC,totalDOF,variable_start)
global Newton_maxIter;
global miu;
%% calculate starting risidual
%% start iteration
variable = variable_start;
p_start = 1;
Newton_tol = [1e-5,1e-5,1e-5,1e-5];
for iter = 1:Newton_maxIter
    [lhs, rhs] = globalKF(iter,variable,totalDOF,meshData,IEN,vertexData,IBC);
    if iter ==1 
        rhs_start = rhs;
    end
    [dvariable] = linear_solve(lhs,rhs);
    variable = variable+dvariable;
    
    %% check convergence
    
    [ev_rhs,ep_rhs,ev,ep] = Newton_check(dvariable,rhs,variable_start,rhs_start,p_start);
    convergence = [ev_rhs,ep_rhs,ev,ep]
    if ev_rhs< Newton_tol(1) && ep_rhs<Newton_tol(2) ...
            && ev<Newton_tol(3) && ep<Newton_tol(4)
        break
    end
    
end
end

function [dvariable] = linear_solve(LHS,RHS)
dvariable = LHS\RHS;
end

function [ev_rhs,ep_rhs,ev,ep] = Newton_check (dvariable,rhs,variable_start,rhs_start,p_start)
global TotalNode;
vel_start= variable_start(1:2*TotalNode);
if isnan(p_start)
    p_start = variable_start(2*TotalNode+1:end);
end
dvel = dvariable(1:2*TotalNode);
dp = dvariable(2*TotalNode+1:end);

vres_start = rhs_start(1:2*TotalNode);
pres_start = rhs_start(2*TotalNode+1:end);

vres = rhs(1:2*TotalNode);
pres = rhs(2*TotalNode+1:end);

ev_rhs = norm(vres)/norm(vres_start);
ep_rhs = norm(pres)/norm(pres_start);

ev = norm(dvel)/norm(vel_start);
ep = norm(dp)/norm(p_start);
end
