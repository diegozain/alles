function [p,pp,EE,ss] = gd_dc(p,s_dc,Mdc,d_dc_o,Ndc,max_iter)
% ..............................................................................
pp = p;
E_E = 1;
count_ = 0;
EE = [];
ss = [];
% ..............................................................................
while E_E > 1e-25
  % forward run & objective calculator
  %
  [u_dc,e_dc,Ldc] = fwd_dc(p,s_dc,Mdc,d_dc_o);
  E_E = Edc_obj(e_dc,Ndc);
  % ............................................................................
  % calculate & store gradient 
  % (- descent direction)
  %
  g_dc = Jte_dc(u_dc,s_dc,p,e_dc,Mdc,Ldc);
  g_dc_norm = norm(g_dc);
  % ............................................................................
  % calculate step size
  %
  step_dc = ( g_dc_norm / norm( Jg_dc(u_dc,s_dc,p,g_dc,Mdc,Ldc) ) )^2;
  % ............................................................................
  % update parameter
  %
  p  = p - step_dc*g_dc;
  % ............................................................................
  % store for plotting
  %
  pp = [pp p];
  EE = [EE E_E];
  ss = [ss step_dc];
  % ............................................................................
  count_ = count_ + 1;
  if count_ > max_iter
      break
  end
end
% ..............................................................................
end
