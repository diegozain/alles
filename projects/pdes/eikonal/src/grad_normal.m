function grad_u = grad_normal(u)

[nx,nz] = size(u);

% gradient of u
[Dx,Dz] = Dx_Dz(nx,nz);
grad_u = [Dx*u(:), Dz*u(:)];

% normalize
grad_u_ = sqrt( grad_u(:,1).^2 + grad_u(:,2).^2 );
grad_u(:,1) = grad_u(:,1) ./ grad_u_;
grad_u(:,2) = grad_u(:,2) ./ grad_u_;

end