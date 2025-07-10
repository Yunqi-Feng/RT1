function [z_trsh] = adjust_trsh(scatterer_height_max,z_trsh,confg_safety,int_safety)
scatterer_safety=confg_safety+int_safety;
if scatterer_height_max<0.5*scatterer_safety(3)+0.25*z_trsh
    error('Please increase the scatterer_height_max')
elseif scatterer_height_max<0.5*scatterer_safety(3)+0.5*z_trsh && scatterer_height_max>0.5*scatterer_safety(3)+0.25*z_trsh
    z_trsh=2*(scatterer_height_max-0.5*scatterer_safety(3));
end
end