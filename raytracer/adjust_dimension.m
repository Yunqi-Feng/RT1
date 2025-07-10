function [d_max,d_min]=adjust_dimension(d_max,d_min,dim)
if d_max>dim/2 && d_min>dim/2
    d_max=dim;
    d_min=dim/2;
elseif d_max>dim/2 && d_min<dim/2
    d_max=dim*0.75;
    d_min=dim*0.25;
elseif d_max<dim/2 && d_min<dim/2
    d_max=dim/2;
    d_min=0;
end
end