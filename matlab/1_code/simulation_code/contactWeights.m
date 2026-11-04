function w = contactWeights(cPoints, dr)
%CONTACTWEIGHTS Area weights for the pressure under the disk, as the fluid applies it.
%   w(j) is the area of the ring the fluid equations assign to contact node j
%   (2*pi*r_j*dr, and pi*(dr/2)^2 at the centre), in units of the disk radius squared.
%   These are the weights under which the discrete Laplacian is symmetric, so the disk
%   receives exactly the load the water feels. The edge node's pressure includes the
%   surface-tension pull of the meniscus on the rim (its Laplacian row reaches the first
%   free node), so no separate contact-line term is added to the disk force.
    r = (0:cPoints-1) * dr;
    w = 2 * pi * r * dr;
    w(1) = pi * (dr / 2)^2;
end
