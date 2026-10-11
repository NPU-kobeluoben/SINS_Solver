function gn = egm(GM, Re, wie, f, C, S, lat, lon, hgt)
    N = length(C) - 1; % order
    C = [C(:,1) * sqrt(2), C(:,2:end) * 2]; 
    S = [S(:,1) * sqrt(2), S(:,2:end) * 2]; 
    
    slat = sin(lat); clat = cos(lat); 
    slon = sin(lon); clon = cos(lon);
    e2 = 2 * f - f^2; 
    RN = Re / sqrt(1 - e2 * slat^2);
    
    x = (RN + hgt) * clat * clon; 
    y = (RN + hgt) * clat * slon; 
    z = (RN * (1 - e2) + hgt) * slat;
    
    theta = atan2(sqrt(x^2 + y^2),z); 
    phi = pi/2-theta;
    r = sqrt(x^2 + y^2 + z^2); 
    dL = lat - phi;
    sth = sin(theta); 
    cth = cos(theta);
    
    Pnm = zeros(1,N+1); 
    Pn_2m = Pnm; 
    Pn_1m = Pnm; 
    dPnm = Pnm;
    P = zeros(N+1); 
    dP = P;
    
    for n = 0:N % Legendre
        n1 = n + 1;
        if n == 0
            Pnm(n1) = 1/sqrt(2);
            dPnm(n1) = 0;
        else
            Pnm(n1) = sqrt((2*n+1)/(2*n)) * Pn_1m(n1-1) * sth;
            Pnm(n1-1) = sqrt(2*n+1) * cth * Pn_1m(n1-1);
            if n >= 2
                m = 0:n-2; 
                k = 1:n1-2;
                Pnm(k) = sqrt((2*n+1)./(n^2-m.^2)) .* ...
                    (sqrt(2*n-1) * cth * Pn_1m(k) - ...
                    sqrt(((n-1)^2-m.^2)./(2*n-3)) .* Pn_2m(k));
            end
            dPnm(n1) = -n * cth/sth^2 * Pnm(n1);
            m = 0:n-1; 
            k = 1:n1-1;
            dPnm(k) = (-n * cth * Pnm(k) + ...
                sqrt((2*n+1)./(2*n-1)) .* sqrt(n^2-m.^2) .* Pn_1m(k))./sth^2;
        end
        P(n1, :) = Pnm;
        dP(n1, :) = dPnm;
        Pn_2m = Pn_1m; 
        Pn_1m = Pnm;
    end
    
    dP = -sth.*dP;
    
    N0 = (0:N);
    Rr = (Re/r).^N0';

    clon = cos(N0 * lon); 
    slon = sin(N0 * lon);
 
    CP = C .* P; 
    SP = S .* P; 
    
    xx = Rr .* N0 .*(- slon .* CP+ clon.* SP);
    yy = Rr .*(clon .* C + slon .* S) .* dP;
    zz = (N0 + 1)' .* Rr .* (clon .* CP+slon .* SP);
    
    gE0 = GM/r^2/sth * sum(sum(xx));
    gN0 = -GM/r^2 * sum(sum(yy)) - wie^2/2 * r * sin(2 * theta);
    gU0 = -GM/r^2 * sum(sum(zz)) + wie^2 * r * sth^2;
    
    gn = [gE0, gN0 * cos(dL) - gU0 * sin(dL), gN0 * sin(dL) + gU0 * cos(dL)]';
end