function [T0,T1,T2] = coupleCartesianSpherically(vec1,vec2)
        % Couples two Cartesian vectors and
        % returns spherical tensors of rank 0, 1 and 2
        % see Brown and Carrington Eq. 5.113-5.118
        ux = vec1(1);uy = vec1(2);uz = vec1(3);
        vx = vec2(1);vy = vec2(2);vz = vec2(3);

        T0 = cell(1,1); %rank 0 tensor -> {T0_0}
        T1 = cell(3,1); %rank 1 tensor -> {T1_-, T1_0, T1_+}
        T2 = cell(5,1); %rank 2 tensor -> {T2_-2, T2_-1, T2_0, T2_+1, T2_+2}

        T0{1} = -(1/sqrt(3))*(ux*vx + uy*vy + uz*vz); %normal dot product but smaller
        
        T1{1} = (1/2)*((ux*vz - uz*vx) + 1i*(uy*vz - uz*vy));
        T1{2} = (1/sqrt(2))*(ux*vy - uy*vx);
        T1{3} = (1/2)*((ux*vz - uz*vx) - 1i*(uy*vz - uz*vy));

        T2{1} = (1/2)*(ux*vx - uy*vy - 1i*(ux*vy + uy*vx)); %T2_-2
        T2{2} = (1/2)*(ux*vz + uz*vx - 1i*(uy*vz + uz*vy)); %T2_-1
        T2{3} = (1/sqrt(6))*(2*uz*vz - ux*vx - uy*vy);      %T2_0
        T2{4} = -(1/2)*(ux*vz + uz*vx + 1i*(uy*vz + uz*vy));%T2_1
        T2{5} = (1/2)*(ux*vx - uy*vy + 1i*(ux*vy + uy*vx)); %T2_2
end

