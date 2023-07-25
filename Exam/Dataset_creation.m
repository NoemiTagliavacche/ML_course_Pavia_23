%program to create dataset for machine learning exam
%I want to generate about 500 matrices, 250 that are physical matrices
%(density matrices) i.e. that represent a quantum state, 250 that are not
%physical matrices

clear all;


%I start from physical matrices (Hilbert space dimension D=4)
Dim=4;


for j=1:300

%First of all I create a random complex matrix of dimension (DimxDim) by
%extracting real and imaginary part from a uniform distribution in the
%interval (-1,1)
    Matrix =(-1+(1+1)*rand(Dim,Dim)) + i*(-1+(1+1)*rand(Dim,Dim));

%Make symmetric
    rho = tril(Matrix,-1)+tril(Matrix,-1)'+diag(diag(Matrix));

%Make diagonal values real
    rho = rho-diag(diag(rho))+real(diag(diag(rho)));

%Make positive definite
    [L, D, P, D0] = modchol_ldlt(rho); 
    rho_chol = P'*L*D*L'*P;

%Make Hermitian
    rho_herm = round(tril(rho_chol,-1)+tril(rho_chol,-1)'+diag(diag(rho_chol)),10);
    Tr=trace(rho_herm);
    
%Normalize trace
    density_matrix = rho_herm./Tr;

    rho_line=reshape(density_matrix,[1,Dim*Dim]);
    if j==1
        memory=rho_line;
    elseif j==2
        matrix_to_save=cat(1,memory,rho_line); %concatenation in the vertical direction
    else
        matrix_to_save=cat(1,matrix_to_save,rho_line);
    end
            
end

writematrix(real(matrix_to_save),'300_Physical_matrices_Dim=4_RealPart.txt','Delimiter','tab')
writematrix(imag(matrix_to_save),'300_Physical_matrices_Dim=4_ImagPart.txt','Delimiter','tab')
clear matrix_to_save;
clear memory;


% Non Physical matrices

%75 random matrices
for k=1:75
    N_physical_mat=(-1+(1+1)*rand(Dim,Dim)) + i*(-1+(1+1)*rand(Dim,Dim));
    N_physical_mat=reshape(N_physical_mat,[1,Dim*Dim]);
    if k==1
        memory=N_physical_mat;
    elseif k==2
        matrix_to_save=cat(1,memory,N_physical_mat); %concatenation in the vertical direction
    else
        matrix_to_save=cat(1,matrix_to_save,N_physical_mat);
    end
end

%75 positive definite matrices
for k=1:75
    Matrix=(-1+(1+1)*rand(Dim,Dim)) + i*(-1+(1+1)*rand(Dim,Dim));
    %Make symmetric
    rho = tril(Matrix,-1)+tril(Matrix,-1)'+diag(diag(Matrix));

    %Make diagonal values real
    rho = rho-diag(diag(rho))+real(diag(diag(rho)));

    %Make positive definite
    [L, D, P, D0] = modchol_ldlt(rho); 
    rho_chol = P'*L*D*L'*P;

    rho_chol=reshape(rho_chol,[1,Dim*Dim]);

    matrix_to_save=cat(1,matrix_to_save,rho_chol);  

end
           


%75 positive definite and hermitian matrices
for k=1:75
    Matrix=(-1+(1+1)*rand(Dim,Dim)) + i*(-1+(1+1)*rand(Dim,Dim));
    %Make symmetric
    rho = tril(Matrix,-1)+tril(Matrix,-1)'+diag(diag(Matrix));

    %Make diagonal values real
    rho = rho-diag(diag(rho))+real(diag(diag(rho)));

    %Make positive definite
    [L, D, P, D0] = modchol_ldlt(rho); 
    rho_chol = P'*L*D*L'*P;

    %Make Hermitian
    rho_herm = round(tril(rho_chol,-1)+tril(rho_chol,-1)'+diag(diag(rho_chol)),10);
    
    rho_herm=reshape(rho_herm,[1,Dim*Dim]);

    matrix_to_save=cat(1,matrix_to_save,rho_herm);  

end

%75 positive definite matrices, trace=1
for k=1:75
    Matrix=(-1+(1+1)*rand(Dim,Dim)) + i*(-1+(1+1)*rand(Dim,Dim));
    %Make symmetric
    rho = tril(Matrix,-1)+tril(Matrix,-1)'+diag(diag(Matrix));

    %Make diagonal values real
    rho = rho-diag(diag(rho))+real(diag(diag(rho)));

    %Make positive definite
    [L, D, P, D0] = modchol_ldlt(rho); 
    rho_chol = P'*L*D*L'*P;
    Tr=trace(rho_chol);
    rho_chol=rho_chol./Tr;

    rho_chol=reshape(rho_chol,[1,Dim*Dim]);

    matrix_to_save=cat(1,matrix_to_save,rho_chol);  

end
writematrix(real(matrix_to_save),'300_Non_Physical_matrices_Dim=4_RealPart.txt','Delimiter','tab')
writematrix(imag(matrix_to_save),'300_Non_Physical_matrices_Dim=4_ImagPart.txt','Delimiter','tab')
