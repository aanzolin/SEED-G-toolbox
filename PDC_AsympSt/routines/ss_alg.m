function SS = ss_alg(A, e_cov, nf)
%    '''Calculates the Spectral density (SS)
%         A -> autoregressive matrix
%         e_cov -> residues
%         nf -> number of frequencies
%         '''
    [n, n, r] = size(A);
    AL = A_to_f(A, nf);
    ss = zeros(size(AL));
    for i = 1:nf,
        H = inv(reshape(AL(i,:,:),n,n));
        ss(i,:,:) = H*e_cov*H';
    end;    
    %print ss[5]
    SS=permute(ss,[2,3,1]);
%    return ss.transpose(1,2,0)