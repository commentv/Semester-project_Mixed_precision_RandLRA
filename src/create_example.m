function A = create_example(example,n,q,G,ksi,p)
%%%Function that creates the matrices for the examples of the figures in Section 3%%%

%%% Parameters %%%
%Example : to choose the example
%n : size of the matrix
%q : decay of the exponetial decay
%G : random gaussian matrix
%ksi : value of the noise parameter
%p : decay of the polynomial decay

switch example
    case 'expdecay'
        val = 10.^(-q * (0:(n-1)) );
        A = diag(val);

    case 'psdNoise'
        A = ksi*(1/n)*(G*G');

    case 'poldecay'
        val = (1:n).^(-p);
        A = diag(val);
        
    case 'stairdecay'
        val = [1,0.99,0.98];
        remainder = rem(n,3);
        vect = val(1:remainder);
        iter = floor(n/3);
        for i = 1:iter
            vect = [val vect/10];
        end
        
        A=diag(vect);

end
end