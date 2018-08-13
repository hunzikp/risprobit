function LL_evd = sprobit_llf_tscs_SAM_const(p,Z,rand_mat,n,nt,WNT,XNT,TL,i_nt,r)

G = (i_nt-p(1,3)*WNT-p(1,4)*TL);
vcov = Z*(G\(inv(G))')*Z';
ACH = chol(inv(vcov));
BCH = inv(ACH);
GX = G\(p(1,1) + p(1,2)*XNT);
V = -Z*GX;

for j = 1:r   
       
    nu=zeros(nt,1);
    
        for z = 1:nt
            zz = nt-(z-1); 
            if zz == nt 
                sumterm(zz,j)=0;     
            else 
                sumterm(zz,j)= BCH(zz,:)*nu;    
            end 
            
                nu0(zz,j) = (1/BCH(zz,zz))*(V(zz,1)-sumterm(zz,j));

                
                if nu0(zz,j) < -8
                    nu0(zz,j) = -8;
                else
                end
                if nu0(zz,j) > 8
                    nu0(zz,j) = 8;
                else
                end
                
                
                ln_prob(zz,j) = log(normcdf(nu0(zz,j)));
                nu(zz,1) = norminv(rand_mat(zz,j)*(normcdf(nu0(zz,j))));     
                
                if nu(zz,1) < -8
                    nu(zz,1) = -8;
                else
                end
                if nu(zz,1) > 8
                    nu(zz,1) = 8;
                else
                end 
      
        end
        
end 

jnt_prob = sum(ln_prob(n+1:nt));
LL_evd = -(mean(jnt_prob'));
