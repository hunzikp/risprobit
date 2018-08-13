function LL_evd  = sprobit_llf_hbksg_RR_rho_phi(p,Z,rand_mat,nt,WNT,XNT,TL,i_nt,r)

G = (i_nt-p(1,10)*TL-p(1,11)*WNT);
vcov = Z*(G\(inv(G))')*Z';
ACH = chol(inv(vcov));
BCH = inv(ACH);
GX = G\(p(1,1) + p(1,2)*XNT(:,1) + p(1,3)*XNT(:,2)+ p(1,4)*XNT(:,3)+ p(1,5)*XNT(:,4) + p(1,6)*XNT(:,5) + p(1,7)*XNT(:,6)+ p(1,8)*XNT(:,7) + p(1,9)*XNT(:,8));
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

jnt_prob = sum(ln_prob(31:nt));
LL_evd = -(mean(jnt_prob'));