function surrogate=gen_surrogate(DATA,k)

% % % GENERATION OF phase incoherent SURROGATE DATA

fdata=fft(DATA);
mag=abs(fdata);
pdata=angle(fdata);

    rand('state',k);
    [m,n]=size(pdata);
    r=rand(m,n);
    mn=mean(mean(r));
    r=r-mn;
    pdata=r*pi; % generation of random phase values between -pi and pi. 
    fdata=mag.*exp(i*pdata);% phase of the data is replaced by random 
                           % values to remove any phase coherence.
    srdata=ifft(fdata);
    surrogate(:,:)=real(srdata); % this is a surrogate data with random phase.
 %    surrogate(k,:,l)=sur_gs1(surDATA); % this is a surrogate data 
 %  end
