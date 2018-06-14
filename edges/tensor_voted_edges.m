% Author: Emmanuel Maggiori. March 2014.
% 
% Complimentary material for the literature review:
% "Perceptual grouping by tensor voting: a comparative survey of recent approaches". E Maggiori, HL Manterola, M del fresno. To be published in IET Computer Vision.
% 
% Implementation of Steerable Tensor Voting as published in:
% "An efficient method for tensor voting using steerable filters", Franken et. al. ECCV 2006.

function grad_im = tensor_voted_edges(grad_im, sigma)
        s = normd(grad_im, 3);
        be = pi/2 + mod(atan2(grad_im(:,:,2), grad_im(:,:,1)), pi);

        [height, width] = size(s);

        c0=c(0,s,be);
        c2=c(2,s,be);
        c4=c(4,s,be);
        c6=c(6,s,be);
        c2bar=conj(c2);

        w0=w(0,height,width,sigma);
        w2=w(2,height,width,sigma);
        w4=w(4,height,width,sigma);
        w6=w(6,height,width,sigma);
        w8=w(8,height,width,sigma);

        c0_f=fft2(c0);
        c2_f=fft2(c2);	
        c4_f=fft2(c4);		
        c6_f=fft2(c6);
        c2bar_f=fft2(c2bar);

        w0_f=fft2(w0);
        w2_f=fft2(w2);	
        w4_f=fft2(w4);		
        w6_f=fft2(w6);
        w8_f=fft2(w8);

        w0_c2bar=w0_f.*c2bar_f; %eight convolutions required

        w2_c0=w2_f.*c0_f;  

        w4_c2=w4_f.*c2_f;

        w6_c4=w6_f.*c4_f;

        w8_c6=w8_f.*c6_f;

        w0_c0=w0_f.*c0_f;

        w2_c2=w2_f.*c2_f;

        w4_c4=w4_f.*c4_f;


        U_minus2= ifftshift( ifft2( (w0_c2bar) + 4*(w2_c0) + 6*(w4_c2) + 4*(w6_c4) + (w8_c6) ) );
        U_2=conj(U_minus2);
        U_0=  real( ifftshift( ifft2( 6*(w0_c0) + 8*(w2_c2) + 2*(w4_c4) )));

        saliency = abs(U_minus2);
        orientation = 0.5 * angle(U_minus2);

        grad_im = cat(3, saliency .* sin(orientation), saliency .* cos(orientation));
end

function im = c(m,s,be) %s: stickness field, be:orientation field

	im=s.*exp(-1i*m*be);


end

function kernel = w(m,h,w,sigma)

	w2=ceil(w/2.0);
	h2=ceil(h/2.0);



	[x, y] = meshgrid(1:w,1:h);

	kernel = exp(-(((x-w2).^2+(y-h2).^2)/(2*sigma^2))).*(((x-w2)+1i*(y-h2))./(sqrt((x-w2).^2+(y-h2).^2))).^m;

	kernel(h2,w2)=1;

end

        