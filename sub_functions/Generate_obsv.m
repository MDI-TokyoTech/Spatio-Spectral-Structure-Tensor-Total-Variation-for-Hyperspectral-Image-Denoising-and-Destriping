function[HSI_noisy, deg] = Generate_obsv(HSI_clean, deg)
[n1, n2, n3] = size(HSI_clean);

%% Generating stripe noise
stripe_rate             = deg.stripe_rate;
stripe_intensity        = deg.stripe_intensity;

sparse_stripe = 2*(imnoise(0.5*ones(1, n2, n3), "salt & pepper", stripe_rate) - 0.5).* ...
    rand(1, n2, n3).*ones(n1, n2, n3);
stripe_noise = stripe_intensity.*sparse_stripe./max(abs(sparse_stripe), [], "all");

HSI_noisy = HSI_clean + stripe_noise;

deg.stripe_noise = stripe_noise;


%% Generating Gaussian noise
Gaussian_sigma = deg.Gaussian_sigma;

Gaussian_noise = Gaussian_sigma*randn(n1, n2, n3);

HSI_noisy = HSI_noisy + Gaussian_noise;
deg.Gaussian_noise = Gaussian_noise;


%% Generating sparse noise
sparse_rate = deg.sparse_rate;

HSI_tmp = HSI_noisy;

Sp = 0.5*ones(n1, n2, n3);
Sp = imnoise(Sp,'salt & pepper',sparse_rate);

HSI_noisy(Sp==0) = 0;
HSI_noisy(Sp==1) = 1;

deg.sparse_noise = HSI_noisy - HSI_tmp;