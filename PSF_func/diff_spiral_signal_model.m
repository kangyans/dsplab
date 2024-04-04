function seq = diff_spiral_signal_model(seq)

        
    seq = seq_time_design(seq);

    if strcmp(seq.decaymodel,'bloch') == 1
        seq.decay = seq_sim_decay_generator(seq);
    end

    if strcmp(seq.decaymodel,'epg') == 1
        seq.decay = seq_sim_decay_epg_generator(seq);
    end
    
    % simulate k-space signal.
    
    FT      = NUFFT(seq.k,1,seq.imsize);
    kdata   = FT*seq.image;
    
    if seq.offres_mode == 1
        seq.offres_tmap = repmat(seq.tADC,1, size(seq.k,2))*1e-3; % convert to second
        A     = CPR(seq.offres_fmap, seq.offres_tmap, seq.k, seq.offres_method, seq.offres_nseg);
        kdata = A*seq.image;
    else 
        A = FT;
    end


    for i = 1:seq.Nint
        seq.kdata(:,i) = kdata(:,i).*seq.decay;
    end

    if (seq.addnoise == 1)
        noiselevel = 0.005;
        noise = noiselevel*randn( size(seq.kdata)) + 1i*noiselevel*randn(size(seq.kdata));
    else 
        noise = 0;
    end

    seq.kdata = seq.kdata + noise;

    
    % nufft recon
    seq.psf = FT'*(seq.w.*seq.kdata);

    if (seq.B0correction == 1) && strcmp(seq.B0correction_method, 'cpr') == 1
        seq.psf2 = A'*(seq.w.*seq.kdata);
    elseif (seq.B0correction == 1) && strcmp(seq.B0correction_method,'icpr') == 1
        seq.psf2 = A'*(seq.w.*seq.kdata);
        R  = FiniteDifference(size(seq.offres_fmap));
        lambda  = 0.5;
        maxIter = 30;
        seq.psf2 = naive_qpwls_cg(A, seq.kdata, seq.psf2, R, lambda, maxIter);
    end

    % if strcmp(seq.traj_type, 'spaio') == 1
    %     seq.psf2 = pocs(seq,30);
    % end

    %[result.specificity, result.sharpening] = psf_2D_analysis(seq);

    %disp('TODO LIST...');
end
