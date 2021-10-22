function alpha = SOCA_CFAR_Alpha(PFA,RefWindow)
    syms a;

    summation = 0;
    for k = 0:1:(RefWindow/2-1)
        summation = summation + nchoosek((RefWindow/2 -1 +k),k)*(2+a)^(-k);
    end

    eqn = PFA == 2*(2+a)^(-RefWindow/2) * summation;

    alpha = double(vpasolve(eqn, a, [0 10^5]));
end