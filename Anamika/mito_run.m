for n = 1:1:steps
    if(m_run)
        p = rand;
        if(p>1-exp(-k_s*G(m(n))))
            m_stop = 1;
            m_run = 0;
            m(n+1) = m(n);
        else
            m(n+1) = m(n)+v;
        end
    else
        p = rand;
        if(p>(1-exp(-k_u)))
            m_stop = 0;
            m_run = 1;
            m(n+1) = m(n)+v;
        else
            m(n+1) = m(n);
        end
    end
end
