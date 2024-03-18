function N_coh = Normalization_warren(Q, I, f_sqrd, f_av_sqrd, I_compton, rho, N_inc, spur_flag)

if strcmp(spur_flag, 'spur')
    int_1 = Q.^2.*(f_sqrd + I_compton + N_inc.*Q.^2);
else
    %int_1 = Q.^2.*(f_sqrd + I_compton);
    int_1 = Q.^2.*(f_sqrd + I_compton);
    %int_1 = Q.^2.*(f_sqrd);
end

int_2 = Q.^2.*I;

integral_1 = trapz(Q, int_1)  - 2*pi^2*rho;
integral_2 = trapz(Q, int_2);

N_coh = integral_1/integral_2;
end


% 
% 
%     "def Normalization_factor(Q, I, f_sqrd, f_av_sqrd, i_atom, electron_density, Z_average, convergence):\n",
%     "    from scipy import integrate\n",
%     "    # Normalization factor that makes f_sqrd and I comparable. Taken from: \n",
%     "    # Warren \n",
%     "    \n",
%     "    P_e = 2*np.pi**2*electron_density/Z_average\n",
%     "    \n",
%     "    # TEST\n",
%     "    #convergence = 1\n",
%     "    \n",
%     "    \n",
%     "    #int_1 = Q**2*convergence*(f_sqrd + i_atom - f_av_sqrd)/f_av_sqrd - P_e \n",
%     "    int_1 = Q**2*convergence*(f_sqrd + i_atom)/f_av_sqrd - P_e\n",
%     "    int_2 = Q**2*convergence*I/f_av_sqrd\n",
%     "    \n",
%     "    nom   = integrate.trapz(int_1, Q)\n",
%     "    den   = integrate.trapz(int_2, Q)\n",
%     "    \n",
%     "    return abs(nom/den)\n",
%     
%     %electron_density_tot.append(   particle_density_tot[-1]*average_Z_dummy)\n",
%     
%       "def Normalization_factor(Q, I, f_sqrd, i_atom, particle_density):\n",
%     "    int_1 = Q**2*(f_sqrd + i_atom) #/f_sqrd\n",
%     "    int_2 = I*Q**2 #/f_sqrd\n",
%     "    nom   = integrate.trapz(int_1, Q) - 2*np.pi**2*particle_density**2\n",
%     "    den   = integrate.trapz(int_2, Q)\n",
%     "    return nom/den\n",
%     "\n",