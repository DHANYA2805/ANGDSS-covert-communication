function [centers, h] = phase_hist4(r, nBins)
%PHASE_HIST4  Modulation-stripping detector for QPSK (base paper Sec. 4.2).
%   Raises the signal to the 4th power, histograms the phase over [-pi, pi)
%   (rotated by pi so that the QPSK spur, s^4 = -1, sits at 0 rad)
%   and normalises the peak to 1. A fixed QPSK constellation collapses onto
%   one phase (bell-shaped histogram); Gaussian-like signals stay flat.
edges   = linspace(-pi, pi, nBins + 1);
centers = (edges(1:end-1) + edges(2:end)) / 2;
ph = angle(-r(:).^4);                    % unit QPSK: s^4 = -1 -> spur moved to 0 rad
ph(ph >= pi) = -pi;                      % keep +pi in the first bin
c  = histc(ph, edges);                   % histc works in MATLAB and Octave
c  = c(1:nBins);
h  = c(:).' / max(c);
end
