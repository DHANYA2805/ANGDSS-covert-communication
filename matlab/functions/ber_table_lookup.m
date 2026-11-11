function ber = ber_table_lookup(tbl, snr_lin)
%BER_TABLE_LOOKUP  Interpolate a precomputed BER-vs-chip-SNR table.
%   tbl.snr_db, tbl.ber (from ber_gdss_numeric or ber_dsss_theory).
%   Interpolation is linear in (dB, log10 BER); values outside the table are
%   clamped to its end points.
s = 10 * log10(max(snr_lin, 1e-30));
s = min(max(s, tbl.snr_db(1)), tbl.snr_db(end));
ber = 10.^interp1(tbl.snr_db, log10(max(tbl.ber, 1e-300)), s, 'linear');
end
