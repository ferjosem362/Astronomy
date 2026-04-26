"""
Diagnostico de los residuos HARPS03 post-Pourbaix2016.
"""
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

df = pd.read_csv('data_interim/qc/master_v2_nightly_residuals.csv')
m = df['sub_inst'] == 'HARPS03'
sub = df[m].copy()
sub['rv_resid_centered'] = sub['rv_resid'] - sub['rv_resid'].median()

print(f"N HARPS03: {len(sub)}")
print(f"\nResidual centered stats:")
print(f"  std:  {sub['rv_resid_centered'].std():.2f} m/s")
print(f"  MAD:  {(sub['rv_resid_centered'] - sub['rv_resid_centered'].median()).abs().median():.2f} m/s")
print(f"  std robusto (1.4826*MAD): {1.4826 * (sub['rv_resid_centered'] - sub['rv_resid_centered'].median()).abs().median():.2f} m/s")

# Distribucion de residuos
print(f"\nQuantiles de |resid_centered|:")
abs_resid = sub['rv_resid_centered'].abs()
for q in [50, 75, 90, 95, 99]:
    print(f"  {q}% -> {np.percentile(abs_resid, q):.2f} m/s")

# Cuantos puntos a >5 sigma robusto?
sigma_robust = 1.4826 * (sub['rv_resid_centered'] - sub['rv_resid_centered'].median()).abs().median()
n_5sig = (abs_resid > 5*sigma_robust).sum()
n_3sig = (abs_resid > 3*sigma_robust).sum()
print(f"\nPuntos con |resid| > 3*sigma_robust ({3*sigma_robust:.1f} m/s): {n_3sig}")
print(f"Puntos con |resid| > 5*sigma_robust ({5*sigma_robust:.1f} m/s): {n_5sig}")

# Std descontando outliers
clean = sub[abs_resid <= 5*sigma_robust]
print(f"\nStd tras quitar outliers >5*sigma_robust: {clean['rv_resid_centered'].std():.2f} m/s")
print(f"  N tras corte: {len(clean)} (perdimos {len(sub)-len(clean)})")

# Plot
fig, axes = plt.subplots(2, 1, figsize=(13, 7))

# Panel 1: serie temporal con outliers marcados
mask_in = abs_resid <= 5*sigma_robust
axes[0].errorbar(sub.loc[mask_in, 'bjd'], sub.loc[mask_in, 'rv_resid_centered'],
                 yerr=sub.loc[mask_in, 'rv_err'],
                 fmt='.', ms=3, alpha=0.6, color='tab:blue', label='Inliers')
axes[0].errorbar(sub.loc[~mask_in, 'bjd'], sub.loc[~mask_in, 'rv_resid_centered'],
                 yerr=sub.loc[~mask_in, 'rv_err'],
                 fmt='o', ms=6, alpha=0.8, color='red',
                 label=f'>5*sigma outliers ({(~mask_in).sum()})')
axes[0].axhline(0, color='k', lw=0.5)
axes[0].axhline(5*sigma_robust, color='gray', ls='--', lw=0.5)
axes[0].axhline(-5*sigma_robust, color='gray', ls='--', lw=0.5)
axes[0].set_ylabel('rv_resid centrado (m/s)')
axes[0].set_xlabel('BJD')
axes[0].legend()
axes[0].set_title(f'HARPS03 residuos post-Pourbaix2016 - std={sub["rv_resid_centered"].std():.1f}, std_clean={clean["rv_resid_centered"].std():.1f} m/s')
axes[0].grid(alpha=0.3)

# Panel 2: histograma
axes[1].hist(sub['rv_resid_centered'], bins=80, color='tab:blue', alpha=0.7, edgecolor='black')
axes[1].axvline(0, color='k', lw=0.5)
axes[1].axvline(5*sigma_robust, color='red', ls='--', lw=1, label=f'+/-5*sigma_robust = +/-{5*sigma_robust:.0f}')
axes[1].axvline(-5*sigma_robust, color='red', ls='--', lw=1)
axes[1].set_xlabel('rv_resid centrado (m/s)')
axes[1].set_ylabel('N')
axes[1].set_yscale('log')
axes[1].legend()
axes[1].set_title('Distribucion residuos HARPS03 (escala log)')
axes[1].grid(alpha=0.3)

fig.tight_layout()
out = 'data_interim/qc/harps03_residuals_diagnostic.png'
fig.savefig(out, dpi=120)
print(f"\n[OK] {out}")
