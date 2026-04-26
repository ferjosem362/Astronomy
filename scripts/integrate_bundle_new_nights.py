"""
Integra las 70 noches NUEVAS del bundle SERVAL en master_v3_v2.

Procedimiento:
1. Identifica noches del bundle no en master (within 0.5d)
2. Aplica Pourbaix2016 IDENTICO al de subtract_binary_pourbaix2016.py
3. Calcula offset zero-point usando las 45 noches que SI solapan
4. Aplica offset a las 70 nuevas
5. Centra residuos por sub_inst (igual que master)
6. Genera CSV master_v3_v3 = master_v3_v2 + 70 noches nuevas
"""
import os, sys
sys.path.insert(0, '.')

import numpy as np
import pandas as pd
from astropy.time import Time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ============================================================
# Parametros Pourbaix2016 - LITERAL del script original
# ============================================================
P_yr        = 79.91
P_s         = P_yr * 365.25 * 86400.0
P_days      = P_yr * 365.25
e_orb       = 0.524
i_deg       = 79.32
i_rad       = np.deg2rad(i_deg)
omega_B_deg = 232.3
omega_B     = np.deg2rad(omega_B_deg)
T_jyear     = 1955.66
V0_kms      = -22.390
V0_ms       = V0_kms * 1000.0
dVB_ms      = 329.0
M_A         = 1.133
M_B         = 0.972

T_obj = Time(T_jyear, format='jyear', scale='tdb')
T0_BJD = T_obj.jd

G = 6.67430e-11
M_sun = 1.98892e30

K_B = (2*np.pi*G * (M_A * M_sun)**3 / (P_s * ((M_A + M_B)*M_sun)**2))**(1/3) \
      * np.sin(i_rad) / np.sqrt(1 - e_orb**2)

print(f"[INFO] K_B = {K_B:.3f} m/s, T0_BJD = {T0_BJD:.4f}")

# ============================================================
# Funcion Pourbaix - LITERAL del script original
# ============================================================
def true_anomaly(t, P_days, e, T0_jd):
    M = 2*np.pi * (t - T0_jd) / P_days
    M = np.mod(M, 2*np.pi)
    E = M.copy()
    for _ in range(100):
        dE = (E - e*np.sin(E) - M) / (1 - e*np.cos(E))
        E -= dE
        if np.max(np.abs(dE)) < 1e-12:
            break
    nu = 2 * np.arctan2(np.sqrt(1+e)*np.sin(E/2), np.sqrt(1-e)*np.cos(E/2))
    return nu

def pourbaix_rv(bjd):
    """Modelo Pourbaix2016 RV de B en m/s para BJD dado."""
    nu = true_anomaly(bjd, P_days, e_orb, T0_BJD)
    return V0_ms + K_B * (e_orb*np.cos(omega_B) + np.cos(omega_B + nu)) + dVB_ms

# ============================================================
BUNDLE = "release/espresso_sflA_independent_SERVAL_bundle/data/espresso_sflA_drs_nightly_group_detrended_clean.csv"
MASTER = "data_interim/qc/master_v3_v2_nightly_residuals_clean.csv"
OUT_NEW   = "data_interim/qc/bundle_new_nights_processed.csv"
OUT_MASTER= "data_interim/qc/master_v3_v3_with_bundle.csv"
OUT_DIAG  = "data_interim/qc/bundle_integration_diagnostic.png"

print("\n" + "="*72)
print("INTEGRACION BUNDLE -> MASTER_V3_V3")
print("="*72)

bundle = pd.read_csv(BUNDLE)
master = pd.read_csv(MASTER)
print(f"\n  Bundle: N={len(bundle)}")
print(f"  Master: N={len(master)}")
print(f"  Master sub_inst:")
for s, n in master['sub_inst'].value_counts().items():
    print(f"    {s}: {n}")

# ============================================================
# Identificar nuevas noches (no en master, within 0.5 d)
# ============================================================
master_bjd_round = master.bjd.round(1).values
bundle['bjd_round'] = bundle.bjd.round(1)
bundle['in_master'] = bundle.bjd_round.isin(master_bjd_round)

new_nights = bundle[~bundle.in_master].copy()
overlap_nights = bundle[bundle.in_master].copy()

print(f"\n  Solape: {len(overlap_nights)} noches")
print(f"  Nuevas: {len(new_nights)} noches")

# ============================================================
# Aplicar Pourbaix2016 a TODO el bundle
# ============================================================
print(f"\n  Aplicando Pourbaix2016 al bundle...")
bundle['rv_pourbaix_ms'] = pourbaix_rv(bundle.bjd.values)
bundle['rv_resid_raw'] = bundle.rv_rel_ms - bundle.rv_pourbaix_ms

print(f"    Pourbaix range: {bundle.rv_pourbaix_ms.min():.1f} a {bundle.rv_pourbaix_ms.max():.1f} m/s")
print(f"    Residuos raw range: {bundle.rv_resid_raw.min():.1f} a {bundle.rv_resid_raw.max():.1f}")
print(f"    Std residuos raw: {bundle.rv_resid_raw.std():.2f}")

# ============================================================
# Calcular offset bundle vs master usando OVERLAP
# ============================================================
print(f"\n  Calculando offset bundle vs master usando {len(overlap_nights)} noches solapadas...")

overlap_nights['rv_pourbaix_ms'] = pourbaix_rv(overlap_nights.bjd.values)
overlap_nights['rv_resid_raw'] = overlap_nights.rv_rel_ms - overlap_nights.rv_pourbaix_ms

# Match con master
master_for_match = master[['bjd','rv_resid','sub_inst']].copy()
master_for_match['bjd_round'] = master_for_match.bjd.round(1)
overlap_nights['bjd_round'] = overlap_nights.bjd.round(1)

merged = overlap_nights.merge(
    master_for_match[['bjd_round','rv_resid','sub_inst']],
    on='bjd_round', how='inner', suffixes=('','_master')
)

if len(merged) > 0:
    diff = merged['rv_resid_raw'] - merged['rv_resid']
    offset = float(np.median(diff))
    offset_std = float(np.std(diff))
    print(f"    Mediana diff: {offset:+.3f} m/s")
    print(f"    Std diff:     {offset_std:.3f} m/s")
    print(f"    Range diff:   {diff.min():+.3f} / {diff.max():+.3f}")
    
    if abs(offset_std) > 5.0:
        print(f"    [WARN] Std de diferencias > 5 m/s - hay variacion ruidosa entre extracciones")
    
    if abs(offset) < 1.0:
        print(f"    [OK] Offset pequeño - usar zero offset")
        offset = 0.0
    elif abs(offset) < 10.0:
        print(f"    [OK] Offset moderado - aplicar al integrar")
    else:
        print(f"    [ALERT] Offset > 10 m/s - extracciones muy distintas")
else:
    print(f"    [ERROR] No hay solape para offset")
    sys.exit(1)

# ============================================================
# Procesar las NUEVAS noches con offset y centrar
# ============================================================
print(f"\n  Procesando {len(new_nights)} noches nuevas...")
new_nights['rv_pourbaix_ms'] = pourbaix_rv(new_nights.bjd.values)
new_nights['rv_resid_raw'] = new_nights.rv_rel_ms - new_nights.rv_pourbaix_ms
new_nights['rv_resid_offset'] = new_nights.rv_resid_raw - offset

# Centrar residuos a la mediana del nuevo sub-inst
median_new = float(new_nights.rv_resid_offset.median())
new_nights['rv_resid_centered'] = new_nights.rv_resid_offset - median_new

print(f"    Mediana centro: {median_new:+.3f} m/s")
print(f"    Std post-centrado: {new_nights.rv_resid_centered.std():.2f} m/s")

# Para comparar: que std tiene master ESPRESSO_EXT?
m_ext = master[master.sub_inst.str.contains('ESPRESSO_EXT|2023', case=False, na=False)]
if len(m_ext) > 0:
    centered_master = m_ext.rv_resid - m_ext.rv_resid.median()
    print(f"    Std master ESPRESSO_EXT centrado: {centered_master.std():.2f} m/s (referencia)")

# ============================================================
# Sanity check: si std_new > 2 * std_master_ext, hay problema
# ============================================================
if len(m_ext) > 0:
    std_new = new_nights.rv_resid_centered.std()
    std_ref = centered_master.std()
    if std_new > 2 * std_ref:
        print(f"\n  [ALERT] Std nuevas ({std_new:.2f}) > 2x master_ext ({std_ref:.2f})")
        print(f"  Las nuevas noches tienen ruido sistemático extra. Diagnosticar antes de integrar.")
    else:
        print(f"\n  [OK] Std nuevas consistente con master_ext")

# ============================================================
# Construir master_v3_v3
# ============================================================
print(f"\n  Construyendo master_v3_v3...")

# Preparar nuevas en formato master
new_to_add = pd.DataFrame({
    'bjd': new_nights.bjd.values,
    'rv': new_nights.rv_rel_ms.values,            # rv original
    'rv_err': new_nights.erv_ms.values,
    'sub_inst': 'ESPRESSO_EXT_BUNDLE',
    'rv_resid': new_nights.rv_resid_offset.values, # post-Pourbaix con offset
})

# Asegurar que las columnas comunes existen en master
required_cols = ['bjd', 'rv', 'rv_err', 'sub_inst', 'rv_resid']
master_subset = master[required_cols].copy()
master_v3_v3 = pd.concat([master_subset, new_to_add], ignore_index=True)
master_v3_v3 = master_v3_v3.sort_values('bjd').reset_index(drop=True)

master_v3_v3.to_csv(OUT_MASTER, index=False)
new_nights.to_csv(OUT_NEW, index=False)
print(f"    [OK] {OUT_MASTER}: N={len(master_v3_v3)}")
print(f"    [OK] {OUT_NEW}: N={len(new_nights)}")

# ============================================================
# Plot diagnostico
# ============================================================
fig, axes = plt.subplots(3, 1, figsize=(14, 10))

# Panel 1: rv_rel del bundle vs Pourbaix model en el rango bundle
t_dense = np.linspace(bundle.bjd.min()-30, bundle.bjd.max()+30, 2000)
rv_model_dense = pourbaix_rv(t_dense)

axes[0].scatter(bundle[~bundle.in_master].bjd, bundle[~bundle.in_master].rv_rel_ms,
                color='red', s=15, alpha=0.7, label=f'Bundle NEW ({(~bundle.in_master).sum()})')
axes[0].scatter(bundle[bundle.in_master].bjd, bundle[bundle.in_master].rv_rel_ms,
                color='blue', s=15, alpha=0.7, label=f'Bundle OVERLAP ({bundle.in_master.sum()})')
axes[0].plot(t_dense, rv_model_dense, 'k-', lw=1.2, alpha=0.7, label='Pourbaix2016')
axes[0].set_ylabel('rv_rel_ms')
axes[0].set_title('Bundle data vs Pourbaix2016 model')
axes[0].legend(fontsize=9)
axes[0].grid(alpha=0.3)

# Panel 2: residuos overlap (bundle vs master)
if len(merged) > 0:
    axes[1].scatter(merged.bjd, merged.rv_resid_raw, color='red', s=20, alpha=0.7,
                    label=f'Bundle resid (raw)')
    axes[1].scatter(merged.bjd, merged.rv_resid, color='blue', s=20, alpha=0.7,
                    label=f'Master resid')
    axes[1].axhline(offset, color='red', ls='--', alpha=0.5, label=f'Offset={offset:+.2f}')
    axes[1].axhline(0, color='blue', ls='--', alpha=0.5)
    axes[1].set_ylabel('rv_resid (m/s)')
    axes[1].set_title(f'Comparacion overlap: bundle vs master en {len(merged)} noches')
    axes[1].legend(fontsize=9)
    axes[1].grid(alpha=0.3)

# Panel 3: residuos finales master_v3_v3 (todos)
sub_colors = {'HARPS03':'tab:blue', 'HARPS15':'tab:cyan',
              'ESPRESSO19_HR11':'tab:purple', 'ESPRESSO_EXT_2023plus':'tab:orange',
              'ESPRESSO_EXT_BUNDLE':'red',
              'ES_zhao2017':'tab:green', 'CHIRON_zhao2017':'tab:olive',
              'ESPRESSO18_HR11':'tab:red', 'ESPRESSO18_HR21':'pink',
              'ESPRESSO18_UHR':'tab:brown'}
for sub in master_v3_v3['sub_inst'].unique():
    m = master_v3_v3['sub_inst'] == sub
    color = sub_colors.get(sub, 'gray')
    centered = master_v3_v3.loc[m, 'rv_resid'] - master_v3_v3.loc[m, 'rv_resid'].median()
    axes[2].scatter(master_v3_v3.loc[m, 'bjd'], centered,
                    color=color, s=8, alpha=0.5, label=f'{sub} (n={m.sum()})')
axes[2].axhline(0, color='k', lw=0.5)
axes[2].set_xlabel('BJD')
axes[2].set_ylabel('rv_resid centered (m/s)')
axes[2].set_title('Master v3_v3 residuos centrados por sub_inst')
axes[2].legend(fontsize=7, ncol=4, loc='lower left')
axes[2].set_ylim(-30, 30)
axes[2].grid(alpha=0.3)

fig.tight_layout()
fig.savefig(OUT_DIAG, dpi=120)
print(f"    [OK] {OUT_DIAG}")

# ============================================================
# Resumen final
# ============================================================
print(f"\n{'='*72}")
print(f"RESUMEN FINAL")
print(f"{'='*72}")
print(f"  master_v3_v2 (antes):  N={len(master)}")
print(f"  master_v3_v3 (despues): N={len(master_v3_v3)}  (+{len(new_nights)} nuevas)")
print(f"  ESPRESSO_EXT_BUNDLE std: {new_nights.rv_resid_centered.std():.2f} m/s")
print(f"  Offset aplicado:        {offset:+.3f} m/s")
print(f"  Output: {OUT_MASTER}")
print(f"\n  Si std_BUNDLE < 5 m/s, podemos lanzar MCMC V5 sobre master_v3_v3.")
