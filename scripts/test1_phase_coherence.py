"""
Test 1: ¿el pico a 270 d en ESPRESSO 2023-2024 tiene fase coherente
con master_v3 (19.6 yr)?

Procedimiento:
  A. Ajustar Kepleriano circular (P~271, K, phi0) en ESPRESSO solo
  B. Congelar P y phi0
  C. Predecir RV en master_v3 con esos parametros
  D. Calcular chi2_red en master_v3
  E. Comparar con master_v3 sin modelo (chi2 baseline)

Si chi2_red(model) << chi2_red(baseline) → fase coherente → Kepleriano probable.
Si chi2_red(model) ~ chi2_red(baseline) o peor → fase NO coherente → actividad.

Repetir para 232 d y 271 d (ambos sospechosos).
"""

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

INPUT_CCF    = "data_interim/qc/espresso_eso_ccf_curated_hr.csv"
INPUT_MASTER = "data_interim/qc/master_v3_v2_nightly_residuals_clean.csv"
OUTDIR       = "data_interim/qc"

# Periodos a testear
TARGETS = [
    {'P': 271.45, 'label': '271d (RV pol1 peak)'},
    {'P': 239.73, 'label': '240d (Halpha peak)'},
    {'P': 288.82, 'label': '289d (Halpha peak banda 270)'},
    {'P': 232.26, 'label': '232d (master_v3 GLS peak)'},
]

# ============================================================
# A. Cargar ESPRESSO 2023-2024 nightly, detrendear, ajustar Kepleriano
# ============================================================
print("="*72)
print("TEST 1: COHERENCIA DE FASE 1.4 yr ESPRESSO vs 19.6 yr MASTER_V3")
print("="*72)

df_ccf = pd.read_csv(INPUT_CCF)
df_ccf['bjd_int'] = df_ccf['time'].astype(int)
agg_ccf = df_ccf.groupby(['bjd_int','model_inst']).agg(
    bjd=('time','mean'),
    rv=('rv_centered','median'),
    rv_err=('rv_err', lambda x: np.sqrt(np.sum(x**2))/len(x)),
).reset_index()

# Detrend pol1 por sub-segmento (porque es donde el pico es mas fuerte)
detrended = []
for seg in agg_ccf['model_inst'].unique():
    m = agg_ccf['model_inst'] == seg
    t_seg = agg_ccf.loc[m, 'bjd'].values
    y_seg = agg_ccf.loc[m, 'rv'].values
    coef = np.polyfit(t_seg, y_seg, 1)
    agg_ccf.loc[m, 'rv_dt'] = y_seg - np.polyval(coef, t_seg)

agg_ccf = agg_ccf.sort_values('bjd').reset_index(drop=True)
t_e = agg_ccf['bjd'].values
y_e = agg_ccf['rv_dt'].values
dy_e = np.maximum(agg_ccf['rv_err'].values, 0.1)
print(f"\n  ESPRESSO nightly: N={len(agg_ccf)}, "
      f"baseline {(t_e.max()-t_e.min())/365.25:.2f} yr")

# Cargar master_v3 (sin las noches del externo, mantengamos ESPRESSO19_HR11 original)
df_v3 = pd.read_csv(INPUT_MASTER)
# Quitar las noches externas, dejar solo lo que era master_v2 + ESPRESSO19_HR11
EXT_SUBS = ['ESPRESSO_EXT_2023e', 'ESPRESSO_EXT_2023plus']
df_master_pre = df_v3[~df_v3['sub_inst'].isin(EXT_SUBS)].copy()
df_master_pre['rv_resid_centered'] = df_master_pre.groupby('sub_inst')['rv_resid'].transform(
    lambda x: x - x.median())

t_m = df_master_pre['bjd'].values
y_m = df_master_pre['rv_resid_centered'].values
dy_m = np.maximum(df_master_pre['rv_err'].values, 0.1)
print(f"  Master pre-merge: N={len(df_master_pre)}, "
      f"baseline {(t_m.max()-t_m.min())/365.25:.2f} yr")

# ============================================================
# Funcion modelo: Kepleriano circular
# ============================================================
def keplerian_circular(t, K, P, phi, offset):
    return K * np.sin(2*np.pi*t/P + phi) + offset

# ============================================================
# Loop sobre periodos target
# ============================================================
print(f"\n{'-'*72}")
print(f"{'Banda':<25} {'K_e (m/s)':>10} {'phi':>6} {'chi2_m_base':>12} "
      f"{'chi2_m_model':>13} {'mejora':>8}")
print(f"{'-'*72}")

results = []
fig, axes = plt.subplots(len(TARGETS), 2, figsize=(15, 3.5*len(TARGETS)))

for i, target in enumerate(TARGETS):
    P_target = target['P']
    label = target['label']
    
    # Ajuste sobre ESPRESSO solamente
    try:
        p0 = [2.0, P_target, 0.0, 0.0]
        bounds = ([0, P_target*0.99, -2*np.pi, -5],
                  [20, P_target*1.01, 2*np.pi, 5])
        popt_e, pcov_e = curve_fit(keplerian_circular, t_e, y_e, p0=p0,
                                    sigma=dy_e, absolute_sigma=True, bounds=bounds)
        K_e, P_e, phi_e, off_e = popt_e
        K_err = np.sqrt(pcov_e[0,0])
        
        # Chi2 baseline en master (sin modelo, solo std)
        chi2_m_base = np.sum((y_m / dy_m)**2) / (len(y_m) - 1)
        
        # Predecir master con modelo de ESPRESSO
        y_pred_master = keplerian_circular(t_m, K_e, P_e, phi_e, 0.0)  # offset 0, ya centrado
        # Ajustar solo offset master (libre porque puede haber zero point distinto)
        offset_m = np.median(y_m - y_pred_master)
        y_pred_master_offset = y_pred_master + offset_m
        
        chi2_m_model = np.sum(((y_m - y_pred_master_offset) / dy_m)**2) / (len(y_m) - 1)
        
        improvement = chi2_m_base / chi2_m_model
        
        verdict = ""
        if improvement > 1.5:
            verdict = "→ COHERENTE (mejora >1.5x)"
        elif improvement > 1.1:
            verdict = "→ marginal"
        elif improvement > 0.9:
            verdict = "→ no coherente"
        else:
            verdict = "→ EMPEORA (modelo de ESPRESSO no encaja en master)"
        
        print(f"{label:<25} {K_e:>9.3f}±{K_err:.2f} {phi_e:>+6.2f} "
              f"{chi2_m_base:>12.3f} {chi2_m_model:>13.3f} {improvement:>7.2f}x  {verdict}")
        
        results.append({
            'label': label, 'P': P_e, 'K_e': K_e, 'K_err': K_err,
            'phi_e': phi_e, 'chi2_m_base': chi2_m_base,
            'chi2_m_model': chi2_m_model, 'improvement': improvement,
            'verdict': verdict,
        })
        
        # Plot
        # Panel izquierdo: ESPRESSO con ajuste
        ax = axes[i, 0]
        t_fine = np.linspace(t_e.min(), t_e.max(), 1000)
        y_fit_e = keplerian_circular(t_fine, *popt_e)
        ax.errorbar(t_e, y_e, yerr=dy_e, fmt='.', ms=3, alpha=0.5, label='ESPRESSO dt')
        ax.plot(t_fine, y_fit_e, 'r-', lw=1.5, label=f'Fit: K={K_e:.2f} m/s')
        ax.set_ylabel('rv_dt (m/s)')
        ax.set_xlabel('BJD')
        ax.set_title(f'{label} - ESPRESSO fit')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        ax.set_ylim(-12, 12)
        
        # Panel derecho: master_v3 con prediccion
        ax = axes[i, 1]
        ax.errorbar(t_m, y_m, yerr=dy_m, fmt='.', ms=2, alpha=0.4, color='gray',
                    label=f'Master (N={len(t_m)})')
        # Phase folded master
        phase_m = ((t_m - t_e[0]) / P_e) % 1
        ax2 = ax.twinx()
        # No, mejor mostrar la prediccion temporal
        t_master_fine = np.linspace(t_m.min(), t_m.max(), 5000)
        y_pred_fine = keplerian_circular(t_master_fine, K_e, P_e, phi_e, offset_m)
        ax.plot(t_master_fine, y_pred_fine, 'r-', lw=0.5, alpha=0.7,
                label=f'Pred: chi2_red={chi2_m_model:.2f}')
        ax.axhline(0, color='k', lw=0.3)
        ax.set_xlabel('BJD')
        ax.set_ylabel('rv_resid (m/s)')
        ax.set_title(f'{label} - Master prediction (improvement {improvement:.2f}x)')
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
        ax.set_ylim(-15, 15)
        
    except Exception as e:
        print(f"{label}: error - {e}")

fig.suptitle('Test 1: coherencia de fase ESPRESSO 1.4yr → master_v3 19.6 yr', y=1.001)
fig.tight_layout()
out_plot = f"{OUTDIR}/test1_phase_coherence.png"
fig.savefig(out_plot, dpi=120)
print(f"\n[OK] {out_plot}")

df_res = pd.DataFrame(results)
df_res.to_csv(f"{OUTDIR}/test1_phase_coherence.csv", index=False)
print(f"[OK] {OUTDIR}/test1_phase_coherence.csv")

# ============================================================
# Veredicto final
# ============================================================
print(f"\n{'='*72}")
print("VEREDICTO TEST 1")
print('='*72)
print()
print("INTERPRETACION:")
print("  - improvement > 1.5x: la fase de ESPRESSO predice master → COHERENTE → Kepleriano")
print("  - improvement ~ 1.0:  fase incoherente → actividad de medio plazo")
print("  - improvement < 1.0:  modelo empeora → señal NO Kepleriana")
print()
if any(r['improvement'] > 1.5 for r in results):
    print("→ AL MENOS UN PERIODO PASA test 1. Pasar a test 2 (GP+Kepleriano).")
elif any(r['improvement'] > 1.1 for r in results):
    print("→ Algunos periodos marginales. Pasar a test 2 con cuidado.")
else:
    print("→ NINGUN PERIODO PASA test 1. La señal a 232-271 d NO es Kepleriano coherente")
    print("  Refuerza el descarte por chi2_red=602 que ya teniamos en master_v3.")
