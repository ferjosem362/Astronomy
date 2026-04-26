"""
Debug de K_B y T0 para Pourbaix2016.

Calcula K_B por dos vias independientes y verifica T0.
"""
import numpy as np
from astropy.time import Time
import astropy.units as u
import astropy.constants as const

# Parametros Pourbaix2016 Tabla 1
P_yr   = 79.91
e      = 0.524
i_deg  = 79.32
a_arcsec = 17.66        # separacion A-B en el plano del cielo
parallax_mas = 743      # paralaje
M_A    = 1.133
M_B    = 0.972
M_tot  = M_A + M_B
kappa  = 0.4617         # M_B/(M_A+M_B), publicado

# === Verificacion #1: kappa publicado vs calculado ===
kappa_calc = M_B / (M_A + M_B)
print(f"kappa publicado: {kappa}")
print(f"kappa calculado: {kappa_calc:.4f}")
print(f"  -> consistente: {abs(kappa-kappa_calc)<0.001}")

# === Calculo de a en AU desde paralaje ===
distance_pc = 1000.0 / parallax_mas
a_AU = a_arcsec * distance_pc   # parsec * arcsec = AU exactamente
print(f"\nDistancia: {distance_pc:.4f} pc")
print(f"a (semi-eje) total = {a_AU:.4f} AU")

# === Verificacion #2: a vs Kepler 3a ===
a_kepler = (P_yr**2 * M_tot)**(1/3)  # a en AU si P en yr y M en Msun
print(f"a por Kepler 3a:    {a_kepler:.4f} AU")
print(f"  -> deberia ser identico: diff = {abs(a_AU-a_kepler):.4f}")

# === Calculo K_B - via 1: formula directa de SB2 ===
# K_B = (2*pi*a_B*sin(i)) / (P*sqrt(1-e^2))
# donde a_B = a_total * (M_A / M_tot) = a_total * (1-kappa)
a_B_AU = a_AU * (M_A / M_tot)
print(f"\na_B (semieje de B alrededor del baricentro) = {a_B_AU:.4f} AU")

a_B_m = a_B_AU * const.au.to(u.m).value
P_s = P_yr * u.yr.to(u.s)
K_B_v1 = (2*np.pi*a_B_m * np.sin(np.deg2rad(i_deg))) / (P_s * np.sqrt(1-e**2))
print(f"K_B via #1 (a_B*sin(i)/P): {K_B_v1:.2f} m/s")

# === Calculo K_B - via 2: formula con masas ===
# K_B = (2*pi*G*M_A^3 / (P*M_tot^2))^(1/3) * sin(i)/sqrt(1-e^2)
G = const.G.value
M_sun = const.M_sun.value
K_B_v2 = (2*np.pi*G * (M_A*M_sun)**3 / (P_s*((M_tot)*M_sun)**2))**(1/3) \
         * np.sin(np.deg2rad(i_deg)) / np.sqrt(1-e**2)
print(f"K_B via #2 (formula masas):  {K_B_v2:.2f} m/s")

print(f"\nDiferencia entre v1 y v2: {abs(K_B_v1-K_B_v2):.2f} m/s")
print(f"  (deberian ser identicas si a viene de Kepler 3a)")

# === Conversion T0 ===
print("\n" + "="*60)
print("VERIFICACION T0 = 1955.66 jyear")
print("="*60)

# Opcion 1: jyear (Julian year)
t1 = Time(1955.66, format='jyear', scale='tdb')
print(f"jyear scale=tdb: BJD = {t1.jd:.4f}  -> calendario: {t1.iso}")

t1b = Time(1955.66, format='jyear', scale='utc')
print(f"jyear scale=utc: BJD = {t1b.jd:.4f}  -> calendario: {t1b.iso}")

# Opcion 2: byear (Besselian year)
t2 = Time(1955.66, format='byear', scale='tdb')
print(f"byear scale=tdb: BJD = {t2.jd:.4f}  -> calendario: {t2.iso}")

# Opcion 3: decimalyear
t3 = Time(1955.66, format='decimalyear', scale='tdb')
print(f"decimalyear:     BJD = {t3.jd:.4f}  -> calendario: {t3.iso}")

# Opcion 4: 1955 + 0.66*365.25 dias desde JD inicio del anio
print(f"\nDiferencias entre convenciones:")
print(f"  jyear vs byear:       {abs(t1.jd-t2.jd):.3f} dias")
print(f"  jyear vs decimalyear: {abs(t1.jd-t3.jd):.3f} dias")
