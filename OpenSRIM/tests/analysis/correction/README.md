# Untersuchung des `correction`-Integranden

Datum: 13.11.2026. Untersucht wurden die aktuellen lokalen Arbeitskopien.
Die Produktionsdateien wurden nicht verändert. Python-/Numba-Varianten wurden
im Speicher erstellt, Fortran-Varianten unter `/tmp/correction-study` kompiliert.

## Bedeutung

In tools/pytrim/cm_scatter.py gelten

```
b = sqrt(2-u²)
rho = r0/(e*p²)
g = sqrt(b² + rho*chi)
correction = rho*chi / (b*g*(g+b)) = 1/b - 1/g
```

Der Term ist exakt die Differenz zwischen freiem und wechselwirkendem
Winkelintegranden. Die rationalisierte Schreibweise vermeidet die Subtraktion
fast gleicher Zahlen bei schwacher Ablenkung. Er ist kein empirischer Fit.
Mit `A=2*acos(p/R)`, `B=2*acos(r0/R)` und
`Q=sum(w/g)/sum(w/b)` ist die endliche Winkelgleichung

```
theta = (A-B) + B*(1-Q) = A-B*Q.
```

`correction` berechnet den Anteil `1-Q`. Einfaches Nullsetzen lässt nur
`A-B` übrig und entfernt damit einen wesentlichen Teil des Streuwinkels.
Numba verwendet `h=p*g`; seine Formel ist algebraisch identisch. Fortran
verwendet in `f_theta_deflection_legendre` den vierfachen Integranden und
eine dazu passende Normierung. Untersucht wurde `LegendreChiCorr`.
Der zusätzliche Integrand wird nur für endliche Reichweite im Winkel verwendet:
ZBL mit unendlicher Reichweite sowie die Tau-Gleichung bleiben unverändert.

## Versuchsaufbau

- NLHlin, Elementpaare B–Si, H–H, Au–Au und H–Au, jeweils native Normierung.
- Reduzierte Energien `10**(-4), ..., 10**2` (7 Werte).
- Regulär: 201 äquidistante Stoßparameter von 0.01 bis 30, einschließlich
  der Nullstreuung oberhalb des jeweiligen Cutoffs.
- Rand: 14 Punkte `p=R*(1-10**(-k))`, k=2,...,15.
- Python und echtes Numba-JIT: 4 und 32 Gauss-Legendre-Punkte.
- Fortran: frisch kompilierte Varianten, 4 Gauss-Legendre-Punkte, gfortran -O0.
- Je Backend/Quadraturordnung 6020 Eingaben; gleiche Eingaben je Variante.
- Bezug für Abweichungen ist die aktuelle rationale Formel derselben
  Implementierung und Quadraturordnung, keine unabhängige exakte Lösung.
- Relative Abweichung: `abs(variante-aktuell)/abs(aktuell)`; Referenznullen
  ausgeschlossen. Absolute und relative Maxima können woanders auftreten.

## Reguläres Raster, 4 Quadraturpunkte

| Backend | Nullsetzen: max. absolut [rad] | Nullsetzen: max. relativ | Direkte Differenz 1/b−1/g: max. absolut [rad] | Direkte Differenz: max. relativ |
|---|---:|---:|---:|---:|
| Python | 3.024913736 | 99.9344865 % | 8.882e-16 | 6.580e-10 |
| Numba | 3.024913736 | 99.9344865 % | 8.882e-16 | 9.559e-10 |
| Fortran | 3.024913736 | 99.9344865 % | 8.882e-16 | 6.261e-10 |

Absolutes Maximum beim Nullsetzen: Au–Au, e=1, p=0.01:
aktuell theta=3.076701594 rad, ohne Beitrag 0.051787858 rad.
Relatives Maximum: Au–Au, e=100, p=0.01.
Für B–Si allein: maximal 2.933021590 rad bzw. 99.8200699 %.
Mit 32 Punkten bleiben die Größenordnungen gleich: maximal 3.024909127 rad
beim Nullsetzen und 1.333e-15 rad beim direkten Subtrahieren (Python/Numba).
Tau ist in allen endlichen Python-/Numba-Vergleichen identisch.

## Weitere Abgrenzungen

1. **Direkte Differenz im Integranden:** Auch im Randraster sehr kleine absolute
   Unterschiede (Python/Fortran kleiner als 6e-17 rad, Numba höchstens
   4.441e-16 rad bei 4 Punkten). Der robuste Term bleibt sinnvoll, weil die
   relative Kondition der direkten Subtraktion mit kleiner Ablenkung schlechter
   wird. Numba gibt `pi-theta` zurück; die Rückumrechnung verliert ihrerseits
   Auflösung für sehr kleine Winkel.
2. **Gesamte Winkelgleichung direkt als `A-B*Q`:** Python-Vergleich mit identischer
   Apsis und Quadratur: regulär höchstens 1.777e-14 rad; an den extremen
   Randpunkten bis 2.332e-9 rad. Maximale relative Randabweichung ca. 20.33 %
   bei einem aktuellen Winkel von nur 1.487e-10 rad. Dies umfasst auch den
   Verlust durch `acos(p/R)` und `acos(r0/R)`, nicht allein den Integranden.
3. **Historische Formel `A*(1-Q)`:** Dies ist eine andere Änderung als die
   Rationalisierung. Ihr fehlt die korrekte Trennung von p und r0.
   Im regulären Raster bis 0.343286021 rad bzw. 49.7856 % Abweichung
   (4 Punkte). Diese Formel sollte nicht als gleichwertiger Ersatz gelten.
4. **Vorhandenes Randproblem:** 18 Python-Auswertungen über beide Ordnungen
   liefern bereits im unveränderten Code NaN (9 je Ordnung). Die Apsis liegt
   dabei jenseits von R, sodass die Quadratwurzel für umax ungültig wird.
   Fortran liefert ebenfalls 9 nichtendliche Winkel pro Variante im Randraster;
   Numba liefert hier durchgehend endliche Ergebnisse.
   Alle Varianten des jeweiligen Backends sind dort betroffen; diese Werte sind aus den Maxima
   ausgeschlossen. Das ist unabhängig vom correction-Term. Die extremen
   Randzahlen sind daher keine Aussage über die physikalische Genauigkeit
   der Apsis oder der Gesamtlösung.

## Reproduzieren

Vom Repository-Stamm:

```
.venv/bin/python OpenSRIM/tests/analysis/correction/analyze.py
.venv/bin/python OpenSRIM/tests/analysis/correction/fortran.py
```

Benötigt die vorhandene Python-Umgebung und gfortran. Arbeitsdateien und
Rohdaten werden unter `/tmp/correction-study` abgelegt. `summary.txt` und
`fortran-summary.txt` enthalten die Ergebnisse dieses Laufs einschließlich
Koordinaten der Maxima. In `summary.txt` bedeutet Spalte 8 für Python/Numba
Nullsetzen, Spalte 9 direkte Integrandendifferenz; für `Formula` bedeutet
Spalte 8 `A-B*Q`, Spalte 9 die historische Formel `A*(1-Q)`.
